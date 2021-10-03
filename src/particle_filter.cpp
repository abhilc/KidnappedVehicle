/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>
#include <float.h>

#include "helper_functions.h"

using std::string;
using std::vector;
using std::normal_distribution;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  
  num_particles = 750;  // TODO: Set the number of particles
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  weights = std::vector<double>(static_cast<unsigned long>(num_particles), 1.0);
  
  double sample_x, sample_y, sample_theta;
  sample_x = dist_x(gen);
  sample_y = dist_y(gen);
  sample_theta = dist_theta(gen);
  

  for(int i=0;i<num_particles;i++)
  {
    //Create a Particle P
    Particle P;
    P.x = sample_x;
    P.y = sample_y;
    P.theta = sample_theta;
    P.weight = 1.0;
    P.id = i;
    
    //Push to particles vector
    particles.push_back(P);
  }
  is_initialized = true; 
  
}



void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  double x, y, theta, x_pred, y_pred, theta_pred;
  double v = velocity, ydot = yaw_rate;
  
  //Noise components to account for prediction that are later added
  std::default_random_engine gen;
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);

  
  for(int i=0;i<num_particles;i++)
  {
      x = particles[i].x;
      y = particles[i].y;
      theta = particles[i].theta;
      theta_pred = 0;
    
      if(fabs(ydot) < 0.00001) //Check to avoid divide by zero
      {
        //Predict without yaw_rate
        x_pred = x + velocity * delta_t * cos(theta);
        y_pred = y + velocity * delta_t * sin(theta);
      }
      else
      {
        //Predict with yaw_rate
        double v_div_yawrate = velocity/yaw_rate;
        double delta_t_times_yaw_rate = yaw_rate * delta_t;
        x_pred = x + v_div_yawrate * (sin(theta + delta_t_times_yaw_rate) - sin(theta));
        y_pred = y + v_div_yawrate * (cos(theta) - cos(theta + delta_t_times_yaw_rate));
        theta_pred = theta + delta_t_times_yaw_rate;
      }
      
      particles[i].x = x_pred + dist_x(gen);
      particles[i].y = y_pred + dist_y(gen);
      particles[i].theta = theta_pred + dist_theta(gen);
  }
  

}


void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  for(auto &observation: observations)
  {
    double min_dist = 1e18;
    
    observation.id = -1;
    
    for(auto &pred_meas: predicted)
    {
      double cur_dist = dist(pred_meas.x, pred_meas.y, observation.x, observation.y);
      if(min_dist > cur_dist)
      {
        min_dist = cur_dist;
        observation.id = pred_meas.id;
      }
    }
  }

}




void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  
  /** Steps in this algorithm are as follows 
  * 1. First transform each observation from vehicle co-ordinate system to 
  *    co-ordinate system of the particle under consideration
  * 2. Create a vector transformed observations that contains a list of all  
  *    results from point 1
  * 3. Create a vector of landmark positions from the map data
  *
  * 4. Using the two vectors in points 2 and 3, call dataAssociation function.
  *    This function associates each transformed observation, a nearest landmark
  *    from the map 
  * 5. Use Multivariate Guassian Probability Density function to finally 
  *.   update the weights of the particles
  */  
  
  
  
  //Part I. Translation of co-ordinate system(Vehicle to Map)
  
  
  for(int j=0;j<particles.size();j++)
  {
    Particle const &p = particles[j];
    
    
    /**
      * Convert the Map instance to LandmarkObs vector that are only within the sensor range
      * Here, the 'predicted' is the vector of possible landmarks detected by each particle
      * The pos of each particle is considered to be the pos of the car
      * So for each particle, predicted measurements is computed below
      **/
      

    vector<LandmarkObs> landmarks;
    for (auto const &landmark : map_landmarks.landmark_list) {
      if (dist(p.x, p.y, landmark.x_f, landmark.y_f) <= sensor_range) {
        LandmarkObs Obs_lm = {
          .id = landmark.id_i,
          .x = (double)(landmark.x_f),
          .y = (double)(landmark.y_f),
        };
        landmarks.push_back(Obs_lm);
      }
    }
    
    /* Observations:
     * Observations are sensor data. Eg: We detect a landmark 1m ahead
     * These observations are from the perspective of the real pos 
     * of the car. Hence these need to be converted from car co-ordinate
     * system to the particle's coordinate system because in this case,
     * we are considering particle's pos as car's pos
     *
    */
    
    vector<LandmarkObs> transformed_observations(observations.size());
    
    for (auto i = 0; i < observations.size(); i++) {
      
      LandmarkObs observation = observations[i];
      double cos_theta = cos(p.theta);
      double sin_theta = sin(p.theta);
      transformed_observations[i].x = p.x + cos_theta * observation.x - sin_theta * observation.y;
      transformed_observations[i].y = p.y + sin_theta * observation.x + cos_theta * observation.y;
      transformed_observations[i].id = -1;
    }

    //Data Association
    dataAssociation(landmarks, transformed_observations);
    
    vector<double> observation_probabilities(transformed_observations.size());
    particles[j].weight = 1.0;  
    
    for (auto i = 0; i < observations.size(); i++) {
      LandmarkObs tobs = transformed_observations[i];
      LandmarkObs nearest_landmark = {
          .id = -1,  
          .x = static_cast<double>(map_landmarks.landmark_list[tobs.id - 1].x_f), 
          .y = static_cast<double>(map_landmarks.landmark_list[tobs.id - 1].y_f),
      };

      double x_diff_2 = pow(tobs.x - nearest_landmark.x, 2.0);
      double y_diff_2 = pow(tobs.y - nearest_landmark.y, 2.0);
      double std_x_2 = pow(std_landmark[0], 2.0);
      double std_y_2 = pow(std_landmark[1], 2.0);

      
      observation_probabilities[i] = (1 / (2 * M_PI * std_landmark[0] * std_landmark[1])) *
                                     exp(-(x_diff_2 / (2 * std_x_2) + y_diff_2 / (2 * std_y_2)));

      
      particles[j].weight *= observation_probabilities[i];
    }
    weights[j] = particles[j].weight;
  }
  
}



void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  std::default_random_engine gen;
  std::discrete_distribution<size_t> dist_index(weights.begin(), weights.end());

  vector<Particle> resampled_particles(particles.size());

  for (auto i = 0; i < particles.size(); i++) {
    resampled_particles[i] = particles[dist_index(gen)];
  }

  particles = resampled_particles;
  
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}