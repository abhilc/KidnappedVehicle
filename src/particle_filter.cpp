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
  num_particles = 1000;  // TODO: Set the number of particles
  std::default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  
  double sample_x, sample_y, sample_theta;
  sample_x = dist_x(gen);
  sample_y = dist_y(gen);
  sample_theta = dist_theta(gen);
  
  for(int i=0;i<num_particles;i++)
  {
    particles[i].x = sample_x;
    particles[i].y = sample_y;
    particles[i].theta = sample_theta;
    particles[i].weight = 1.0;
  }
  
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
  std::default_random_engine gen;
  double x, y, theta, x_pred, y_pred, theta_pred;
  double v = velocity, ydot = yaw_rate;

  
  for(int i=0;i<num_particles;i++)
  {
      x = particles[i].x;
      y = particles[i].y;
      theta = particles[i].theta;
      
      x_pred = x + v/ydot * (sin(theta + ydot * delta_t) - sin(theta));
      y_pred = y + v/ydot * (cos(theta) - cos(theta + ydot * delta_t))
      theta_pred = theta + ydot * delta_t;
      
      normal_distribution<double> dist_x(x_pred, std_pos[0]);
      normal_distribution<double> dist_y(y_pred, std_pos[1]);
      normal_distribution<double> dist_theta(theta_pred, std_pos[2]);
      
      particles[i].x = dist_x(gen);
      particles[i].y = dist_y(gen);
      particles[i].theta = dist_theta(gen);
  }

}

void ParticleFilter::distance(double x1, double y1, double x2, double y2)
{
  return sqrt(((x2 - x1)^2 + (y2 - y1)^2));
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
  
  for(int i=0;i<predicted.size();i++)
  {
    double dist = DBL_MAX;
    LandmarkObs best_observation = observations[0];
    for(int j=0;j<observations.size();j++)
    {
      double cur_dist = distance(predicted[i].x, predicted[i].y, observations[j].x, observations[j].y);
      if(dist > cur_dist)
      {
        dist = cur_dist;
        best_observation = observations[j];
      }
    }
    predicted[i] = best_observation;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  double x_obs, y_obs, x_part, y_part;
  double x_trans, y_trans;
  vector<LandmarkObs> predicted;
   
  //First transform the observations
  for(int i=0;i<observations.size();i++)
  {
    x_obs = observations[i].x;
    y_obs = observations[i].y;
    
    for(int j=0;j<num_particles;j++)
    {
      x_part = particles[i].x;
      y_part = particles[i].y;
      
      //Translated and Rotated co-ordinates in (x_trans, y_trans)
      rotate(x_part, y_part, x_obs, y_obs, &x_trans, &y_trans);
      LandmarkObs obs = {x_trans, y_trans};
      predicted.push_back(obs);
    }
  }
  
  //Data Association
  

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

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