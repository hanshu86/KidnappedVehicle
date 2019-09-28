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
  num_particles = 1000;  // Set the number of particles to 1000
  std::default_random_engine gen;

  //Create Gaussian (normal) distribution
  normal_distribution<double> dist_x(x,std[0]);// for x
  normal_distribution<double> dist_y(y,std[1]);// for y
  normal_distribution<double> dist_theta(theta,std[2]);// for theta

  //crete particles and add it to vector
  for (int i = 0; i < num_particles; i++) {

    double sample_x, sample_y, sample_theta;

    sample_x = dist_x(gen);
    sample_y = dist_y(gen);
    sample_theta = dist_theta(gen);

    //Push back new particle created with default constructor
    particles.push_back(Particle());
    //modify it with correct value
    particles[i].id = i;
    particles[i].x = sample_x;
    particles[i].y = sample_y;
    particles[i].theta = sample_theta;
    particles[i].weight = 1;
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
  std::default_random_engine gen;

  for (int i = 0; i < num_particles; i++) {
    //Add measurement
    double x_f = getPredictedX(particles[i], delta_t, velocity, yaw_rate);
    double y_f = getPredictedY(particles[i], delta_t, velocity,  yaw_rate);
    double theta_f = getPredictedTheta(particles[i], delta_t, velocity, yaw_rate);
    
    //Add random gaussian noise
    normal_distribution<double> dist_x(x_f,std[0]);// for x
    normal_distribution<double> dist_y(y_f,std[1]);// for y
    normal_distribution<double> dist_theta(theta_f,std[2]);// for theta

    double sample_x, sample_y, sample_theta;

    noise_x = dist_x(gen);
    noise_y = dist_y(gen);
    noise_theta = dist_theta(gen);

    particles[i].x = noise_x;
    particles[i].y = noise_y;
    particles[i].theta = noise_theta;
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

double ParticleFilter::getPredictedX(Particle particle, double delta_t,  double velocity, double yaw_rate) {
  double x_new = 0;

  x_new = particle.x + (velocity / yaw_rate) * (sin(particle.theta + (yaw_rate * delta_t)) - sin(particle.theta));

  return x_new;
}

double ParticleFilter::getPredictedY(Particle particle, double delta_t,  double velocity, double yaw_rate) {
  double y_new = 0;

  y_new = particle.y + (velocity / yaw_rate) * (cos(particle.theta) - cos(particle.theta + (yaw_rate * delta_t)));

  return y_new;
}

double ParticleFilter::getPredictedTheta(Particle particle, double delta_t,  double velocity, double yaw_rate) {
  double theta_new = 0;

  theta_new = particle.theta + (yaw_rate * delta_t);

  return theta_new;
}

