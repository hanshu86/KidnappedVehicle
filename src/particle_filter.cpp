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
#include <limits>

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
  num_particles = 10;  // Set the number of particles to 1000
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
    particles[i].weight = 1.0f;
  }

  is_initialized = true;
  std::cout << "Particle Vector Size is: " << particles.size() << std::endl;
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
    normal_distribution<double> dist_x(0,std_pos[0]);// for x
    normal_distribution<double> dist_y(0,std_pos[1]);// for y
    normal_distribution<double> dist_theta(0,std_pos[2]);// for theta

    double noise_x, noise_y, noise_theta;

    noise_x = dist_x(gen);
    noise_y = dist_y(gen);
    noise_theta = dist_theta(gen);

    particles[i].x = x_f + noise_x;
    particles[i].y = y_f + noise_y;
    particles[i].theta = theta_f + noise_theta;
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
  // run loop for each observation so that we can assign it to its nearest
  // landmark

  // for(int a = 0; a < predicted.size(); a++)
  // {
  //     std::cout << "predicted valid landmark id: " << predicted[a].id << std::endl;
  //     std::cout << "predicted valid landmark x: " << predicted[a].x << std::endl;
  //     std::cout << "predicted valid landmark y: " << predicted[a].y << std::endl;
  // }

  for (int i = 0; i < observations.size(); i++) {
    double min_distance = std::numeric_limits<double>::max();
    int idx = -1;

    for (int p = 0; p < predicted.size(); p++) {
      double curr_distance = dist(predicted[p].x, predicted[p].y, observations[i].x, observations[i].y);
      if (curr_distance < min_distance){
          min_distance = curr_distance;
          idx = predicted[p].id;
      }
    } // predicted landmark loop

    observations[i].id = idx; // assign the closest landmark id to observation
    // std::cout << "Observation " << i << " is associated with landmark  " << idx << std::endl;
    // std::cout << "Observation X " << observations[i].x << std::endl;
    // std::cout << "Observation Y " << observations[i].y << std::endl;
    // for(int a = 0; a < predicted.size(); a++) {
    //   if(predicted[a].id == idx) {
    //     std::cout << "Landmark X " << predicted[a].x << std::endl;
    //     std::cout << "Landmark Y " << predicted[a].y << std::endl;
    //   }
    // }

  } // observation loop

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

  double weight_normalizer = 0.0;
  for (int i = 0; i < num_particles; i++) {
    
    vector<LandmarkObs> map_transformed_observations;

    for (int j = 0; j < observations.size(); j++) {
        LandmarkObs obs;
        std::vector<double> map_coordinate = transformedObservation(particles[i], observations[j].x, observations[j].y);
        obs.id = j;
        obs.x = map_coordinate[0];
        obs.y = map_coordinate[1];
        map_transformed_observations.push_back(obs);
    } // observation loop

    // take only landmarks in the sensor_range for data association step
    vector<LandmarkObs> valid_landmarks;

    for (int k = 0; k < map_landmarks.landmark_list.size(); k++) {
        
        Map::single_landmark_s curr_landmark = map_landmarks.landmark_list[k];
        LandmarkObs with_in_range_landmark;

        if( (fabs(curr_landmark.x_f - particles[i].x) <= sensor_range) && (fabs(curr_landmark.y_f - particles[i].y) <= sensor_range)) {
            with_in_range_landmark.id = curr_landmark.id_i;
            with_in_range_landmark.x = curr_landmark.x_f;
            with_in_range_landmark.y = curr_landmark.y_f;
            valid_landmarks.push_back(with_in_range_landmark);
            // std::cout << "valid landmark id: " << with_in_range_landmark.id << std::endl;
            // std::cout << "valid landmark x: " << with_in_range_landmark.x << std::endl;
            // std::cout << "valid landmark y: " << with_in_range_landmark.y << std::endl;
          }
    } // valid landmark loop
    // std::cout << "Valid landmark Size: " << valid_landmarks.size() << std::endl;

    dataAssociation(valid_landmarks, map_transformed_observations);
    particles[i].weight = 1.0f;

    for (int l = 0; l < map_transformed_observations.size(); l++) {
      double x_obs = map_transformed_observations[l].x;
      double y_obs = map_transformed_observations[l].y;
      int nearest_landmark = map_transformed_observations[l].id;
      double mu_x = 0.0f;
      double mu_y = 0.0f;
      double obsWeight = 0.0f;
      bool found = false;
      for (int b = 0; b < valid_landmarks.size(); b++) {
        if(valid_landmarks[b].id == nearest_landmark) {
          mu_x = valid_landmarks[b].x;
          mu_y = valid_landmarks[b].y;
          found = true;
          break;
        }//if condition to find matching id with observation 
      } // landmark loop to find matching id for observation
      if (found == false)
      {
        std::cout << "Could not find landmark..assert" << std::endl;
        while(true); 
      }
      obsWeight = multiv_prob(std_landmark[0], std_landmark[1], x_obs, y_obs, mu_x, mu_y);
      particles[i].weight *= obsWeight;
      // std::cout << "X_obs: " << x_obs << std::endl;
      // std::cout << "Y_obs: "<< y_obs << std::endl;
      // std::cout << "mu_x: " << mu_x << std::endl;
      // std::cout << "mu_y: " << mu_y  << std::endl;
      // std::cout << "nearest landmark " << nearest_landmark << std::endl;
      // std::cout << "weight_obs: "<< obsWeight << std::endl;
    } // weight update loop
    // std::cout << "particle weight: " << particles[i].weight << std::endl;
    weight_normalizer += particles[i].weight;
  } // each particle loop

  // std::cout << "total weight sum: " << weight_normalizer << std::endl;
  double newSum = 0;
  weights.clear();
  //normalize particle weight so that they will be in range 0 to 1
  for (int i = 0; i < num_particles; i++) {
    // particles[i].weight /= weight_normalizer;
    newSum += particles[i].weight;
    weights.push_back(particles[i].weight);
  } //weight normalizer loop
  // std::cout << "total new weight sum: " << newSum << std::endl;
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  //This code is mostly converted from Sebastian Python code (from class) to C++
  // Vector for new particles
  vector<Particle> new_particles (num_particles);
  std::default_random_engine gen;

  // generate random starting index for resampling wheel
  std::uniform_int_distribution<int> unidist(0, num_particles-1);
  auto index = unidist(gen);
  // std::cout << "Uniformly distributed index: " << index << std::endl;

  // get max weight
  double max_weight = *max_element(weights.begin(), weights.end());
  // std::cout << "max weight : " << max_weight << std::endl;

  std::uniform_real_distribution<double> unirealdist(0.0, max_weight);

  double beta = 0.0;

  // spin the resample wheel!
  for (int i = 0; i < num_particles; i++) {
    beta += unirealdist(gen) * 2.0;
    // std::cout << "New Beta: " << beta << std::endl;
    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    // std::cout << "choosen index: " << index << std::endl;
    new_particles.push_back(Particle());
    new_particles[i].id = particles[index].id;
    new_particles[i].x = particles[index].x;
    new_particles[i].y = particles[index].y;
    new_particles[i].theta = particles[index].theta;
    new_particles[i].weight = particles[index].weight;
    new_particles[i].associations = particles[index].associations;
    new_particles[i].sense_x = particles[index].sense_x;
    new_particles[i].sense_y = particles[index].sense_y;
    // std::cout << "X = :" << particles[index].x << "Y = : " << particles[index].y << std::endl;
    // std::cout << "X->:" << new_particles[i].x << "Y->: " << new_particles[i].y << std::endl;
  }

  //assign resampled particle to particles
  particles = new_particles;

  // print particles to see what all particles has been resamples
  // for (int i = 0; i < num_particles; i++) {
  //   std::cout << "X: " << particles[i].x << "Y: " << particles[i].y << std::endl;
  // }
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
  double x_new = 0.0f;

  if (fabs(yaw_rate) < 0.0001){
    x_new = particle.x + (velocity * cos(particle.theta) * delta_t);
  }
  else {
    x_new = particle.x + ((velocity / yaw_rate) * (sin(particle.theta + (yaw_rate * delta_t)) - sin(particle.theta)));
  }

  return x_new;
}

double ParticleFilter::getPredictedY(Particle particle, double delta_t,  double velocity, double yaw_rate) {
  double y_new = 0.0f;

  if (fabs(yaw_rate) < 0.0001){
    y_new = particle.y + (velocity * sin(particle.theta) * delta_t);
  }
  else {
    y_new = particle.y + ((velocity / yaw_rate) * (cos(particle.theta) - cos(particle.theta + (yaw_rate * delta_t))));
  }

  return y_new;
}

double ParticleFilter::getPredictedTheta(Particle particle, double delta_t,  double velocity, double yaw_rate) {
  double theta_new = 0.0f;

  if (fabs(yaw_rate) < 0.0001){
    theta_new = particle.theta;
  }
  else {
    theta_new = particle.theta + (yaw_rate * delta_t);
  }

  return theta_new;
}

std::vector<double> ParticleFilter::transformedObservation(Particle particle, double x_obs, double y_obs) {
  std::vector<double> map_coordinate_vector;

  double x_map, y_map;
  x_map = particle.x + (cos(particle.theta) * x_obs) - (sin(particle.theta) * y_obs);
  y_map = particle.y + (sin(particle.theta) * x_obs) + (cos(particle.theta) * y_obs);

  map_coordinate_vector.push_back(x_map);
  map_coordinate_vector.push_back(y_map);

  return map_coordinate_vector;  
}

//This function is taken from UDACITY Particle weight solution
double ParticleFilter::multiv_prob(double sig_x, double sig_y, double x_obs, double y_obs,
                   double mu_x, double mu_y) {
  // calculate normalization term
  double gauss_norm;
  gauss_norm = 1 / (2 * M_PI * sig_x * sig_y);

  // calculate exponent
  double exponent;
  exponent = (pow(x_obs - mu_x, 2) / (2 * pow(sig_x, 2)))
               + (pow(y_obs - mu_y, 2) / (2 * pow(sig_y, 2)));
    
  // calculate weight using normalization terms and exponent
  double weight;
  weight = gauss_norm * exp(-exponent);
    
  return weight;
}
