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
using std::default_random_engine;
using std::normal_distribution;
using std::discrete_distribution;
using std::max_element;

// Define a random engine generator to use throughout
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Set up Particle vector
  num_particles = 100;

  // Create normal (Gaussian) distributions for x, y, theta to account for noise
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  // Initialize each particle with approximations from the GPS location data
  for (int i = 0; i < num_particles; ++i) { 
    // Initialize particle
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1;

    // Add Particle and weight to respective vectors
    particles.push_back(p);

    // Print the first 10 particles -- set this to `true` when debugging.
    if (false && i < 10) {
      std::cout << "ParticleFilter::init| Particle " << i + 1 << " " << particles[i].x << " " 
                << particles[i].y << " " << particles[i].theta << std::endl;
    }
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  // Apply the prediction step to each particle
  for (int i = 0; i < num_particles; ++i) {
    // Obtain old values for x, y, theta
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;

    // Prepare to calculate predictions based on whether or not the vehicle is turning
    double theta_f, x_f, y_f;

    if (abs(yaw_rate) < 0.00001) {
      // Update yaw and positions using the equations for a vehicle going straight
      x_f = x + velocity * delta_t * cos(theta);
      y_f = y + velocity * delta_t * sin(theta);
      theta_f = theta;
    } else {
      // Update yaw and positions using the equations for a turning vehicle
      theta_f = theta + (yaw_rate*delta_t);
      x_f = x + (velocity/yaw_rate) * (sin(theta_f) - sin(theta));
      y_f = y + (velocity/yaw_rate) * (cos(theta) - cos(theta_f));
    }

    // Create normal (Gaussian) distributions for velocity, yaw_rate to account for noise
    normal_distribution<double> dist_x(x_f, std_pos[0]);
    normal_distribution<double> dist_y(y_f, std_pos[1]);
    normal_distribution<double> dist_theta(theta_f, std_pos[2]);

    // Set the updated Particle values based on the noisy distribution
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);

    // Print a subset of particles -- set this to `true` when debugging.
    if (false && particles[i].id < 20) {
      std::cout << "ParticleFilter::prediction| Particle " << particles[i].id + 1 << " " 
                << particles[i].x << " " << particles[i].y << " " << particles[i].theta << std::endl;
    }
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  // Cycle through each observation to find its closest prediction
  for (auto & obs : observations) {

    // Define a placeholder to represent the minimum distance between predicted<->observed so far
    double min_dist = std::numeric_limits<double>::max();

    // Calculate the distance for every prediction to find the nearest neighbor
    for (auto & pred : predicted) {
      double curr_dist = dist(obs.x, obs.y, pred.x, pred.y);
      if (curr_dist < min_dist) {
        min_dist = curr_dist;
        obs.id = pred.id;
      }
    }
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  // Cycle through each particle
  for (int i = 0; i < num_particles; ++i) {

    // Isolate the landmarks that are in range of the sensor
    vector<LandmarkObs> map_landmarks_ir = vector<LandmarkObs>();
    for (int j = 0; j < map_landmarks.landmark_list.size(); ++j) {
      double distance = dist(particles[i].x, particles[i].y, 
                        map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);
      if (distance <= sensor_range) {
        LandmarkObs landmark_ir = LandmarkObs();
        landmark_ir.x = map_landmarks.landmark_list[j].x_f;
        landmark_ir.y = map_landmarks.landmark_list[j].y_f;
        landmark_ir.id = map_landmarks.landmark_list[j].id_i;
        map_landmarks_ir.push_back(landmark_ir);
      }
    }

    // Transform the vehicle observations into landmark coordinates for each observation
    vector<LandmarkObs> trans_observations = vector<LandmarkObs>();
    for (int j = 0; j < observations.size(); ++j) {
      // Transform sensor landmark observations to map coordinates
      LandmarkObs trans_obs = LandmarkObs(); 
      trans_obs.x = particles[i].x + cos(particles[i].theta)*observations[j].x 
                    - sin(particles[i].theta)*observations[j].y;
      trans_obs.y = particles[i].y + sin(particles[i].theta)*observations[j].x 
                    + cos(particles[i].theta)*observations[j].y;
      trans_obs.id = observations[j].id;
      trans_observations.push_back(trans_obs);
    }

    // Associate the observations to the nearest landmark on the map
    dataAssociation(map_landmarks_ir, trans_observations);

    // Prepare to calculate the particle weight
    double particle_weight = 1.0;

    // Cycle through the associated observations to update the respective particle weights
    for (int j = 0; j < trans_observations.size(); ++j) {
      double mu_x, mu_y;

      for (int k = 0; k < map_landmarks_ir.size(); ++k) {
        // Find the observation <-> landmark associated pair
        if (trans_observations[j].id == map_landmarks_ir[k].id) {
          mu_x = map_landmarks_ir[k].x;
          mu_y = map_landmarks_ir[k].y;
          break;
        }
      }

      // Update weights using the Gaussian probability density function
      double gauss_norm = 1 / (2*M_PI * std_landmark[0] * std_landmark[1]);
      double exponent = (pow(trans_observations[j].x - mu_x, 2) / (2 * pow(std_landmark[0], 2)))
                        + (pow(trans_observations[j].y - mu_y, 2) / (2 * pow(std_landmark[1], 2)));
      double obs_weight = gauss_norm * exp(-1*exponent);
      if (obs_weight < 0.00001) {
        obs_weight = 0.00001;
      }
      particle_weight *= obs_weight;
    }

    // Set the particle weight
    particles[i].weight = particle_weight;
  }

  // Normalize the particle weights
  double normalize_factor = 0.0;
  for (int i = 0; i < num_particles; i++) {
    normalize_factor += particles[i].weight;
  }

  for (int i = 0; i < num_particles; i++) {
    // Prevent division by zero
    particles[i].weight /= normalize_factor + std::numeric_limits<double>::epsilon();
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */

  // Retrieve weights for simpler access; store max weight
  weights = vector<double>();
  double max_weight = 0.0;

  for (auto & part : particles) {
    weights.push_back(part.weight);
    if (part.weight > max_weight) {
      max_weight = part.weight;
    }
  }

  // Set up the resampled particle vector and random number generator
  vector<Particle> resampled_particles = vector<Particle>();
  discrete_distribution<> rand_num_particles(0, num_particles);

  int index = rand_num_particles(gen);
  double beta = 0.0;

  for (int i = 0; i < num_particles; ++i) {
    beta += 2.0 * ((double) rand()/(RAND_MAX)) * max_weight;
    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    resampled_particles.push_back(particles[index]);
  }

  particles = resampled_particles;

  // Print the resampled particles -- set this to `true` when debugging.
  if (false) {
    for (int i = 0; i < num_particles; ++i) {
      std::cout << "ParticleFilter::resample| Particle " << i + 1 << " " << particles[i].x << " " 
                << particles[i].y << " " << particles[i].theta << " " << particles[i].weight << std::endl;
    }
  }
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