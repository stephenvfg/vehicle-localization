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

// Define a random engine generator to use throughout
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Set up Particle vector
  num_particles = 1000;
  particles = std::vector<Particle>();
  weights = std::vector<double>(); 

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
    weights.push_back(p.weight);

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
  for (auto & part : particles) {
    // Obtain old values for x, y, theta
    double x = part.x;
    double y = part.y;
    double theta = part.theta;

    // Prepare to calculate predictions based on whether or not the vehicle is turning
    double theta_f, x_f, y_f;

    if (yaw_rate < 0.001) {
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
    part.x = dist_x(gen);
    part.y = dist_y(gen);
    part.theta = dist_theta(gen);

    // Print a subset of particles -- set this to `true` when debugging.
    if (false && part.id < 20) {
      std::cout << "ParticleFilter::prediction| Particle " << part.id + 1 << " " << part.x 
                << " " << part.y << " " << part.theta << std::endl;
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
  // Create normal (Gaussian) distributions for landmark measurements to account for noise
  //normal_distribution<double> dist_x(velocity, std_pos[0]);
  //normal_distribution<double> dist_y(velocity, std_pos[1]);

  // Apply the prediction step to each particle
  //for (int i = 0; i < num_particles; ++i) {

  //}
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