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

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // Set up Particle vector
  num_particles = 1000;
  particles = std::vector<Particle>();
  weights = std::vector<double>(); 

  // Create normal (Gaussian) distributions for x, y, theta to account for noise
  default_random_engine gen;
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
      std::cout << "ParticleFilter::init| Particle " << i + 1 << " " << particles[i].x << " " << particles[i].y << " " 
                << particles[i].theta << std::endl;
    }
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  // Create normal (Gaussian) distributions for velocity, yaw_rate to account for noise
  default_random_engine gen;
  normal_distribution<double> dist_v_x(velocity, std_pos[0]);
  normal_distribution<double> dist_v_y(velocity, std_pos[1]);
  normal_distribution<double> dist_yaw(yaw_rate, std_pos[1]);

  // Apply the prediction step to each particle
  for (int i = 0; i < num_particles; ++i) {
    // Obtain noisey values for velocity, yaw_rate
    double velocity_x_n = dist_v_x(gen);
    double velocity_y_n = dist_v_y(gen);
    double yaw_rate_n = dist_yaw(gen);

    // Update yaw and positions
    double theta_f = particles[i].theta + (yaw_rate_n*delta_t);
    double x_f = particles[i].x + (velocity_x_n/yaw_rate_n) * (sin(theta_f) - sin(particles[i].theta));
    double y_f = particles[i].y + (velocity_y_n/yaw_rate_n) * (cos(particles[i].theta) - cos(theta_f));

    // Set the updated Particle values
    particles[i].x = x_f;
    particles[i].y = y_f;
    particles[i].theta = theta_f;

    // Print the first 10 particles -- set this to `true` when debugging.
    if (false && i < 10) {
      std::cout << "ParticleFilter::prediction| Particle " << i + 1 << " " << x_f << " " << y_f << " " << theta_f << std::endl;
    }
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