/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */
#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
  //   x, y, theta and their uncertainties from GPS) and all weights to 1. 
  // Add random Gaussian noise to each particle.
  // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    default_random_engine gen;
    // initial number of perticles
    num_particles = 100;

    // normal distribution for x,y, theta
    normal_distribution<double> N_x(x, std[0]);
    normal_distribution<double> N_y(y, std[1]);
    normal_distribution<double> N_theta(theta, std[2]);

    // generating particles
    for (int i = 0; i < num_particles; i++) {
      Particle particle;
      particle.id = i;
      particle.x = N_x(gen);
      particle.y = N_y(gen);
      particle.theta = N_theta(gen);
      particle.weight = 1.0;

      particles.push_back(particle);
      weights.push_back(1);
    }

    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
  // TODO: Add measurements to each particle and add random Gaussian noise.
  // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
  //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
  //  http://www.cplusplus.com/reference/random/default_random_engine/

  // random number generation
  default_random_engine gen;
  double n_x;
  double n_y;
  double n_theta;

  for (int i = 0; i < num_particles; i++) {
    if (yaw_rate == 0.0) {  
      n_x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
      n_y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
      n_theta = particles[i].theta;
    } 
    else {
      n_x = particles[i].x + velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
      n_y = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
      n_theta =  particles[i].theta + yaw_rate * delta_t;
    }

    // Adding noise
    normal_distribution<double> N_x(n_x, std_pos[0]);
    normal_distribution<double> N_y(n_y, std_pos[1]);
    normal_distribution<double> N_theta(n_theta, std_pos[2]);

    particles[i].x = N_x(gen);
    particles[i].y = N_y(gen);
    particles[i].theta = N_theta(gen);
  }

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
  // TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
  //   observed measurement to this particular landmark.
  // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
  //   implement this method and use it as a helper during the updateWeights phase.
  double min_dist;
  double distance;

  for(int i = 0; i < observations.size(); i++) {
    min_dist = 99999;
    for(int j = 0; j < predicted.size(); j++) {
      distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
      if (min_dist > distance) {
        min_dist = distance;
        observations[i].id = predicted[j].id;
      }
    }
  }
    
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
    std::vector<LandmarkObs> observations, Map map_landmarks) {
  // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
  //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
  // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
  //   according to the MAP'S coordinate system. You will need to transform between the two systems.
  //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
  //   The following is a good resource for the theory:
  //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
  //   and the following is a good resource for the actual equation to implement (look at equation 
  //   3.33
  //   http://planning.cs.uiuc.edu/node99.html
  const double sigma_x = std_landmark[0];
  const double sigma_y = std_landmark[1];
  for (int i = 0; i < num_particles; ++i){
    vector<LandmarkObs> trans_observations;
    LandmarkObs obs;

    // transformed observations
    for(int j = 0; j < observations.size(); j++){
      LandmarkObs trans_obs;
      obs = observations[j];
      // coord transformation from vehicle coord to map coord
      trans_obs.x = cos(particles[i].theta) * obs.x - sin(particles[i].theta) * obs.y + particles[i].x;
      trans_obs.y = sin(particles[i].theta) * obs.x + cos(particles[i].theta) * obs.y + particles[i].y;
      trans_observations.push_back(trans_obs);
    }
    
    // filtering predictions
    vector<LandmarkObs> predicted_landmarks;
    for(int j = 0; j < map_landmarks.landmark_list.size(); j++){
      double distance = dist(particles[i].x, particles[i].y, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f);
      if(distance <= sensor_range){
        LandmarkObs landmark = {map_landmarks.landmark_list[j].id_i, map_landmarks.landmark_list[j].x_f, map_landmarks.landmark_list[j].y_f};
        predicted_landmarks.push_back(landmark);
      }
    }
        
    dataAssociation(predicted_landmarks,trans_observations);

    // updating weights
    particles[i].weight = 1.0;
    for (int m = 0; m < trans_observations.size(); m++){
      LandmarkObs assoc;
      LandmarkObs obs = trans_observations[m];
      
      for (int n=0 ; n < predicted_landmarks.size(); n++)
      {
        if (trans_observations[m].id == predicted_landmarks[n].id){
          assoc = predicted_landmarks[n];
      } 

    }
      double x_diff = pow(obs.x - assoc.x, 2) /  (2 * pow(sigma_x, 2));
      double y_diff = pow(obs.y - assoc.y, 2) / (2 * pow(sigma_y, 2));
      particles[i].weight *= 1 / (2 * M_PI * sigma_x * sigma_y) * exp(-(x_diff + y_diff));
  }
    weights[i] = particles[i].weight;

  }             
}

void ParticleFilter::resample() {
  // TODO: Resample particles with replacement with probability proportional to their weight. 
  // NOTE: You may find std::discrete_distribution helpful here.
  //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
  default_random_engine gen;
  discrete_distribution<int> distribution(weights.begin(), weights.end());

  vector<Particle> resample_particles;

    for(int i = 0; i < num_particles; i++){
      resample_particles.push_back(particles[distribution(gen)]);
    }
    particles = resample_particles;

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
  //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates

  //Clear the previous associations
  particle.associations.clear();
  particle.sense_x.clear();
  particle.sense_y.clear();

  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;

  return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
  vector<int> v = best.associations;
  stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
  vector<double> v = best.sense_x;
  stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
  vector<double> v = best.sense_y;
  stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}