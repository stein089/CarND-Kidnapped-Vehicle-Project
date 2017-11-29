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

 // std::cout << "Init " << std::endl;

  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);

  num_particles = 10;

  // resize vectors
  weights.resize(num_particles);
  particles.resize(num_particles);

  default_random_engine gen;

  for (int i = 0; i < num_particles; ++i) {
    particles[i].id = i;
    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
    particles[i].weight = 1;
    weights[i] = particles[i].weight;
  }

  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http:www.cplusplus.com/reference/random/default_random_engine/

 // std::cout << "Prediction " << std::endl;

  default_random_engine gen;

  for (int i = 0; i < num_particles; ++i) {
    if (fabs(yaw_rate) > 0.001) {
      particles[i].x = particles[i].x + velocity / yaw_rate * (sin(particles[i].theta  + yaw_rate * delta_t) - sin(particles[i].theta ));
      particles[i].y = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta ) - cos(particles[i].theta  + yaw_rate * delta_t));
    } else {
      particles[i].x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
      particles[i].y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
    }

    particles[i].theta = particles[i].theta + yaw_rate * delta_t;

    normal_distribution<double> dist_x(particles[i].x, std_pos[0]);
    normal_distribution<double> dist_y(particles[i].y, std_pos[1]);
    normal_distribution<double> dist_theta(particles[i].theta, std_pos[2]);

    particles[i].x = dist_x(gen);
    particles[i].y = dist_y(gen);
    particles[i].theta = dist_theta(gen);
  }


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
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


 // std::cout << "UpdateWeights " << std::endl;

  //std::cout << "map_landmarks.landmark_list.size(): " << map_landmarks.landmark_list.size() << std::endl;
  //std::cout << "observations.size(): " << observations.size() << std::endl;
  //for (int j = 0; j < observations.size(); j++) {
  //    std::cout << "observations "<< j << ": " << observations[j].id << std::endl;
  //}

  for (int i = 0; i < particles.size(); i++) {

    //std::cout << "particle "<< i << "ID=" << particles[i].id << "("<< particles[i].x << "," << particles[i].y << ")" << std::endl;

    particles[i].associations.clear();
    particles[i].sense_x.clear();
    particles[i].sense_y.clear();

    for (int j = 0; j < observations.size(); j++) {

      double x_map_obs = particles[i].x + (cos(particles[i].theta) * observations[j].x) - (sin(particles[i].theta) * observations[j].y);
      double y_map_obs = particles[i].y + (sin(particles[i].theta) * observations[j].x) + (cos(particles[i].theta) * observations[j].y);

      double dist = -1;
      int NN_index = 0;

      for (int k = 0; k <= map_landmarks.landmark_list.size(); k++)
      {
        double dx = x_map_obs - map_landmarks.landmark_list[k].x_f;
        double dy = y_map_obs - map_landmarks.landmark_list[k].y_f;

        double new_dist = sqrt(dx*dx + dy*dy);

        if (new_dist < dist || dist < 0)
        {
          dist = new_dist;
          NN_index = k;
        }
      }
      double x_map = map_landmarks.landmark_list[NN_index].x_f;
      double y_map = map_landmarks.landmark_list[NN_index].y_f;

      particles[i].associations.push_back(map_landmarks.landmark_list[NN_index].id_i);
      particles[i].sense_x.push_back(x_map_obs);
      particles[i].sense_y.push_back(y_map_obs);


    //  for(int idx = 0; idx< particles[i].associations.size(); idx++){
   //     std::cout << "ass " << idx << " " << particles[i].associations[idx] << std::endl;
   //   }


      //std::cout << "observations[j].x: " << observations[j].x << " observations[j].y: " << observations[j].y << std::endl;
      //std::cout << "x_map: " << x_map << " y_map: " << y_map << std::endl;
      //std::cout << "x_map_obs: " << x_map_obs << " y_map_obs: " << y_map_obs << std::endl;
      //std::cout << "dist: " << dist << std::endl;


      double gauss_norm = (1.0/(2.0 * M_PI * std_landmark[0] * std_landmark[1]));
      double exponent = ((x_map_obs- x_map)*(x_map_obs- x_map))/(2.0 * std_landmark[0]*std_landmark[0]) +
              ((y_map_obs - y_map)*(y_map_obs - y_map))/(2.0 * std_landmark[1]*std_landmark[1]);

      particles[i].weight = particles[i].weight * gauss_norm * exp(-exponent);
    }
  }

  double weight_sum = 0.0;

  for (int i = 0; i < particles.size(); i++) {
    weight_sum += particles[i].weight;
  }

  for (int i = 0; i < particles.size(); i++) {
    particles[i].weight = particles[i].weight/weight_sum;
    weights[i] = particles[i].weight;
  }



}


//Resample particles with replacement with probability proportional to their weight.
void ParticleFilter::resample() {

  default_random_engine gen;
  std::discrete_distribution<int> weights_dist (weights.begin(), weights.end());

  std::vector<Particle> new_particles;
  new_particles.resize(num_particles);

 // for (int i = 0; i < weights.size(); ++i) {
 //   std::cout << "weight "<< i << ": " << weights[i] << std::endl;
 // }

  for (int i = 0; i < num_particles; ++i) {
    int sampled_index = weights_dist(gen);
  //  std::cout << "new index" << i << "; sampled_index: " << sampled_index << std::endl;
    new_particles[i] = particles[sampled_index];
  }


  particles = new_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
