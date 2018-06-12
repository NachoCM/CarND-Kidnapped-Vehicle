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
	num_particles=100;

	particles = std::vector<Particle>(num_particles);
	// normal distributions for x, y and theta
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);


	for (auto& p:particles){
		p.x=dist_x(gen);
		p.y=dist_y(gen);
		p.theta=dist_theta(gen);
		p.weight=1;
	}
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	for (auto& p:particles){
		
		if (yaw_rate==0){
			//TODO Check if cos and sin are not inverted
			p.x=p.x+velocity*delta_t*cos(p.theta);
			p.y=p.y+velocity*delta_t*sin(p.theta);
		} else {
			p.x=p.x+(velocity/yaw_rate)*(sin(p.theta+yaw_rate*delta_t)-sin(p.theta));
			p.y=p.y+(velocity/yaw_rate)*(cos(p.theta)-cos(p.theta+yaw_rate*delta_t));
		}
		
		p.theta=p.theta+yaw_rate*delta_t;

		p.x=normal_distribution<double>(p.x, std_pos[0])(gen);
		p.y=normal_distribution<double>(p.y, std_pos[1])(gen);
		p.theta=normal_distribution<double>(p.theta, std_pos[2])(gen);
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
	for (auto& p:particles){
		std::vector<LandmarkObs> transformed_obs;
		//Initialize current particle weight
		p.weight=1;

		//Transform observations to map coordinate system
		for (auto& car_obs:observations){
			LandmarkObs map_obs;
			map_obs.x=p.x + cos(p.theta)*car_obs.x - sin(p.theta)*car_obs.y;
			map_obs.y=p.y + sin(p.theta)*car_obs.x + cos(p.theta)*car_obs.y;
			transformed_obs.push_back(map_obs);
		}

		for (auto& o:transformed_obs){
			//Initialize landmark associated to current observation and distance to 
			//nearest neighbour
			o.id=0;
			double min_dist=std::numeric_limits<double>::max();

			//Traverse map landmarks (filter using sensor range would imply additional computation)
			for (int i = 0;i < map_landmarks.landmark_list.size();i++){
				Map::single_landmark_s l=map_landmarks.landmark_list[i];
				double distance=dist(o.x,o.y,l.x_f,l.y_f);
				if (distance<min_dist){
					o.id=i;
					min_dist=distance;
				}
			}

			//Extract landmark coordinates
			double landmark_x=map_landmarks.landmark_list[o.id].x_f;
			double landmark_y=map_landmarks.landmark_list[o.id].y_f;

			//Calculate weight for current particle
			long double gauss_norm = 1/(2*M_PI*std_landmark[0]*std_landmark[1]);
			long double exponent = pow(o.x - landmark_x,2)/(2*pow(std_landmark[0],2)) + pow(o.y - landmark_y,2)/(2*pow(std_landmark[1],2));
			long double weight = gauss_norm * exp(-exponent);
			if (weight>0){
				p.weight*=weight;
			}

			p.associations.push_back(o.id+1);
			p.sense_x.push_back(o.x);
			p.sense_y.push_back(o.y);
		}
	}
	

}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	discrete_distribution<int> distribution(weights.begin(), weights.end());
	vector<Particle> resampled;
	for (auto& p:particles){
		resampled.push_back(particles[distribution(gen)]);
	}

	particles=resampled;

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
