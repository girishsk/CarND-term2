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
	num_particles = 100;
    std::default_random_engine generator;
    std::normal_distribution<double> dist_x(x,std[0]);
    std::normal_distribution<double> dist_y(y,std[1]);
    std::normal_distribution<double> dist_theta(theta,std[2]);
	for(int i = 0; i < num_particles; i++){
		Particle p = Particle();

		p.x = dist_x(generator);
		p.y = dist_y(generator);
		p.theta = dist_theta(generator);
		p.id = i;
		p.weight = 1.0;
		particles.push_back(p);
		weights.push_back(1.);
	}
    is_initialized = true;



}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/


    std::default_random_engine generator;
    for(int i = 0; i < num_particles; i++){
        Particle p = particles[i];
        double x;
        double y;
        double theta;
        if(fabs(yaw_rate) < 0.0001){
            x = p.x + velocity*delta_t*cos(p.theta);
            y = p.y + velocity*delta_t*sin(p.theta);
            theta = p.theta;
        }
        else
        {
            x = p.x + velocity/yaw_rate*(sin(p.theta+yaw_rate*delta_t) - sin(p.theta));
            y = p.y + velocity/yaw_rate*(cos(p.theta) - cos(p.theta+yaw_rate*delta_t));
            theta = p.theta + yaw_rate*delta_t;

        }
        std::normal_distribution<double> dist_x(x,std_pos[0]);
        std::normal_distribution<double> dist_y(y,std_pos[1]);
        std::normal_distribution<double> dist_theta(theta,std_pos[2]);

        particles[i].x = dist_x(generator);
        particles[i].y = dist_y(generator);
        particles[i].theta = dist_theta(generator);


    }


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.


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
    std::default_random_engine generator;
    for(int i =0; i < num_particles; i++){
        Particle p = particles[i];
        double dist = 0.0;
        double prob = 1.0;
        vector<LandmarkObs> trans_observations;

        for(int j=0; j < observations.size(); j++) {
            LandmarkObs obs = observations[j];
            LandmarkObs trans_ob;


            trans_ob.x = particles[i].x + (obs.x * cos(particles[i].theta) - obs.y * sin(particles[i].theta));
            trans_ob.y = particles[i].y + (obs.x * sin(particles[i].theta) + obs.y * cos(particles[i].theta));
            trans_ob.id = obs.id;

            trans_observations.push_back(trans_ob);
        }

        particles[i].weight = 1.;

        // find all the landmarks that are in the range of the particle
        std::vector<LandmarkObs> predicted;
        for(auto landmark:map_landmarks.landmark_list){
            double dist = sqrt((landmark.x_f-p.x)*(landmark.x_f-p.x) +(landmark.y_f-p.y) *(landmark.y_f-p.y) );
            if(dist <= sensor_range)   {
                LandmarkObs landmarkObs;
                landmarkObs.x = landmark.x_f;
                landmarkObs.y = landmark.y_f;
                landmarkObs.id = landmark.id_i;
                predicted.push_back(landmarkObs) ;

            }

            
        }

        LandmarkObs nearest;


        for (auto obs: trans_observations){
            bool found = false;

            double shortest = 1E10;
            for(auto pred: predicted){
                double dist = sqrt((pred.x-obs.x)*(pred.x-obs.x) +(pred.y-obs.y) *(pred.y-obs.y));
                if(dist < shortest){
                    shortest = dist;
                    nearest = pred;
                    found = true;

                }
            }
            if(found) {
                double dx = obs.x - nearest.x;
                double dy = obs.y - nearest.y;

                prob *= 1.0 / (2 * M_PI * std_landmark[1] * std_landmark[0]) *
                        exp(-dx * dx / (2 * std_landmark[0] * std_landmark[0])) *
                        exp(-dy * dy / (2 * std_landmark[1] * std_landmark[1]));
            }



        }
        particles[i].weight = prob;
        weights[i] = prob;

        


    }
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    std::default_random_engine generator;
    std::discrete_distribution<int> distribution(weights.begin(),weights.end());

    vector<Particle> resampled;
    for(int i = 0; i < num_particles; i++){
        Particle p = particles[distribution(generator)];
        resampled.push_back(p);
    }
    particles = resampled;

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
