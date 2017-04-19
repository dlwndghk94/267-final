#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h> 
#include <time.h>
#include <cmath>
#include <string.h>
#include "robot.cpp"
#include "common.cpp"


// float bearing_noise = 0.1;
// float steering_noise = 0.1;
// float distance_noise = 5.0;

void get_position(Robot* p, int N, float* rtn){
	float x = 0.0;
	float y = 0.0;
	float orientation = 0.0;
	for (int i = 0; i < N; i++){
		x = x + p[i].x;
		y = y + p[i].y;
		orientation = orientation + fmod(p[i].orientation - p[0].orientation + M_PI, 2.0 * M_PI) + p[0].orientation - M_PI;
	}
	rtn[0] = x/N;
	rtn[1] = y/N;
	rtn[2] = orientation/N;
	return;
}

// float *output = particle_filter(motion, measurement, num_particles, particles);

float *particle_filter(float motion[2], float measurement[4], int N, Robot * p, FILE* fp, bool output){
	
	// Motion update (prediction)
	// Robot p2[N];

	for (int i = 0; i<N; i++){
		p[i] = p[i].move(motion);
	}

	// Measurement update
	float w[N];
	float mw = -999999999.0;

	for (int i = 0; i < N; i++){
		w[i] = p[i].measurement_prob(measurement);
		if (mw < w[i]){
			mw = w[i];
		}
	}

	// Resampling
	Robot p3[N];
	int index = rand() % N;
	float beta = 0.0;
	for (int i = 0; i < N; i++){
		float rand_num = (double)rand() / (double)RAND_MAX;
		beta = beta + rand_num * 2.0 * mw;

		while( beta > w[index]){
			beta = beta - w[index];
			index = (index +1) % N;
		}
		p3[i] = p[index];
	}

	for (int i = 0; i < N; i++) {
		p[i] = p3[i];
		if (output) {
			fprintf(fp,"%f, %f, ",p[i].x,p[i].y);
		}	
	}

	if (output){
		fprintf(fp,"\n");
	}

	float *rtn = (float *) malloc(3* sizeof(float));
	get_position(p, N, rtn);
	return rtn;
}


int main(){

	int num_motions = 1000;
	int num_particles = 5000;
	float length = 20.0;
	char* filename = "output.csv";
	FILE* fp = fopen(filename, "a+"); 

	double simulation_time = read_timer();
	// Initializing particles array
	Robot particles[num_particles];
	for (int i = 0; i < num_particles; i++){
		Robot r;
		r.initialize(length);
		r.set_noise(bearing_noise, steering_noise, distance_noise);
		particles[i] = r;
	}


	Robot car;
	car.initialize(length);
	car.set(20.0, 20.0, 0);
	car.set_noise(bearing_noise, steering_noise, distance_noise);

	float measurement[4];
	float motion[2];

	// clock_t start = clock(), diff;
	for (int i = 0; i < num_motions; i++) {
		
		motion[0] = 2.0 *M_PI / 10.0;
		motion[1] = 20.0;
		car = car.move(motion);
		// comment back in
		fprintf(fp,"%f, %f, ",car.x,car.y);
		car.sense(measurement, 1);
		float *output = particle_filter(motion, measurement, num_particles, particles, fp, false);
	}
	simulation_time = read_timer() - simulation_time;
	printf("Time taken %f seconds\n", simulation_time);
	fclose(fp);
}