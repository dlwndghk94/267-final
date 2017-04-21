#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h> 
#include <time.h>
#include <cmath>
#include <string.h>
#include "robot.cpp"
#include "common.cpp"
#include "omp.h"


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


// float *particle_filter(float motion[2], float measurement[4], int N, Robot * p, FILE* fp, bool output){

// 	#pragma omp parallel for
// 	for (int i = 0; i<N; i++){
// 		p[i] = p[i].move(motion);
// 	}

// 	// Measurement update
// 	float w[N];
// 	float mw = -999999999.0;

// 	#pragma omp parallel for reduction(max:mw)
// 	for (int i = 0; i < N; i++){
// 		w[i] = p[i].measurement_prob(measurement);
// 		if (mw < w[i]){
// 			mw = w[i];
// 		}
// 	}

// 	// Resampling
// 	Robot p3[N];
// 	#pragma omp parallel for
// 	for (int i = 0; i < N; i++){
// 		int index = rand() % N;
// 		float beta = 0.0;
// 		float rand_num = (double)rand() / (double)RAND_MAX;
// 		beta = beta + rand_num * 2.0 * mw;

// 		while( beta > w[index]){
// 			beta = beta - w[index];
// 			index = (index +1) % N;
// 		}
// 		p3[i] = p[index];
// 	}

// 	#pragma omp parallel for
// 	for (int i = 0; i < N; i++) {
// 		p[i] = p3[i];
// 	}

// 	if (output) {
// 		for (int i = 0; i < N; i++) { 
// 			fprintf(fp,"%f, %f, ",p[i].x,p[i].y);
// 		}
// 		fprintf(fp,"\n");
// 	}

// 	float *rtn = (float *) malloc(3* sizeof(float));
// 	get_position(p, N, rtn);
// 	return rtn;
// }


int main(){
	srand(1);
	int num_motions = 1;
	int num_particles = 100000;
	float length = 20.0;
	char* filename = (char*) "output.csv";
	FILE* fp = fopen(filename, "a+");

	double simulation_time = read_timer();
	// Initializing particles array
	Robot p[num_particles];
	#pragma omp parallel for
	for (int i = 0; i < num_particles; i++){
		Robot r;
		r.initialize(length);
		r.set_noise(bearing_noise, steering_noise, distance_noise);
		p[i] = r;
	}

	Robot car;
	car.initialize(length);
	car.set(20.0, 20.0, 0);
	car.set_noise(bearing_noise, steering_noise, distance_noise);


	float measurement[4];
	float motion[2];
	int N = num_particles;
	Robot p3[N];
	float w[N];
	float mw;

	float movement_time = 0;
	float measurement_time = 0;
	float resampling_time = 0;
	float reassignment_time = 0;

	#pragma omp parallel 
	{
		for (int t = 0; t < num_motions; t++) {
			
			#pragma omp master 
			{
				movement_time -= read_timer();
			}

			#pragma omp master
			{
				motion[0] = 2.0 *M_PI / 10.0;
				motion[1] = 20.0;
				car = car.move(motion);
				//fprintf(fp,"%f, %f, ",car.x,car.y);
				car.sense(measurement, 1);
				float mw = -999999999.0;
			}
			#pragma omp barrier

			// -----------------------------------//
			// 	  PARTICLE FILTER STARTS HERE 	  //
			// -----------------------------------//

			// move all particles
			#pragma omp for
			for (int i = 0; i<N; i++){
				p[i] = p[i].move(motion);
			}
			#pragma omp barrier

			#pragma omp master 
			{
				movement_time += read_timer();
			}

			#pragma omp master 
			{
				measurement_time -= read_timer();
			}
			// Measurement update
			#pragma omp for reduction(max:mw)
			for (int i = 0; i < N; i++){
				w[i] = p[i].measurement_prob(measurement);
				if (mw < w[i]){
					mw = w[i];
				}
			}
			#pragma omp barrier

			#pragma omp master 
			{
				measurement_time += read_timer();
			}

			#pragma omp master 
			{
				resampling_time -= read_timer();
			}
			// Resampling
			#pragma omp for
			for (int i = 0; i < N; i++){
				int index = rand() % N;
				float beta = 0.0;
				float rand_num = (double)rand() / (double)RAND_MAX;
				beta = beta + rand_num * 2.0 * mw;

				while( beta > w[index]){
					beta = beta - w[index];
					index = (index +1) % N;
				}
				p3[i] = p[index];
			}
			#pragma omp barrier
			#pragma omp master 
			{
				resampling_time += read_timer();
			}


			#pragma omp master 
			{
				reassignment_time -= read_timer();
			}
			#pragma omp for
			for (int i = 0; i < N; i++) {
				p[i] = p3[i];
			}
			#pragma omp barrier
			#pragma omp master 
			{
				reassignment_time += read_timer();
			}

			// if (output) {
			// 	for (int i = 0; i < N; i++) { 
			// 		fprintf(fp,"%f, %f, ",p[i].x,p[i].y);
			// 	}
			// 	fprintf(fp,"\n");
			// }

			// float *rtn = (float *) malloc(3* sizeof(float));
			// get_position(p, N, output);
			
			// -----------------------------------//
			// 	  PARTICLE FILTER ENDS HERE 	  //
			// -----------------------------------//	

			//float *output = particle_filter(motion, measurement, num_particles, particles, fp, false);

		}
	}
	simulation_time = read_timer() - simulation_time;
	// float msec  = diff * 1000 / CLOCKS_PER_SEC;
	printf("num_motions = %i, num_particles = %i\n", num_motions, num_particles);
	printf("Time taken %f seconds\n", simulation_time);
	printf("Movement time: %f seconds\n", movement_time);
	printf("Measurement time: %f seconds\n", measurement_time);
	printf("Resampling time: %f seconds\n", resampling_time);
	printf("Reassignment time: %f seconds\n", reassignment_time);
	fclose(fp);    

}