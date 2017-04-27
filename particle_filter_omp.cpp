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

int main(){
	int num_motions = 30;
	int num_particles = 100000;

	float length = 20.0;
	char* filename = (char*) "output.csv";
	FILE* fp = fopen(filename, "a+");

	double simulation_time = read_timer();
	// Initializing particles array
	Robot *p = (Robot *) malloc(num_particles * sizeof(Robot));
	// srand(1);
	// #pragma omp parallel for
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

	int N = num_particles;
	Robot p3[N];
	float w[N];
	float mw;
	float sum_lst[N];

	float movement_time = 0;
	float measurement_time = 0;
	float resampling_time = 0;
	float reassignment_time = 0;
	
	#pragma omp parallel shared(car)
	{
		float motion[2];
		motion[0] = 2.0 *M_PI / 10.0;
		motion[1] = 20.0;
		float measurement[4];
		for (int t = 0; t < num_motions; t++) {
			
			#pragma omp master 
			{
				movement_time -= read_timer();
				car = car.move(motion);
				car.sense(measurement, 1);
			}
			
			#pragma omp barrier

			// -----------------------------------//
			// 	  PARTICLE FILTER STARTS HERE 	  //
			// -----------------------------------//

			// move all particles
			#pragma omp for
			for (int i = 0; i<N; i++){
				p[i].move(motion);
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
			#pragma omp for
			for (int i = 0; i < N; i++){
				w[i] = p[i].measurement_prob(measurement);
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
			#pragma omp master
			{
				sum_lst[0] = w[0];
				for (int i = 1; i < N; i++){
					sum_lst[i] = sum_lst[i-1] + w[i];
				}
			}

			// Can parallelize this
			#pragma omp parallel
			for (int i = 0; i < N-1; i++)
			{
				sum_lst[i] = sum_lst[i]/sum_lst[N-1];
			}

			#pragma omp single
			{
				sum_lst[N-1] = 1;
			}

			#pragma omp barrier

			#pragma omp for
			for(int i =0; i < N; i++){
				unsigned int seed = (unsigned int) i;
				float rand_num = (double)rand_r(&i) / (double)RAND_MAX;
				int index = N/2;
				int step_size = N/2;
				while (1){
					if (step_size >1){
						step_size = step_size /2;
					}
					if (index == 0){
						break;
					}

					else if (sum_lst[index -1] <= rand_num && sum_lst[index] > rand_num){
						break;

					}

					else if (sum_lst[index] > rand_num){
						index = index - step_size;
					}
					else if (sum_lst[index] < rand_num){
						index = index + step_size;
					}
					else { // equality
						break;
					}
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

			// float *rtn = (float *) malloc(3* sizeof(float));
			// get_position(p, N, output);
			
			// -----------------------------------//
			// 	  PARTICLE FILTER ENDS HERE 	  //
			// -----------------------------------//	

			//float *output = particle_filter(motion, measurement, num_particles, particles, fp, false);

		}
	}
	simulation_time = read_timer() - simulation_time;
	printf("num_motions = %i, num_particles = %i\n", num_motions, num_particles);
	printf("Time taken %f seconds\n", simulation_time);
	printf("Movement time: %f seconds\n", movement_time);
	printf("Measurement time: %f seconds\n", measurement_time);
	printf("Resampling time: %f seconds\n", resampling_time);
	printf("Reassignment time: %f seconds\n", reassignment_time);
	fclose(fp);    

}