#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h> 
#include <time.h>
#include <cmath>
#include <string.h>
#include "robot.cpp"
#include "common.cpp"

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
	srand(1);
	int num_motions = 5;
	int num_particles = 100000;
	float length = 20.0;
	char* filename = (char*) "output.csv";
	FILE* fp = fopen(filename, "a+");

	double simulation_time = read_timer();
	// Initializing particles array
	// Robot p[num_particles];
	Robot *p= (Robot *) malloc(num_particles * sizeof(Robot));
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

	// printf("Starting\n");
	for (int t = 0; t < num_motions; t++) {

		motion[0] = 2.0 *M_PI / 10.0;
		motion[1] = 20.0;
		car = car.move(motion);
		//fprintf(fp,"%f, %f, ",car.x,car.y);
		car.sense(measurement, 1);
		float mw = -999999999.0;

		// -----------------------------------//
		// 	  PARTICLE FILTER STARTS HERE 	  //
		// -----------------------------------//

		// move all particles
		movement_time -= read_timer();
		for (int i = 0; i<N; i++){
			p[i] = p[i].move(motion);
		}
		movement_time += read_timer();
		// printf("Measuring\n");

		// Measurement update
		measurement_time -= read_timer();
		for (int i = 0; i < N; i++){
			w[i] = p[i].measurement_prob(measurement);
			if (mw < w[i]){
				mw = w[i];
			}
		}
		measurement_time += read_timer();

		// Resampling
		resampling_time -= read_timer();

		// bin search number higher.
		// printf("binsearch:\n");
		float sum_lst[N];
		sum_lst[0] = w[0];
		for (int i = 1; i < N; i++){
			sum_lst[i] = sum_lst[i-1] + w[i];
		}
		for (int i = 0; i < N; i++){
			sum_lst[i] = sum_lst[i]/sum_lst[N-1];
		}

		for(int i =0; i < N; i++){
			// if (i % 1000 == 0) {
			// 	printf("Particle %i\n", i);
			// }
			// unsigned int seed = i;
			// printf("RAND_MAX IS %i\n", RAND_MAX);
			// float rand_num = (double)rand_r(&seed) / (double)RAND_MAX;
			float rand_sample = (float) rand();
			float rand_num = rand_sample / (float) RAND_MAX;
			// printf("rand_num is %f\n", rand_num);
			int index = N/2;
			int step_size = N/2;
			int count = 0;

			if (rand_num > 1.0) {
				printf("Error! Random num sampled was %f\n", rand_sample);
			}

			while (1){
				// printf("Rand num is %d\n", rand_num);
				count += 1;
				if (count > N) {
					printf("Error on rand_num %f, rand_sample %f, particle number %i, particle weight %d, index %i\n", rand_num, rand_sample, i, w[i], index);
					printf("Max sum lst idx is %f", sum_lst[99999]);
					printf("Current sum_list val is %f", sum_lst[index]);
				}
				// printf("index: %i\n", index);
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
				else { // sum_lst[index] == rand_num
					break;
				}
			}
			p3[i] = p[index];
		}

		// for (int i = 0; i < N; i++){

		// 	unsigned int seed = 1;
		// 	int index = rand_r(&seed) % N;
		// 	float beta = 0.0;
		// 	float rand_num = (double)rand_r(&seed) / (double)RAND_MAX;
		// 	beta = beta + rand_num * 2.0 * mw;

			
		// 	while( beta > w[index]){
		// 		beta = beta - w[index];
		// 		index = (index +1) % N;
		// 	}
			
		// 	p3[i] = p[index];
		// }
		resampling_time += read_timer();

		//reassignment_time
		reassignment_time -= read_timer();
		for (int i = 0; i < N; i++) {
			p[i] = p3[i];
		}
		reassignment_time += read_timer();

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
	simulation_time = read_timer() - simulation_time;
	// float msec  = diff * 1000 / CLOCKS_PER_SEC;
	printf("num_motions = %i, num_particles = %i\n", num_motions, num_particles);
	printf("Time taken %f seconds (time/4 = %f)\n", simulation_time, simulation_time/4);
	printf("Movement time: %f seconds\n", movement_time);
	printf("Measurement time: %f seconds\n", measurement_time);
	printf("Resampling time: %f seconds\n", resampling_time);
	printf("Reassignment time: %f seconds\n", reassignment_time);
	fclose(fp);    

}