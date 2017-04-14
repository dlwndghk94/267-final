#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h> 
#include <time.h>
#include <cmath>
#include <string.h>
#include "robot.cpp"
#include <time.h>
#include "omp.h"


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

float *particle_filter(float motion[2], float measurement[4], int N, Robot * p, FILE* fp){
	
	// Motion update (prediction)
	// Robot p2[N];

	#pragma omp parallel for
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
		fprintf(fp,"%f, %f, ",p[i].x,p[i].y);	
	}

	fprintf(fp,"\n");

	// filename=strcat(filename,".csv");
	// filename = "output.csv";

	// fp=fopen(filename,"w+");
	// fprintf(fp,)
	// for(int i=0;i<N;i++){
	 
	//     fprintf(fp,"\n%d",i+1);
	 
	//     for(int j=0;j<n;j++)
	 
	//         fprintf(fp,",%d ",a[i][j]);
	 
	//     }

	// fclose(fp);
	float *rtn = (float *) malloc(3* sizeof(float));
	get_position(p, N, rtn);
	return rtn;
}


int main(){
	/*for (int i = 0; i < 10; i++) {
		printf("%f\n", drand48());
	}
	*/

	int num_motions = 1000;
	int num_particles = 500;
	float length = 20.0;
	char* filename = "output.csv";
	FILE* fp = fopen(filename, "a+"); 

	// Initializing particles array
	Robot particles[num_particles];
	for (int i = 0; i < num_particles; i++){
		Robot r;
		r.initialize(length);
		r.set_noise(bearing_noise, steering_noise, distance_noise);
		particles[i] = r;
	}

	// 
	// for (int i = 0; i < num_motions; i++){
		
	// }
	// float measurements[8][4] = {{4.746936, 3.859782, 3.045217, 2.04550},
 //              {3.510067, 2.916300, 2.146394, 1.598332},
 //              {2.972469, 2.407489, 1.588474, 1.611094},
 //               {1.906178, 1.193329, 0.619356, 0.807930},
 //              {1.352825, 0.662233, 0.144927, 0.799090},
 //              {0.856150, 0.214590, 5.651497, 1.062401},
 //               {0.194460, 5.660382, 4.761072, 2.471682},
 //               {5.717342, 4.736780, 3.909599, 2.342536}};
    // float **motion_ptr = (float **) motions;
    // float **measurements_ptr = (float **)measurements;

	Robot car;
	car.initialize(length);
	car.set(20.0, 20.0, 0);
	car.set_noise(bearing_noise, steering_noise, distance_noise);

	float measurement[4];
	float motion[2];

	clock_t start = clock(), diff;
	for (int i = 0; i < num_motions; i++) {
		
		motion[0] = 2.0 *M_PI / 10.0;
		motion[1] = 20.0;
		car = car.move(motion);
		fprintf(fp,"%f, %f, ",car.x,car.y);
		car.sense(measurement, 1);
		float *output = particle_filter(motion, measurement, num_particles, particles, fp);

		// for (int i = 0; i < 3; i++){
  //   		printf("Output %i: %f\n", i,output[i]);
  //  		}
		// float estimate = particle_filter(motion, &car, )
	}
	diff = clock() - start;
	float msec  = diff * 1000 / CLOCKS_PER_SEC;
	printf("Time taken %f milliseconds\n", msec);
	fclose(fp);

    // float *output = particle_filter(motions,8,measurements,50000,20);
    

}