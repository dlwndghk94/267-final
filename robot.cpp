#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h> 
#include <time.h>
#include <cmath>



float max_steering_angle = M_PI/4.0;
float bearing_noise = 1; // Noise parameter: should be included in sense function.
float steering_noise = 0.1; // Noise parameter: should be included in move function.
float distance_noise = 5.0; // Noise parameter: should be included in move function.
//float bearing_noise = 0.1; // Noise parameter: should be included in sense function.
//float steering_noise = 0.1; // Noise parameter: should be included in move function.
//float distance_noise = 5.0; // Noise parameter: should be included in move function.

float tolerance_xy = 15.0; // Tolerance for localization in the x and y directions.
float tolerance_orientation = 0.25; // Tolerance for orientation.

float landmarks[4][2] = {{0.0,100.0}, {0.0,0.0}, {100.0,0.0}, {100.0, 100.0}};
float world_size = 100.0;


double rand_normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand()/RAND_MAX - 1;
            y = 2.0*rand()/RAND_MAX - 1;

            r = x*x + y*y;
        }
        while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r)/r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}


class Robot {
	public:
		float x,y,orientation,length, bearing_noise, steering_noise, distance_noise;

    	void initialize (float length_);
    	void set(float new_x, float new_y,float new_orientation);
    	void set_noise(float b_noise, float s_noise, float d_noise);
    	Robot move(float* motion);
    	void sense(float* Z,int add_noise);
    	float measurement_prob(float* measurements);
};

void Robot::initialize (float length_) {
	x = drand48() * world_size;
	y = drand48() * world_size;
	orientation = drand48() * 2.0 * M_PI;
	length = length_;
	bearing_noise = 0.0;
	steering_noise = 0.0;
	distance_noise = 0.0;
}

void Robot::set (float new_x, float new_y,float new_orientation) {
	if (new_orientation < 0 or new_orientation >= 2 * M_PI){
		printf("Orientation must be in [0..2pi]\n");
		exit(1);
	}
	x = new_x;
	y = new_y;
	orientation = new_orientation;
}

void Robot::set_noise (float b_noise, float s_noise, float d_noise) {
	bearing_noise = b_noise;
	steering_noise = s_noise;
	distance_noise = d_noise;
}

Robot Robot::move(float *motion){
	float L = length;
	float theta = orientation;
	float alpha = motion[0] + rand_normal(0.0, steering_noise);
	float d = motion[1] + rand_normal(0.0, distance_noise);
	float beta = (d/L) * tan(alpha);
	if (fabs(beta) < 0.001){
		x = x + d*cos(theta);
		y = y + d*sin(theta);
	}
	else{
		float R = d/beta;
		float cx = x - sin(theta) * R;
		float cy = y + cos(theta) * R;
		x = cx + sin(theta+beta)*R;
		y = cy - cos(theta+beta)*R;
	}
	orientation = fmod(theta + beta, 2*M_PI);
	if (orientation < 0){
		orientation = orientation + 2*M_PI;
	}
	Robot result;
	result.initialize(L);
	result.set(x,y,orientation);
	result.set_noise(bearing_noise,steering_noise,distance_noise);
	return result;
}

void Robot::sense(float *Z, int add_noise){
	// float Z[4];
	float theta = orientation;
	for(int i = 0; i < 4; i++){
		float* L = landmarks[i];
		float dx = L[1] - x;
		float dy = L[0] - y;
		if (add_noise == 0){
			*(Z+i) = fmod(atan2(dy,dx),2*M_PI) - orientation;
		} else{
			*(Z+i) = fmod(atan2(dy,dx),2*M_PI) - orientation+ rand_normal(0.0, bearing_noise);
		}
	}
	return;
}

float Robot::measurement_prob(float *measurements) {
	// calculate the correct measurement
	float predicted_measurements[4];
	sense(predicted_measurements, 0);

	// compute error
	float error = 1.0;
	float error_bearing;
	int num_measurements = 4;
	for (int i = 0; i < num_measurements; i++) {
		// printf("measurements: %f\n",measurements[i]);
		// printf("predicted_measurements: %f\n",predicted_measurements[i]);
		error_bearing = fabs(measurements[i] - predicted_measurements[i]);
		error_bearing = fmod((error_bearing + M_PI), (2.0*M_PI)) - M_PI;
		// printf("%i: %f\n",i,error_bearing);

		// update Gaussian
		error *= exp(-(error_bearing*error_bearing) / (bearing_noise * bearing_noise) / 2.0) / 
				sqrt(2.0 * M_PI * bearing_noise * bearing_noise);
	}

	return error;
}





// int main(){
// 	/*for (int i = 0; i < 10; i++) {
// 		printf("%f\n", drand48());
// 	}
// 	*/
// 	Robot p;
// 	p.initialize(20.0);
// 	printf("initializing...\n");
// 	printf("x = %f \n y = %f \n orientation = %f\n length = %f\n bearing_noise = %f\n steering_noise = %f\n distance_noise=%f\n", p.x, p.y, p.orientation, p.length, p.bearing_noise, p.steering_noise, p.distance_noise);

// 	p.set(20, 20, 0);
// 	p.set_noise(1, 1, 1);
// 	printf("\nsetting...\n");
// 	printf("x = %f \n y = %f \n orientation = %f\n length = %f\n bearing_noise = %f\n steering_noise = %f\n distance_noise=%f\n", p.x, p.y, p.orientation, p.length, p.bearing_noise, p.steering_noise, p.distance_noise);

// 	float motion[2] = {0, 10};
// 	p.move(motion);
// 	printf("\nmoving...\n");
// 	printf("x = %f \n y = %f \n orientation = %f\n length = %f\n bearing_noise = %f\n steering_noise = %f\n distance_noise=%f\n", p.x, p.y, p.orientation, p.length, p.bearing_noise, p.steering_noise, p.distance_noise);

// 	float measurements[4];
// 	p.sense(measurements, 1);
// 	printf("\nsensing...\n");
// 	printf("z1 = %f\n z2 = %f\n z3 = %f\n z4 = %f\n\n", measurements[0], measurements[1], measurements[2], measurements[3]);

// 	float probability = p.measurement_prob(measurements);
// 	printf("\nmeasurement_prob..\n");
// 	printf("probability = %f\n", probability);
// }

