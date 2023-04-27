#pragma once

#include<stdlib.h>
#include<string>
#include<math.h>
#include"imgProcess.h"
using namespace std;

#define RESOLUTION 10

template<typename ty>
class Jar {
private:
	ty* a;
	int max_cnt;
	int cnt;

public:
	Jar(int volume) {
		a = (ty*)malloc(volume * sizeof(ty));
		max_cnt = volume;
	}
	~Jar() {
		free(a);
	}
	void append(const ty x) { //with copying
		/*
		if (cnt == max_cnt) {
			throw "Jar overflow!\n";
		}
		*/
		// dangerous, but less time-consuming
		a[cnt++] = x;
	}
	ty& operator[](const int idx) {
		return a[idx];
	}
	int size() {
		return cnt;
	}
	ty* ptr() {
		return a;
	}
};

typedef unsigned char uchar;
struct Pt {
	int x, y; uchar val;
};
struct DoublePt {
	double r, theta; uchar val;
};

void convert_point(Jar<Pt>* jar, uchar* mat, int line, int col, uchar threshold) {
	for (int i = 0; i < line; i++) {
		for (int j = 0; j < col; j++) {
			uchar val = mat[i * col + j];
			if (val > threshold) {
				jar->append(Pt{ i,j,val });
			}
		}
	}
}

Shape shape; //global cache, dangerous

Jar<Pt>* collect_points(string filename) {
	uchar* buffer = NULL;
	shape = read(filename, buffer);
	Jar<Pt>* jar = new Jar<Pt>(shape.x * shape.y);
	convert_point(jar, buffer, shape.x, shape.y, 0);
	return jar;
}

typedef Jar<Pt> jarpt;

Jar<DoublePt>* r_transform(jarpt* jar, int center_x, int center_y) {
	int n = jar->size();
	Jar<DoublePt>* newjar = new Jar<DoublePt>(n);
	for (int i = 0; i < n; i++) {
		Pt p = (*jar)[i];
		int x = p.x - center_x;
		int y = p.y - center_y;
		double r = sqrt(x * x + y * y);
		// double theta = atan2(y, x);
		double theta = 0; // atan2 cost most cpu time, but unused
		newjar->append(DoublePt{ r,theta,p.val });
	}
	return newjar;
}

typedef Jar<DoublePt> jard;

double* collect_spectrum(jard* jar, int max_r, int resolution) {
	double* spe = new double[max_r * resolution]();
	int n = jar->size();
	for (int i = 0; i < n; i++) {
		DoublePt p = (*jar)[i];
		int idx = (int)(p.r * resolution);
		spe[idx] += p.val;
	}
	return spe;
}

double mean(double* x, int n) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += x[i];
	}
	double mean = sum / n;
	return mean;
}
double mse(double* x, int n) {
	double mean_val = mean(x, n);
	double sum = 0;
	for (int i = 0; i < n; i++) {
		sum += pow(x[i] - mean_val, 2);
	}
	double mse_val = sum / n;
	return mse_val;
}

int argmax(double* arr, int n){
    int min_idx = 0;
    for(int i=1; i<n; i++){
        if(arr[i] > arr[min_idx]){
            min_idx = i;
        }
    }
    return min_idx;
}

struct Coord {
	int x, y;
};
typedef double Act(jarpt* jar, int center_x, int center_y);

Coord gridSearch(Coord c0, Coord c1, Act act, jarpt* jar, double known_value = -1, int known_binary = -1) {
	double x00, x01, x10, x11;
	x00 = known_binary == 0 ? known_value : act(jar, c0.x, c0.y);
	x01 = known_binary == 1 ? known_value : act(jar, c0.x, c1.y);
	x10 = known_binary == 2 ? known_value : act(jar, c1.x, c0.y);
	x11 = known_binary == 3 ? known_value : act(jar, c1.x, c1.y);
	double arr[4] = {x00, x01, x10, x11};
	int arg = argmax(arr, 4);
	cout << arr[arg] << endl;
	if (abs(c0.x - c1.x) <= 1 && abs(c0.y - c1.y) <= 1) {
		// get out of the recursion
		switch (arg)
		{
		case 0:return c0;
		case 1:return Coord{ c0.x,c1.y };
		case 2:return Coord{ c1.x,c0.y };
		case 3:return c1;
		}
	}
	else {
		// continue the recursion
		int cx = (c0.x + c1.x) / 2;
		int cy = (c0.y + c1.y) / 2;
		if (arg == 0) {
			Coord newc1 = { cx,cy };
			return gridSearch(c0, newc1, act, jar, x00, 0);
		}
		else if (arg == 1) {
			Coord newc0 = { c0.x,cy };
			Coord newc1 = { cx,c1.y };
			return gridSearch(newc0, newc1, act, jar, x01, 1);
		}
		else if (arg == 2) {
			Coord newc0 = { cx,c0.y };
			Coord newc1 = { c1.x,cy };
			return gridSearch(newc0, newc1, act, jar, x10, 2);
		}
		else {
			Coord newc0 = { cx,cy };
			return gridSearch(newc0, c1, act, jar, x11, 3);
		}
	}
}

double single_spectrum(jarpt* jar, int center_x, int center_y) {
	jard* jd = r_transform(jar, center_x, center_y);
	int max_r = int(sqrt(shape.x * shape.x + shape.y * shape.y)) + 1;
	double* spe = collect_spectrum(jd, max_r, RESOLUTION);
	double score = mse(spe, max_r * RESOLUTION);
	delete jd; delete[] spe;
	return score;
}

Coord search(jarpt* jar) {
	Coord c0 = { 0,0 };
	Coord c1 = { shape.x,shape.y };
	return gridSearch(c0, c1, single_spectrum, jar);
}

//dull method
Coord regionSearch(jarpt* jar, Coord center, int boxsize) {
	int n = 2 * boxsize + 1;
	Jar<double> score(n * n);
	for (int i = center.x - boxsize; i <= center.x + boxsize; i++) {
		for (int j = center.y - boxsize; j <= center.y + boxsize; j++) {
			score.append(single_spectrum(jar, i, j));
		}
	}
	int arg = argmax(score.ptr(), n * n);
	cout << score[arg] << endl;
	int x = center.x - boxsize + arg / n;
	int y = center.y - boxsize + arg % n;
	return Coord{ x,y };
}

//particle swarm optimization
#define no_more_than(__x, __a) __x=__x>__a?__a:__x
#define no_less_than(__x, __a) __x=__x<__a?__a:__x
struct Particle
{
	int x, y;
	double x_velocity, y_velocity;
};
Coord optimize(Act act, jarpt* jar, int num_particles, int num_iterations,
	double inertia_weight_max, double inertia_weight_min, double cognitive_weight, double social_weight) {
	// Particle Swarm Optimization algorithm implementation
	// Initialize particles with random positions and velocities
	std::vector<Particle> particles(num_particles);
	std::vector<Particle> best_positions(num_particles);
	double global_best_score = std::numeric_limits<double>::lowest();
	Particle global_best_position;

	for (int i = 0; i < num_particles; i++) {
		particles[i].x = rand() % shape.x;
		particles[i].y = rand() % shape.y;
		double score = act(jar, particles[i].x, particles[i].y);
		if (score > global_best_score) {
			global_best_score = score;
			global_best_position = particles[i];
		}
		best_positions[i] = particles[i];
	}

	double k = (inertia_weight_max - inertia_weight_max) / num_iterations;

	// Update particle velocities and positions iteratively
	for (int iter = 0; iter < num_iterations; iter++) {

		// Update meta-parameters
		double inertia_weight = inertia_weight_max - k * iter;

		for (int i = 0; i < num_particles; i++) {
			// Update velocity
			double r1 = static_cast<double>(rand()) / RAND_MAX;
			double r2 = static_cast<double>(rand()) / RAND_MAX;
			particles[i].x_velocity = inertia_weight * particles[i].x_velocity +
				cognitive_weight * r1 * (best_positions[i].x - particles[i].x) +
				social_weight * r2 * (global_best_position.x - particles[i].x);
			particles[i].y_velocity = inertia_weight * particles[i].y_velocity +
				cognitive_weight * r1 * (best_positions[i].y - particles[i].y) +
				social_weight * r2 * (global_best_position.y - particles[i].y);

			// Update position
			particles[i].x += particles[i].x_velocity;
			particles[i].y += particles[i].y_velocity;
			no_more_than(particles[i].x, shape.x);
			no_less_than(particles[i].x, 0);
			no_more_than(particles[i].y, shape.y);
			no_less_than(particles[i].y, 0);

			// Check if new position is better than previous best
			double score = act(jar, particles[i].x, particles[i].y);
			// cout << "# " << score << endl;
			if (score > act(jar, best_positions[i].x, best_positions[i].y)) {
				best_positions[i] = particles[i];
				if (score > global_best_score) {
					global_best_score = score;
					global_best_position = particles[i];
				}
			}
		}
		//Show score
		cout << global_best_score << endl;
	}

	// Return the globally best position found
	return Coord{ global_best_position.x,global_best_position.y };
}

Coord optimizeSearch(jarpt* jar, int num_particles, int num_iterations) {
	double inertia_max = 8;
	double inertia_min = 0;
	double cognitive = 1;
	double social = 0.05;
	return optimize(single_spectrum, jar, num_particles, num_iterations,
		inertia_max, inertia_min, cognitive, social);
}