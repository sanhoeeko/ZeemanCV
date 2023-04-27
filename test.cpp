#include<iostream>
#include"algo.h"
#include"csv.h"
#include<matplotlibcpp.h>
#include<vector>
using namespace std;
namespace plt = matplotlibcpp;

void main_process(string filename) {
	jarpt* points_xy = collect_points(filename);
	Coord center = search(points_xy);
	cout << center.x << ", " << center.y << endl;
	center = regionSearch(points_xy, center, 25); // 2601 searches
	cout << center.x << ", " << center.y << endl;
	/*
	*  Experiments shows that PSO algorithm cannot find the best solution,
	*  if we eager the best and we have some priori knowledge (i.e. the center of
	*  circles is near the center of the original picture.
	* 
	srand(time(0));
	Coord center = optimizeSearch(points_xy, 10, 30);
	cout << center.x << ", " << center.y << endl;
	*/
	jard* points_rt = r_transform(points_xy, center.x, center.y);
	int max_r = int(sqrt(shape.x * shape.x + shape.y * shape.y)) + 1;
	double* spe = collect_spectrum(points_rt, max_r, RESOLUTION);
	vector<double> Spe(spe, spe + max_r * RESOLUTION);
	array_to_csv(filename + ".csv", Spe);
	/*
	plt::plot(Spe);
	plt::show();
	*/
	delete points_xy, points_rt, spe;
}

int main() {
	string files[3] = { "data/1.bmp","data/2.bmp", "data/3.bmp" };
	for (int i = 0; i < 3; i++) {
		cout << files[i] << endl;
		main_process(files[i]);
	}
}