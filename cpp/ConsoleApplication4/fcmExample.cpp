// ConsoleApplication2.cpp : Defines the entry point for the console application.
//

#include "fcm.h"
#include "point.h"
#include <chrono>
int main()
{
	const int numOfPoints = 8;
	int dimensional = 2;
	Point points[numOfPoints];
	double c1[2] = { 1, 1 };
	Point p1(c1, dimensional);
	points[0] = p1;
	double c2[2] = { 2, 1 };
	Point p2(c2, dimensional);
	points[1] = p2;
	double c3[2] = { 2, 2};
	Point p3(c3, dimensional);
	points[2] = p3;
	double c4[2] = { 10, 9};
	Point p4(c4, dimensional);
	points[3] = p4;
	double c5[2] = { 10, 10 };
	Point p5(c5, dimensional);
	points[4] = p5;
	double c6[2] = { 9, 10 };
	Point p6(c6, dimensional);
	points[5] = p6;
	double c7[2] = { 9, 9};
	Point p7(c7, dimensional);
	points[6] = p7;
	double c8[2] = { 20, 20 };
	Point p8(c8, dimensional);
	points[7] = p8;
	int numOfClusters = 3;
	FCM fcm(numOfClusters);
	// Record start time
	auto start = std::chrono::high_resolution_clock::now();
	fcm.algorithm(points, numOfPoints, dimensional);
	// Record end time
	auto finish = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elapsed = finish - start;
	std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    return 0;
}

