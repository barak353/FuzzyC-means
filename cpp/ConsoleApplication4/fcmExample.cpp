// ConsoleApplication2.cpp : Defines the entry point for the console application.
//

#include "fcm.h"
#include "point.h"


int main()
{
	const int numOfPoints = 8;
	int dimensional = 3;
	Point points[numOfPoints];
	double c1[3] = { 1, 1, 1 };
	Point p1(c1, dimensional);
	points[0] = p1;
	double c2[3] = { 1, 2, 2 };
	Point p2(c2, dimensional);
	points[1] = p2;
	double c3[3] = { 2, 2, 3 };
	Point p3(c3, dimensional);
	points[2] = p3;
	double c4[3] = { 9, 10, 11 };
	Point p4(c4, dimensional);
	points[3] = p4;
	double c5[3] = { 10, 10, 10 };
	Point p5(c5, dimensional);
	points[4] = p5;
	double c6[3] = { 10, 9, 9 };
	Point p6(c6, dimensional);
	points[5] = p6;
	double c7[3] = { 9, 9, 9 };
	Point p7(c7, dimensional);
	points[6] = p7;
	double c8[3] = { 20, 20, 20 };
	Point p8(c8, dimensional);
	points[7] = p8;
	int numOfClusters = 15;
	FCM fcm(numOfClusters);
    fcm.algorithm(points, numOfPoints, dimensional);
    return 0;
}

