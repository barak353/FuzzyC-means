// stdafx.h : include file for standard system include files,
// or project specific include files that are used frequently, but
// are changed infrequently
//

#pragma once
#include "point.h"
#include <stdio.h>
#include <iostream>

using namespace std;

class FCM {
    /*
    Pseudo code :
    1) initialized W matrix
    2) compute the means(cluster centers)
		3) update the membership matrix
		4) compare the new membership with the old one, is difference is less than a threshold, stop.otherwise
		5) return to step 2)
    */
public:
    FCM(int numOfClusters = 2, int dimenstional = 3, int p = 9, size_t threshold = 30);
    ~FCM();
	void algorithm(Point* points, int n, int m);
private:
	static const double Epsilon;
	static const double Threshold;

    int m_numOfClusters;
	int m_numOfPoints;
    int m_p;
    double** m_W;
	double** m_sW;
	size_t m_threshold;
    Point* m_clustersCenters;
	double ComputeWeight(Point* points, size_t i, size_t c);
	void InitMembershipMatrix(int numOfPoints = 2);
	void ComputeCentroids(Point* points);
	void updateWeights(Point* points);
	void saveWeights(void);
	void printPoints(Point* points, int n, int m) const;
	void printW() const;
	bool checkCondition(void);
};

