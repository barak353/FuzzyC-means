// stdafx.cpp : source file that includes just the standard includes
// ConsoleApplication2.pch will be the pre-compiled header
// stdafx.obj will contain the pre-compiled type information
#include "fcm.h"
#include "point.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <iomanip>

const double FCM::epsillon = 0.00001;

FCM::FCM(int numOfClusters, int dimenstional, int p, int threshold) : m_clustersCenters(nullptr), m_numOfPoints(0){
    m_numOfClusters = numOfClusters;
	m_clustersCenters = new Point[m_numOfClusters];
	Point* temp = new Point(dimenstional);
	for (int i = 0; i < numOfClusters; i++) { m_clustersCenters[i] = *(temp); };
	delete temp;

    // The fuzziness.
    m_p = p;
    // The maximum number of iteration of the algorithm.
    m_threshold = threshold;
    // W is the weights matrix.
    m_W = nullptr;
	srand(time(NULL));
}

FCM::~FCM() {
    //delete all data allocation.

}

void FCM::InitMembershipMatrix(int numOfPoints) {
    /*  This function initialize the weights matrix.
    :param num_of_points: number of points define the size of the W matrix.*/
    //Initialize m_W ti matrix numOfPoints X m_numOfClusters with zeros.
    int rows = m_numOfPoints = numOfPoints, cols = m_numOfClusters;
	//delete privious.
	if (m_W != nullptr) {
		for (int i = 0; i < rows; ++i) {
			delete[] m_W[i];
		}
		delete[] m_W;
	}
	//Init new.
	m_W = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        m_W[i] = new double[cols];
        for (int j = 0; j < cols; j++) {
            m_W[i][j] = 0;
        }
    }

    for (int i = 0; i < numOfPoints; i++) {
        double rowSum = 0.0;
        for (int j = 0; j < m_numOfClusters - 1; j++) {
            double ran = rand();
            double randNum = ran / RAND_MAX; //random number between 0 to 1.
                                                //The sum of the weights can't be larger than 1.
            if (randNum + rowSum <= 1.0) {
                m_W[i][j] = randNum;
                rowSum += m_W[i][j];
            }
        }
        //If it's the last iteration (sum can be more than 1.0).
        m_W[i][m_numOfClusters - 1] = 1.0 - rowSum;
    }
//	this->print(m_W, numOfPoints, m_numOfClusters);
}

Point* FCM::ComputeCentroids(Point* points, int n, int dimenstional){
	/*  :param points: data points.
        :param n: num of points.
        :param dimenstional: dimenstional of point.
        :return: weight of the given element
        */
	Point* centers = new Point[m_numOfClusters];
	// For each cluster
	for (int c = 0; c < m_numOfClusters; c++) {
		//Weight.
		Point WeightsOfXvec(dimenstional);
		double sumOfWeights = 0.0;
		// (numerator\denominator) = sigma[0-n]((W[i][c])^2) * X[i]\ sigma[0-n]((W[i][c])^2)
		for (int i = 0; i < n; i++) {
			double denominator = pow(m_W[i][c], m_p);
			Point numerator = denominator * points[i];
			WeightsOfXvec += numerator;
			sumOfWeights += denominator;
		}
		if (sumOfWeights == 0.0) {
			// To prevent division in zero.
			sumOfWeights = 0.000001;
		}
		m_clustersCenters[c] = WeightsOfXvec / sumOfWeights;
	}
	printPoints(m_clustersCenters, m_numOfClusters, dimenstional);
	return m_clustersCenters;
}

double FCM::euclideanDistance(double* x, int* c, int len) {
	/*Compute the Euclidean distance between X and C
	: param x : set of X points
	: param c : the center or a cluster
	: return : the distance*/
	double dis = 0.0;
	for(int i = 0; i < len; i++) {
		dis += pow(x[i] - c[i], 2);
	}
	return sqrt(dis);//Ask Ravit why she missed here the sqrt.
}

void FCM::algorithm(Point* points, int n , int dimenstional) {
	/*:param X: set of points*/
		this->printPoints(points, n, dimenstional);
	bool computeCenters = false;
	if (m_clustersCenters != nullptr) computeCenters = true;
	if (m_W == nullptr) {
		int numOfPoints = n;
		InitMembershipMatrix(numOfPoints);
	}
	this->printW();

	std::vector<Point> centersList;
	int c = 0;
	while(computeCenters) {
		ComputeCentroids(points, n, dimenstional);
		//Point* centers = ComputeCentroids(points, n, m);
		//centersList.push_back(centers[c]);
		/*            self.UpdateWeights(X)
            if centersList.__len__() >= 3:
                if any(x in centersList[centersList.__len__() - 1]  for x in centersList[centersList.__len__() - 3]):
                    computeCenters = False
		c++;*/
	}
	
}

void FCM::printW() const {
	std::cout << "Whight Matrix is:" << endl << "  | ";
	std::cout << std::fixed;
	std::cout << std::setprecision(2);
	for (int i = 0; i < m_numOfClusters; i++) {
		if(i < 10) std::cout << i << "    | ";
		else std::cout << i << "   | ";
	}
	std::cout << endl << "--";
	for (int i = 0; i < m_numOfClusters; i++) {
		std::cout << "-------";
	}
	std::cout << "|" << endl;
	for (int i = 0; i < m_numOfPoints; i++) {
		cout << i << " | ";
		for (int j = 0; j < m_numOfClusters; j++) {
			cout << m_W[i][j];
			if (j < m_numOfClusters - 1) cout << " | ";
		}
		cout << " |" << endl;
	}
}


void FCM::printPoints(Point* points, int n, int m) const{
	std::cout << "Points:" << endl << "  | ";

	for (int i = 0; i < m; i++) {
		std::cout << i << "  | ";
	}
	std::cout << endl << "---";
	for (int i = 0; i < m; i++) {
		std::cout << "----";
	}
	std::cout << "--|";
	cout << endl;
	for (int i = 0; i < n; i++) {

		cout << i << " | ";
		cout << points[i];
		cout << endl;
	}
	cout << endl;
}
