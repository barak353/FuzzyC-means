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

const double FCM::Epsilon = 0.00001;

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

Point* FCM::ComputeCentroids(Point* points, int n){
	/*  :param points: data points.
        :param n: num of points.
        :param dimenstional: dimenstional of point.
        :return: weight of the given element
        */
	if (n == 0 || points == nullptr) return nullptr;
	Point* centers = new Point[m_numOfClusters];
	// For each cluster
	for (int c = 0; c < m_numOfClusters; c++) {
		//Weight.
		Point WeightsOfXvec(points[0].getDimensional());
		double sumOfWeights = 0.0;
		// (numerator\denominator) = sigma[0-n]((W[i][c])^2) * X[i]\ sigma[0-n]((W[i][c])^2)
		for (int i = 0; i < n; i++) {
			double denominator = pow(m_W[i][c], m_p);
			Point& numerator = denominator * points[i];
			//std::cout << "numerator: " << numerator << endl;
			WeightsOfXvec += numerator;
			//std::cout << "WeightsOfXvec: " << WeightsOfXvec << endl;

			sumOfWeights += denominator;
			//std::cout << "sumOfWeights: " << sumOfWeights << endl;

		}
		if (sumOfWeights == 0.0) {
			// To prevent division in zero.
			sumOfWeights = 0.000001;
		}
		//std::cout << "WeightsOfXvec: " << WeightsOfXvec << endl;
		//std::cout << "sumOfWeights: " << sumOfWeights << endl;
		//std::cout << "WeightsOfXvec / sumOfWeights: " << WeightsOfXvec / sumOfWeights << endl;

		m_clustersCenters[c] = WeightsOfXvec / sumOfWeights;
		//std::cout << "m_clustersCenters[" << c << "]: " << m_clustersCenters[c] << endl;
	}
	return m_clustersCenters;
}

void FCM::updateWeights(Point* points, size_t n) {
	/*
		update the weights matrix
		: param X : data points
	*/
	for (size_t i = 0; i < n; i++) {
		for (size_t c = 0; c < m_numOfClusters; c++) {
			m_W[i][c] = ComputeWeight(points, i, c);
		}
		this->printW();
	}

}



double FCM::ComputeWeight(Point* points, size_t i, size_t c) {
	/*
		:param X : data points
		: param i : specific element in x
		: param c : specific cluster
		: return : weight of the given element
	*/

	double numerator = Point::euclideanDistance(points[i], m_clustersCenters[c]);
	numerator = pow(numerator, 2);
	//To prevent from return infinity.
	if (numerator == 0)
		numerator = 1.0 - Epsilon;
	numerator = pow(1 / numerator, 1.0 / (m_p - 1));
	double sum = 0.0;
	//(numerator/denominator) = (1 / dist(xi,cj)^(1/p-1)) / sigma[0-numOfClusters](1/dist(xi,cj)^(1/p-1))
	for (size_t q = 0; q < m_numOfClusters; q++) {
		double denominator = Point::euclideanDistance(m_clustersCenters[q], points[i]);
		denominator = pow(denominator, 2);
		if (denominator == 0.0) {
			denominator = 1 / Epsilon;
		}
		sum += pow(1 / denominator, 1.0 / (m_p - 1));
	}
	//To prevent from return infinity.
	if (sum == 0) return 1.0 - Epsilon;
	return (numerator / sum);

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
		Point* centers = ComputeCentroids(points, n);
		printPoints(centers, m_numOfClusters, centers[0].getDimensional());
		updateWeights(points, n);

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
	std::cout << std::setprecision(8);
	for (int i = 0; i < m_numOfClusters; i++) {
		if(i < 10) std::cout << i << "          | ";
		else std::cout << i << "         | ";
	}
	std::cout << endl << "--";
	for (int i = 0; i < m_numOfClusters; i++) {
		std::cout << "-------------";
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
	std::cout << "Clusters:" << endl << "  | ";

	for (int i = 0; i < m; i++) {
		if (i < 10) std::cout << i << "           | ";
		else std::cout << i << "          | ";
	}
	std::cout << endl << "--";
	for (int i = 0; i < m; i++) {
		std::cout << "-------------";
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
