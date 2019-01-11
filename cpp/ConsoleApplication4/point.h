#pragma once

#include <iostream>

class Point {
public:
	Point(double* coordinates, int dimensional);
	Point(const Point& point);
	Point(const Point&& point);
	Point(int dimensional);
	Point();
	int getDimensional(void);
	Point& operator+(const  Point& p) const;
	double& operator[](size_t index) const;
	void operator=(const  Point& p);
	void operator=(const  Point&& p);
	Point& operator/(const Point& p) const;
	friend std::ostream& operator<< (std::ostream& stream, const Point& point);
	Point& operator/(double denominator) const;
	Point& operator*(double multiplier) const;
	friend Point& operator*(double multiplier, Point& point);
	Point& Point::operator+(double val);
	friend Point& operator+(double val, Point& point);
	Point& Point::operator+=(const  Point& p);
	static double euclideanDistance(Point& p1, Point& p2);
private:
	double* m_coordinates;
	int m_dimensional;
	int getDimensional() const;
	double* getCoordinate() const;
};

