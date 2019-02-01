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
	inline Point& operator+(const  Point& p) const;
	inline double& operator[](size_t index) const;
	inline void operator=(const  Point& p);
	inline void operator=(const  Point&& p);
	inline Point& operator/(const Point& p) const;
	friend std::ostream& operator<< (std::ostream& stream, const Point& point);
	inline Point& operator/(double denominator) const;
	inline Point& operator*(double multiplier) const;
	friend Point& operator*(double multiplier, Point& point);
	inline Point& Point::operator+(double val);
	inline friend Point& operator+(double val, Point& point);
	inline Point& Point::operator+=(const  Point& p);
	static double euclideanDistance(Point& p1, Point& p2);
private:
	double* m_coordinates;
	int m_dimensional;
	int getDimensional() const;
	double* getCoordinate() const;
};

