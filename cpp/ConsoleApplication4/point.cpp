#include "point.h"
#include <stdexcept>

Point::Point() {
	m_coordinates = nullptr;
	m_dimensional = 0;
}

Point::Point(int dimensional) {
	m_coordinates = new double[dimensional];
	for (int i = 0; i < dimensional; i++) { m_coordinates[i] = 0.0; }
	m_dimensional = dimensional;
}

Point::Point(double* coordinates,int dimensional) {
	m_coordinates = new double[dimensional];
	if (m_dimensional < 0) throw std::invalid_argument("received negative value");
	for (int i = 0; i < dimensional; i++) {
		m_coordinates[i] = coordinates[i];
	}
	m_dimensional = dimensional;
}

Point::Point(const Point& point) {
	*this = point;
}

Point::Point(const Point&& point) {
	*this = std::move(point);
}

Point& Point::operator+(const  Point& p) const{
	if (p.getDimensional() != m_dimensional)  throw std::invalid_argument("received negative value");
	double* coordinates = new double[m_dimensional];
	for (int i = 0; i < m_dimensional; i++) {
		coordinates[i] = m_coordinates[i];
		coordinates[i] += p[i];
	}
	return *(new Point(coordinates, m_dimensional));
}

void Point::operator=(const  Point& p)
{
	if (this != nullptr) {
		m_dimensional = p.getDimensional();
		m_coordinates = new double[m_dimensional];
		for (int i = 0; i < m_dimensional; i++) {
			m_coordinates[i] = p[i];
		}
	}
}

void Point::operator=(const  Point&& p)
{
	if (this != nullptr) {
		m_coordinates = p.getCoordinate();
		m_dimensional = p.getDimensional();
	}
}

double& Point::operator[](size_t index) const{
	return m_coordinates[index];
}

Point& Point::operator/(const Point& p) const {
	Point* point = new Point(m_coordinates, m_dimensional);
	for (int i = 0; i < m_dimensional; i++) {
		(*point)[i] = (m_coordinates[i] / p.getCoordinate()[i]);
	}
	return *point;
}

Point& Point::operator/(double denominator) const {
	Point* point = new Point(m_coordinates, m_dimensional);
	for (int i = 0; i < m_dimensional; i++) {
		(*point)[i] = (m_coordinates[i] / denominator);
	}
	return *point;
}

Point& Point::operator*(double multiplier) const {
	Point* point = new Point(m_coordinates, m_dimensional);
	for (int i = 0; i < m_dimensional; i++) {
		(*point)[i] = (m_coordinates[i] * multiplier);
	}
	return *point;
}

double* Point::getCoordinate() const {
	return m_coordinates;
}

int Point::getDimensional() const {
	return m_dimensional;
}

std::ostream& operator<< (std::ostream& stream, const Point& point) {
	for (int j = 0; j < point.getDimensional(); j++) {
		if(point.getCoordinate()[j] < 10)
			stream << point.getCoordinate()[j] << ' ';
		else
			stream << point.getCoordinate()[j];
		//if (j < point.getDimensional() - 1) 
		stream << " | ";
	}
	return stream;
}

Point& operator*(double multiplier, Point& point) {
	return point * multiplier;
}

Point& Point::operator+(double val) {
	Point* point = new Point(m_coordinates, m_dimensional);
	for (int i = 0; i < m_dimensional; i++) {
		(*point)[i] = (m_coordinates[i] + val);
	}
	return *point;
}

Point& operator+(double val, Point& point) {
	return point + val;
}

Point& Point::operator+=(const  Point& p) {
	if (p.getDimensional() != m_dimensional)  throw std::invalid_argument("received negative value");
	for (int i = 0; i < m_dimensional; i++) {
		m_coordinates[i] += p[i];
	}
	return *(this);
}