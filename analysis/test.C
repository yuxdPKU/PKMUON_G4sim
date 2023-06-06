#include <stdio.h>
#include <math.h>

typedef struct {
    double x;
    double y;
    double z;
} Point3D;

typedef struct {
    Point3D startPoint;
    Point3D direction;
} Line3D;

Point3D calculatePerpendicularMidpoint(const Line3D* line1, const Line3D* line2) {
    Point3D midpoint;
    
    double a = line1->startPoint.x;
    double b = line1->startPoint.y;
    double c = line1->startPoint.z;
    
    double d = line1->direction.x;
    double e = line1->direction.y;
    double f = line1->direction.z;
    
    double g = line2->startPoint.x;
    double h = line2->startPoint.y;
    double i = line2->startPoint.z;
    
    double j = line2->direction.x;
    double k = line2->direction.y;
    double l = line2->direction.z;
    
    double m = a - g;
    double n = b - h;
    double o = c - i;
    
    double t = -(m * j + n * k + o * l) / (d * j + e * k + f * l);
    
    midpoint.x = a + t * d;
    midpoint.y = b + t * e;
    midpoint.z = c + t * f;
    
    return midpoint;
}

int test() {
    Line3D line1 = {{1.0, 1.0, 1.0}, {1.0, 0.0, 0.0}};
    Line3D line2 = {{2.0, 3.0, 0.0}, {0.0, 1.0, 0.0}};
    
    Point3D midpoint = calculatePerpendicularMidpoint(&line1, &line2);
    
    printf("Midpoint: (%f, %f, %f)\n", midpoint.x, midpoint.y, midpoint.z);
    
    return 0;
}

