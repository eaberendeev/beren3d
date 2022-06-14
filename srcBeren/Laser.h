#ifndef LASER_H_
#define LASER_H_
#include "Vec.h"
#include "Read.h"
#include <string>

struct Laser{
    double tau;
    double w0;
    std::string type;
    double vg;
    double3 focus;
    double3 start;
    double angle,sin_angle,cos_angle;
    double sigma0;
    double a0;
    double delay;
    double y0;
    
    Laser(const std::vector<std::string>& vecStringParams);
    void set_params_from_string(const std::string& line);
    bool is_work(long timestep) const;
    double get_Ey(double x, double y, long timestep) const;
    double get_Ez(double y, double times) const;

    double3 force(double3 x, long timestep) const;
    double3 get_field_coll(double3 x, long timestep) const;
    double3 transform_to_loc(const double3& coord) const;
    double3 transform_to_glob(const double3& coord) const;
    double3 rotatate(const double3& coord) const;
    double3 rotatate_back(const double3& coord) const;

};

#endif 
