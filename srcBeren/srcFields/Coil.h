#ifndef COIL_H_
#define COIL_H_
#include "World.h"
struct Coil{
    double z0, R, I;
    Coil(double z0,double R,double I):z0(z0),R(R),I(I){}
};


struct CoilsArray{
    std::vector<Coil> coils;
    const int N = 2000;
    const double hp = 2*PI/N;
    double R, z0, I;
    double* cs;
    CoilsArray(){
        auto nCoils = BCoil[0];
        for (auto k =0; k < nCoils; k++){
            z0 = BCoil[1+3*k];
            R = BCoil[2+3*k];
            I = BCoil[3+3*k];
            coils.emplace_back(Coil(z0,R,I));
        }
        cs = new double[N];
        for (auto i = 0; i < N; i++){
            cs[i] = cos(i*hp);
        }
    }
    ~CoilsArray(){
        delete[] cs;
    }

    double get_Bz(double z, double r);
    double get_Br(double z, double r);
    double get_integ_z(double z, double r, double R);
    double get_integ_r(double z, double r, double R);
};
void set_coils( Array3D<double3>& fieldB,const World& world);
#endif 
