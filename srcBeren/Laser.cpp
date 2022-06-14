#include "Laser.h"
#include "const.h"

inline double pow2(double a){
  return a*a;
}
inline double pow3(double a){
  return a*a*a;
}
inline double pow4(double a){
  return a*a*a*a;
}

Laser::Laser(const std::vector<std::string>& vecStringParams){
    for (const auto& line: vecStringParams){
        set_params_from_string(line);
    }

    double r = sqrt(pow2(start.z() - focus.z() ) + pow2(start.x() - focus.x() )); 
    if (r==0) {
        std::string msg = "Error! The initial position of the laser cannot be the focus position!";
        std::cout << msg << "\n";
        throw msg;
    }
    sin_angle = (start.z() - focus.z() ) / r;
    cos_angle = fabs(start.x() - focus.x() ) / r;

}

// Устанавливаем значения основных параметров в переменной Params в соответствии с содержимым сроки
void Laser::set_params_from_string(const std::string& line){
    std::vector<std::string> strvec;

    strvec = split(line, ' ');

    if(strvec[0]=="tau"){
        tau = stod(strvec[1]);
    }

    if(strvec[0]=="w0"){
        w0 = stod(strvec[1]);
    }
    if(strvec[0]=="type"){
        type = strvec[1];
    }
    if(strvec[0]=="y0"){
        y0 = stod(strvec[1]);
    }
    if(strvec[0]=="vg"){
       vg =  stod(strvec[1]);
    }
    if(strvec[0]=="sigma0"){
       sigma0 =  stod(strvec[1]);
    }
    if(strvec[0]=="focus"){
        focus.x() = stod(strvec[1]);
        focus.y() = stod(strvec[2]);
        focus.z() = stod(strvec[3]);
    }
    if(strvec[0]=="start"){
        start.x() = stod(strvec[1]);
        start.y() = stod(strvec[2]);
        start.z() = stod(strvec[3]);
    }
    if(strvec[0]=="a0"){
       a0 = stod(strvec[1]);
    }
    if(strvec[0]=="delay"){
        delay = stod(strvec[1]);
    }
   
}

bool Laser::is_work(long timestep) const{
  double cTime = timestep*Dt;

  return  (cTime >= delay && cTime <= delay + 2 * tau);

}

double3 rotatate_xz(const double3& coord, double angle){
  double x,z;
  x = coord.x() * cos(angle) - coord.z()*sin(angle);
  z = coord.x() * sin(angle) + coord.z()*cos(angle);
  return double3(x,coord.y(),z);
}

double3 Laser::rotatate(const double3& coord) const{
  double x,z;
  x = coord.x() * cos_angle + coord.z()*sin_angle;
  z = -coord.x() * sin_angle + coord.z()*cos_angle;
  return double3(x,coord.y(),z);
}
double3 Laser::rotatate_back(const double3& coord) const{
  double x,z;
  x = coord.x() * cos_angle - coord.z()*sin_angle;
  z = coord.x() * sin_angle + coord.z()*cos_angle;
  return double3(x,coord.y(),z);
}

double3 Laser::transform_to_loc(const double3& coord) const{

  return rotatate(coord - focus);
}
double3 Laser::transform_to_glob(const double3& coord) const{
	return coord + focus;
}
double3 Laser::force(double3 coord, long timestep) const{
   double3 f;
   double3 coord_loc = transform_to_loc(coord); 
   double3 start_loc = transform_to_loc(start);
   double startLas;
   double cTime = timestep*Dt;
   double x = coord_loc.x();
   double y = coord_loc.y();
   double z = coord_loc.z();
   double r = sqrt(y*y+z*z);
   double fx_loc,fr_loc,fx,fy,fz;
   double cosp,sinp;
    
    startLas = (x - start_loc.x() ) / vg + delay; 

    if( cTime  >= startLas && cTime  <= startLas + 2 * tau) {
      const double RR = 0.5 * w0 * sigma0 * sigma0;
      const double sigma = sigma0 * sqrt(1. + pow2( x / RR) );
      const double w = 0.5*PI*(cTime - startLas) / tau;
      const double Psin = sin(w);
      const double Pcos = cos(w);
      fx_loc =  0.5 * PI * pow2(a0) / tau / vg * pow2(sigma0) / pow2(sigma) * pow3(Psin) * Pcos * exp(-2.*pow2( r / sigma));
      fr_loc = r * pow2(a0) * pow2(sigma0) / pow4(sigma) * pow4(Psin) * exp(-2. * pow2(r / sigma));

      if (r!=0.){
        sinp = z/r; 
        cosp = y/r; 
      } else{
        sinp = 0.;
        cosp = 0.;
      }
      fy = fr_loc * cosp;
      fz = fr_loc * sinp;
    	fx = fx_loc;
      f = double3(fx,fy,fz);	
  } 

   return rotatate_back(f);

}
double3 Laser::get_field_coll(double3 coord, long timestep) const{
   double3 f;
   double3 coord_loc = transform_to_loc(coord); 
   double3 start_loc = transform_to_loc(start);
    double x = coord_loc.x();
   double y = coord_loc.y();
   double z = coord_loc.z();
  double r = sqrt(y*y+z*z);
   double fr,fy,fz;

  double startLas;
  double cTime = timestep*Dt;
   double cosp,sinp;

  startLas = (x - start_loc.x() ) / vg + delay; 

  if( cTime >= startLas && cTime <= startLas + 2 * tau) {
      const double RR = 0.5 * w0 * sigma0 * sigma0;
      const double sigma = sigma0 * sqrt(1. + pow2( x / RR) );
      const double w = 0.5*PI*(cTime - startLas) / tau;

      fr =  a0 * w0 * (sigma0 / sigma) * exp( - pow2(r / sigma) ) * pow2(sin(w));
      if (r!=0.){
        sinp = z/r; 
        cosp = y/r; 
      } else{
        sinp = 0.;
        cosp = 0.;
      }
      fy = fr * cosp;
      fz = fr * sinp;
      f = double3(0.,fy,fz);    
  }
   return rotatate_back(f);

}


double Laser::get_Ez(double y, double t) const{
 
    const double s2 = exp( -(y - y0) * (y - y0) / (sigma0 * sigma0) ) * cos(w0 * t);

    const double s1 = sin(0.5 * PI * t / tau);
    
    return a0 * w0 * s1 * s1 * s2;
}
double Laser::get_Ey(double x, double y, long timestep) const{
    double ctime = timestep*Dt - delay;
    double s2;
    if(focus.x() >= 0.){
      x -= focus.x();
      const double Rlenght = 0.5 * w0 * sigma0 * sigma0;
      const double wx = sigma0 * sqrt(1. + x * x / (Rlenght * Rlenght) );
      const double phi = atan( x / Rlenght);
      const double R = x * (1. + Rlenght * Rlenght / (x * x));
      s2 = sqrt(sigma0 / wx) * exp( - (y - y0) * (y - y0) / ( wx * wx)) 
                        * cos(w0 * ctime + phi - w0 * (y - y0) * (y - y0) * 0.5 / R);
    }
    else{
      s2 = exp( -(y - y0) * (y - y0) / (sigma0 * sigma0) ) * sin(w0 * ctime);

    }

    const double s1 = sin(0.5 * PI * ctime / tau);
    
    return a0 * w0 * s1 * s1 * s2;

}
