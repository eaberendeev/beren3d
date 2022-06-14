#include "Particles.h"
#include "Shape.h"

inline double layer_resumption_left( const Region& region){
	return Dx * (region.dampCells[0].x() + 1);
}
inline double layer_resumption_right( const Region& region){
	return Dx * (region.numCells.x() - region.dampCells[1].x() - 1);
}

inline int iround(double x){
    if ( x < 0 ) x -= 0.5;
        else x += 0.5;
    return (int) x;
}

void CopyToBuf(const Particle& particleSource, Array<Particle>& ParticlesBuf, const Region& domain){
  Particle particle = particleSource;
  particle.set_global(domain);
  ParticlesBuf.push_back(particle);
}

void addFromBuf(Array<Particle>& particlesData,Array<Particle>& ParticlesBuf, const Region& domain){
  	Particle particle;
 	while( ParticlesBuf.size() > 0){
  	   	particle = ParticlesBuf.back();
  	   	ParticlesBuf.pop_back();
		particle.set_local(domain);
		particlesData.push_back(particle);
  }
}

void MPIExchangeParticles(Array<Particle>& ParticlesBufLeft, Array<Particle>& ParticlesBufRight, const MPI_Topology& MPIconf){
  const long maxSize = ParticlesBufLeft.capacity();
  static Array<Particle> RecvBufLeft(maxSize);
  static Array<Particle> RecvBufRight(maxSize);
  long guestLeft, guestRight, sizeLeft,sizeRight;
  int tag = 2;
  long k;
  long sizeParticleStruct = iround(sizeof(Particle) / sizeof(double));
  MPI_Status status;
  
	auto right = MPIconf.next_line();
  
  	auto left = MPIconf.prev_line();


	sizeLeft = ParticlesBufLeft.size();
	sizeRight = ParticlesBufRight.size();

  MPI_Sendrecv(&sizeLeft, 1, MPI_LONG, left, tag, 
                           &guestRight, 1, MPI_LONG, right, tag, 
                           MPIconf.comm_line(), &status);
  MPI_Sendrecv(&sizeRight, 1, MPI_LONG, right, tag, 
                           &guestLeft, 1, MPI_LONG, left, tag, 
                           MPIconf.comm_line(), &status);

  //MPI_Barrier(MPIconf.CommLine);
  MPI_Sendrecv(&ParticlesBufLeft(0), (sizeLeft + 1) *  sizeParticleStruct, MPI_DOUBLE_PRECISION, left, tag, 
                           &RecvBufRight(0), (guestRight + 1) *  sizeParticleStruct, MPI_DOUBLE_PRECISION, right, tag, 
                           MPIconf.comm_line(), &status);
  MPI_Sendrecv(&ParticlesBufRight(0), (sizeRight + 1) *  sizeParticleStruct, MPI_DOUBLE_PRECISION, right, tag, 
                           &RecvBufLeft(0), (guestLeft + 1) *  sizeParticleStruct, MPI_DOUBLE_PRECISION, left, tag, 
                           MPIconf.comm_line(), &status);
   // MPI_Barrier(MPIconf.CommLine);

	ParticlesBufLeft.clear();
	ParticlesBufRight.clear();
  	
  	for(k = 0; k < guestLeft; ++k){
  		ParticlesBufLeft.push_back(RecvBufLeft(k) );
  	}
   
  	for(k = 0; k < guestRight; ++k){
      	ParticlesBufRight.push_back(RecvBufRight(k) ); 
  	}

}
void current_in_cell(double3 r, double3 r1, long indx, long indy, long indz,double koeff, Array3D<double3>& fieldJ){
	double dx = r1.x() - r.x();
	double dy = r1.y() - r.y();
	double dz = r1.z() - r.z();
	double sx =  (0.5*(r1.x() + r.x() ) - (indx-1)*Dx ) / Dx;
	double sy =  (0.5*(r1.y() + r.y() ) - (indy-1)*Dy ) / Dy;
	double sz =  (0.5*(r1.z() + r.z() ) - (indz-1)*Dz ) / Dz;
	//std::cout << r << " " << r1 << " " << long3(indx,indy,indz) << " " << double3(sx,sy,sz) << "\n";
	fieldJ(indx,indy,indz ).x()     += koeff * dx * (  (1-sy ) * (1-sz ) + dy*dz/(12*Dy*Dz) ); 
	fieldJ(indx,indy,indz+1 ).x()   += koeff * dx * (  (1-sy ) * sz      - dy*dz/(12*Dy*Dz) ); 
	fieldJ(indx,indy+1,indz ).x()   += koeff * dx * (  sy      * (1-sz ) - dy*dz/(12*Dy*Dz) ); 
	fieldJ(indx,indy+1,indz+1 ).x() += koeff * dx * (  sy      * sz      + dy*dz/(12*Dy*Dz) ); 

	fieldJ(indx,indy,indz ).y()     += koeff * dy * (  (1-sx ) * (1-sz ) + dx*dz/(12*Dx*Dz) ); 
	fieldJ(indx,indy,indz+1 ).y()   += koeff * dy * (  (1-sx ) * sz      - dx*dz/(12*Dx*Dz) ); 
	fieldJ(indx+1,indy,indz ).y()   += koeff * dy * (  sx      * (1-sz ) - dx*dz/(12*Dx*Dz) ); 
	fieldJ(indx+1,indy,indz+1 ).y() += koeff * dy * (  sx      * sz      + dx*dz/(12*Dx*Dz) ); 

	fieldJ(indx,indy,indz ).z()     += koeff * dz * (  (1-sx ) * (1-sy ) + dx*dy/(12*Dx*Dy) ); 
	fieldJ(indx,indy+1,indz ).z()   += koeff * dz * (  (1-sx ) * sy      - dx*dy/(12*Dx*Dy) ); 
	fieldJ(indx+1,indy,indz ).z()   += koeff * dz * (  sx      * (1-sy ) - dx*dy/(12*Dx*Dy) ); 
	fieldJ(indx+1,indy+1,indz ).z() += koeff * dz * (  sx      * sy      + dx*dy/(12*Dx*Dy) ); 
}
void push_pic(double3& POS, double3& PULS, long q, double mass, double mpw, const Array3D<double3>& fieldE, \
		   const Array3D<double3>& fieldB, Array3D<double3>& fieldJ, double3 E = double3(0.,0.,0.)){
	long indx, indy,indz,indx1, indy1,indz1,indx_n, indy_n,indz_n;
	double gama;
	double3 US,U1,U2,T,C;
	double xn, yn,zn;
	double xx,yy,zz;
	double xx_n, yy_n, zz_n;	
	double a,b;
	double sx0,sy0,sz0,sdx0,sdy0,sdz0;
	double sx1,sy1,sz1,sdx1,sdy1,sdz1;
	double xp,yp,zp,s;
	alignas(64) double3 B;

	double3 POS1;
	const double rdx = 1. / Dx;
	const double rdy = 1. / Dy;
	const double rdz = 1. / Dz;
	const double dtp = 0.5 * Dt;
	const double koeff = q / (Dt*Dx*Dy*Dz) * mpw;
	int dir; 
	xx = POS.x() * rdx;
	yy = POS.y() * rdy;
	zz = POS.z() * rdz;
	
	indx = long(xx)+1.;
	indy = long(yy)+1.;
	indz = long(zz)+1.;


	indx1 = long(xx+0.5);
	indy1 = long(yy+0.5);
	indz1 = long(zz+0.5);	
    sx1 = (xx - indx+1);
    sy1 = (yy - indy+1);
    sz1 = (zz - indz+1);
    sdx1 = (xx - indx1+0.5);
    sdy1 = (yy - indy1+0.5);
    sdz1 = (zz - indz1+0.5);


    sx0 = 1. - sx1;
    sy0 = 1. - sy1;
    sz0 = 1. - sz1;
    sdx0 = 1. - sdx1;
    sdy0 = 1. - sdy1;
    sdz0 = 1. - sdz1;

	E.x() += sdx0 * ( sy0 * ( sz0 * fieldE(indx1,indy,indz).x() + sz1 * fieldE(indx1,indy,indz+1).x() ) 
		          + sy1 * ( sz0 * fieldE(indx1,indy+1,indz).x() + sz1 * fieldE(indx1,indy+1,indz+1).x() ) ) 
			+ sdx1 * ( sy0 * ( sz0 * fieldE(indx1+1,indy,indz).x() + sz1 * fieldE(indx1+1,indy,indz+1).x() ) 
		          + sy1 * ( sz0 * fieldE(indx1+1,indy+1,indz).x() + sz1 * fieldE(indx1+1,indy+1,indz+1).x() ) );

	E.y() += sx0 * ( sdy0 * ( sz0 * fieldE(indx,indy1,indz).y() + sz1 * fieldE(indx,indy1,indz+1).y() ) 
		          + sdy1 * ( sz0 * fieldE(indx,indy1+1,indz).y() + sz1 * fieldE(indx,indy1+1,indz+1).y() ) ) 
			+ sx1 * ( sdy0 * ( sz0 * fieldE(indx+1,indy1,indz).y() + sz1 * fieldE(indx+1,indy1,indz+1).y() ) 
		          + sdy1 * ( sz0 * fieldE(indx+1,indy1+1,indz).y() + sz1 * fieldE(indx+1,indy1+1,indz+1).y() ) );

	E.z() += sx0 * ( sy0 * ( sdz0 * fieldE(indx,indy,indz1).z() + sdz1 * fieldE(indx,indy,indz1+1).z() ) 
		          + sy1 * ( sdz0 * fieldE(indx,indy+1,indz1).z() + sdz1 * fieldE(indx,indy+1,indz1+1).z() ) ) 
			+ sx1 * ( sy0 * ( sdz0 * fieldE(indx+1,indy,indz1).z() + sdz1 * fieldE(indx+1,indy,indz1+1).z() ) 
		          + sy1 * ( sdz0 * fieldE(indx+1,indy+1,indz1).z() + sdz1 * fieldE(indx+1,indy+1,indz1+1).z() ) );

	B.x() += sx0 * ( sdy0 * ( sdz0 * fieldB(indx,indy1,indz1).x() + sdz1 * fieldB(indx,indy1,indz1+1).x() ) 
		           + sdy1 * ( sdz0 * fieldB(indx,indy1+1,indz1).x() + sdz1 * fieldB(indx,indy1+1,indz1+1).x() ) ) 
		   + sx1 * ( sdy0 * ( sdz0 * fieldB(indx+1,indy1,indz1).x() + sdz1 * fieldB(indx+1,indy1,indz1+1).x() ) 
		           + sdy1 * ( sdz0 * fieldB(indx+1,indy1+1,indz1).x() + sdz1 * fieldB(indx+1,indy1+1,indz1+1).x() ) );

	B.y() += sdx0 * ( sy0 * ( sdz0 * fieldB(indx1,indy,indz1).y() + sdz1 * fieldB(indx1,indy,indz1+1).y() ) 
		          + sy1 * ( sdz0 * fieldB(indx1,indy+1,indz1).y() + sdz1 * fieldB(indx1,indy+1,indz1+1).y() ) ) 
			+ sdx1 * ( sy0 * ( sdz0 * fieldB(indx1+1,indy,indz1).y() + sdz1 * fieldB(indx1+1,indy,indz1+1).y() ) 
		          + sy1 * ( sdz0 * fieldB(indx1+1,indy+1,indz1).y() + sdz1 * fieldB(indx1+1,indy+1,indz1+1).y() ) );

	B.z() += sdx0 * ( sdy0 * ( sz0 * fieldB(indx1,indy1,indz).z() + sz1 * fieldB(indx1,indy1,indz+1).z() ) 
		          + sdy1 * ( sz0 * fieldB(indx1,indy1+1,indz).z() + sz1 * fieldB(indx1,indy1+1,indz+1).z() ) ) 
			+ sdx1 * ( sdy0 * ( sz0 * fieldB(indx1+1,indy1,indz).z() + sz1 * fieldB(indx1+1,indy1,indz+1).z() ) 
		          + sdy1 * ( sz0 * fieldB(indx1+1,indy1+1,indz).z() + sz1 * fieldB(indx1+1,indy1+1,indz+1).z() ) );


	U1 = PULS + q * dtp * E;
	a = q * dtp / sqrt(1. + dot(U1,U1) );
	T = a * B;
	b = 2. / (1. + dot(T,T) );
	C = b * T;
		
	US = U1 + cross(U1,T);
		
	U2 = U1 + cross(US,C);
		
	PULS = U2 + q * dtp * E;
		
	gama = 1. / sqrt(mass * mass + dot(PULS,PULS) );
			
	xn = POS.x() + Dt * PULS.x() * gama;
	yn = POS.y() + Dt * PULS.y() * gama;
	zn = POS.z() + Dt * PULS.z() * gama;
	POS1 = double3(xn,yn,zn);
	xx_n = POS1.x() * rdx;
	yy_n = POS1.y() * rdy;
	zz_n = POS1.z() * rdz;

	indx_n = long(xx_n+1.);
	indy_n = long(yy_n+1.);
	indz_n = long(zz_n+1.);
	bool mx = indx - indx_n;
	bool my = indy - indy_n;
	bool mz = indz - indz_n;
	dir = 4*mx + 2*my + mz;
	if (!dir){
		current_in_cell(POS, POS1, indx, indy, indz,koeff, fieldJ);
	}
	else{
		switch (dir){
			case 1: // z cell intersection
				zp = Dz*(0.5*(indz_n + indz) - 0.5);
				s = (zp - POS.z()) / (POS1.z() - POS.z() );
				xp = POS.x() + s*(POS1.x()-POS.x() );
				yp = POS.y() + s*(POS1.y()-POS.y() );
				break;
			case 2: // y cell intersection
				yp = Dy*(0.5*(indy_n + indy) - 0.5);
				s = (yp - POS.y())/(POS1.y() - POS.y() );
				xp = POS.x() + s*(POS1.x()-POS.x() );
				zp = POS.z() + s*(POS1.z()-POS.z() );
				break;
			case 3: //z and y intersection
				zp = Dz*(0.5*(indz_n + indz) - 0.5);
				yp = Dy*(0.5*(indy_n + indy) - 0.5);
				s = ( (POS1.z()-POS.z()) * (zp-POS.z()) + (POS1.y()-POS.y()) * (yp-POS.y()) )
				/ ( (POS1.z()-POS.z()) * (POS1.z()-POS.z()) + (POS1.y()-POS.y()) * (POS1.y()-POS.y()) );
				xp = POS.x() + s*(POS1.x()-POS.x() );
				break;
			case 4: // x cell intersection 
				xp = Dx*(0.5*(indx_n + indx) - 0.5);
				s = (xp - POS.x())/(POS1.x() - POS.x() );
				yp = POS.y() + s*(POS1.y()-POS.y() );
				zp = POS.z() + s*(POS1.z()-POS.z() );
				break;
			case 5: //z and x intersection
				zp = Dz*(0.5*(indz_n + indz) - 0.5);
				xp = Dx*(0.5*(indx_n + indx) - 0.5);
				s = ( (POS1.z()-POS.z()) * (zp-POS.z()) + (POS1.x()-POS.x()) * (xp-POS.x()) )
				/ ( (POS1.z()-POS.z()) * (POS1.z()-POS.z()) + (POS1.x()-POS.x()) * (POS1.x()-POS.x()) );
				yp = POS.y() + s*(POS1.y()-POS.y() );	
				break;		
			case 6: //x and y intersection
				xp = Dx*(0.5*(indx_n + indx) - 0.5);
				yp = Dy*(0.5*(indy_n + indy) - 0.5);
				s = ( (POS1.x()-POS.x()) * (xp-POS.x()) + (POS1.y()-POS.y()) * (yp-POS.y()) )
				/ ( (POS1.x()-POS.x()) * (POS1.x()-POS.x()) + (POS1.y()-POS.y()) * (POS1.y()-POS.y()) );
				zp = POS.z() + s*(POS1.z()-POS.z() );
				break;
			case 7: // x,y,z intersection
				xp = Dx*(0.5*(indx_n + indx) - 0.5);
				zp = Dz*(0.5*(indz_n + indz) - 0.5);
				yp = Dy*(0.5*(indy_n + indy) - 0.5);	
				break;
			default: 
				xp = -1;
				yp = -1;
				zp = -1;
				std::cout << "Error in calculation cell intersection\n";
		}
		current_in_cell(POS, double3(xp,yp,zp), indx, indy, indz,koeff, fieldJ);
		current_in_cell(double3(xp,yp,zp), POS1, indx_n, indy_n, indz_n,koeff, fieldJ);
	}


	POS = POS1;
}




void push(double3& POS, double3& PULS, long q, double mass, double mpw, const Array3D<double3>& fieldE, \
		   const Array3D<double3>& fieldB, Array3D<double3>& fieldJ, double3 E = double3(0.,0.,0.)){
  	constexpr auto SMAX = 2*SHAPE_SIZE;
	long xk, yk, zk, n, m, k,indx, indy,indz;
	double xx, yy, zz,arg;
	double snm1,snm2,snm3,snm12,snm13,snm23,gama;
	double3 US,U1,U2,T,C;
	double xn, yn,zn;
	double a,b;
	alignas(64) double sx[SMAX], sy[SMAX], sz[SMAX];
	alignas(64) double sdx[SMAX], sdy[SMAX], sdz[SMAX];
	alignas(64) double sx_n[SMAX], sy_n[SMAX], sz_n[SMAX];
	alignas(64) double jx[SMAX][SMAX][SMAX];
	alignas(64) double jy[SMAX][SMAX][SMAX];
	alignas(64) double jz[SMAX][SMAX][SMAX];
	alignas(64) double3 B;

	const double rdx = 1. / Dx;
	const double rdy = 1. / Dy;
	const double rdz = 1. / Dz;
	const double dtp = 0.5 * Dt;
	const double conx = Dx / (6*Dt) * mpw;
	const double cony = Dy / (6*Dt) * mpw;
	const double conz = Dz / (6*Dt) * mpw;
			
	xx = POS.x() * rdx;
	yy = POS.y() * rdy;
	zz = POS.z() * rdz;

	xk = long(xx);
	yk = long(yy);
	zk = long(zz);
	
	for(n = 0; n < SMAX; ++n){
		arg = -xx + double(xk - CELLS_SHIFT + n);
		sx[n] = Shape(arg)/ Dx;
		sdx[n] = Shape(arg + 0.5)/ Dx;
		arg = -yy + double(yk - CELLS_SHIFT + n);
		sy[n] = Shape(arg)/ Dy;
		sdy[n] = Shape(arg + 0.5)/ Dy;
		arg = -zz + double(zk - CELLS_SHIFT + n);
		sz[n] = Shape(arg)/ Dz;
		sdz[n] = Shape(arg + 0.5)/ Dz;
	}
		
	for(n = 0; n < SMAX; ++n){
			indx = xk + n;

		for(m = 0; m < SMAX; ++m){
			indy = yk  + m;
			for(k = 0; k < SMAX; ++k){

				jx[n][m][k] = 0.;
				jy[n][m][k] = 0.;
				jz[n][m][k] = 0.;
				
				snm1 = sdx[n] * sy[m] * sz[k];
				snm2 = sx[n] * sdy[m] * sz[k];
				snm3 = sx[n] * sy[m] * sdz[k];
				snm12 = sdx[n] * sdy[m] * sz[k];
				snm13 = sdx[n] * sy[m] * sdz[k];
				snm23 = sx[n] * sdy[m] * sdz[k];
				
				indz = zk  + k;
				E += Dx * Dy * Dz * (double3(snm1,snm2,snm3) * fieldE(indx,indy,indz) );
				B += Dx * Dy * Dz * (double3(snm23,snm13,snm12) * fieldB(indx,indy,indz) );
			}
		}
	}

	U1 = PULS + q * dtp * E;
	a = q * dtp / sqrt(1. + dot(U1,U1) );
	T = a * B;
	b = 2. / (1. + dot(T,T) );
	C = b * T;
		
	US = U1 + cross(U1,T);
		
	U2 = U1 + cross(US,C);
		
	PULS = U2 + q * dtp * E;
		
	gama = 1. / sqrt(mass * mass + dot(PULS,PULS) );
			
	xn = POS.x() + Dt * PULS.x() * gama;
	yn = POS.y() + Dt * PULS.y() * gama;
	zn = POS.z() + Dt * PULS.z() * gama;
	POS = double3(xn,yn,zn);

	for(n = 0; n < SMAX; ++n){
		arg = -xn * rdx + double(xk - CELLS_SHIFT + n);
		sx_n[n] = Shape(arg)/ Dx;
		arg = -yn * rdy + double(yk - CELLS_SHIFT + n);
		sy_n[n] = Shape(arg)/ Dy;
		arg = -zn * rdz + double(zk - CELLS_SHIFT + n);
		sz_n[n] = Shape(arg)/ Dz;
	}

	for(n = 0; n < SMAX; ++n){
		indx = xk  + n;
		for(m = 0; m < SMAX; ++m){
			indy = yk + m;
			for(k = 0; k < SMAX; ++k){
		  
				if(n == 0) jx[n][m][k] = -q * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));

				if(n > 0 && n < SMAX-1) jx[n][m][k] = jx[n-1][m][k] - q * conx * (sx_n[n] - sx[n]) *  (sy_n[m] * (2*sz_n[k] + sz[k]) + sy[m] * (2 * sz[k] + sz_n[k]));
				
				if(m == 0) jy[n][m][k] = -q * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
				if(m > 0 && m < SMAX-1) jy[n][m][k] = jy[n][m-1][k] -q * cony * (sy_n[m] - sy[m]) *(sx_n[n] * (2*sz_n[k] + sz[k]) + sx[n] * (2 * sz[k] + sz_n[k]));
				
				if(k == 0) jz[n][m][k] = -q * conz * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
				if(k > 0 && k < SMAX-1) jz[n][m][k] = jz[n][m][k-1]-q * conz  * (sz_n[k] - sz[k]) * (sy_n[m] * (2*sx_n[n] + sx[n]) + sy[m] * (2 * sx[n] + sx_n[n]));
			
				indz = zk + k;
			
				fieldJ(indx,indy,indz) += double3(jx[n][m][k],jy[n][m][k],jz[n][m][k]);
			}
		}
	}
		
}

void ParticlesArray::bound_resumption(const Particle& particle, const double3& r_new, const double3& p_new){
	if( option.boundResumption != 1) return;

	Particle newParticle;
	double layer;
	layer = layer_resumption_left(_world.region);
	
	bool addFromLeft = _world.region.boundType[0].x() == OPEN 
		  				&& particle.coord.x() <= layer && r_new.x() > layer;
	
	layer = layer_resumption_right(_world.region);
	
	bool addFromRight = _world.region.boundType[1].x() == OPEN 
		  				&& particle.coord.x() >= layer && r_new.x() < layer;
		  
	if (addFromLeft || addFromRight){
		    
		    newParticle = particle;
			newParticle.coord.x() = (addFromLeft ? r_new.x()-Dx : r_new.x()+Dx);
			newParticle.coord.y() = r_new.y();
			newParticle.coord.z() = r_new.z();
			newParticle.pulse = double3(p_new.x(),Gauss(temperature), Gauss(temperature) );
			particlesData.push_back(newParticle);		    
	}		
}

void ParticlesArray::move(Mesh& mesh,long timestep){
	double3 r;
	double3 p;
	double3 r_glob;
	Particle particle;
	bool in_area;
	bool lostXLeft, lostXRight, lostY, lostZ;
	const long ParticlesBufSize = 2*MaxSizeOfParts / (_world.region.numCells.x());
	static Array<Particle> ParticlesBufLeft(ParticlesBufSize);
	static Array<Particle> ParticlesBufRight(ParticlesBufSize);

	if (charge == 0) return;
	
	long k = 0;
	long kmax = size();

	while (k < kmax ) {
		r = particlesData(k).coord;
		p = particlesData(k).pulse;

		double3 ELas;	
		
		for (const auto& las : mesh.lasers){
			r_glob = _world.region.get_coord_glob(r);
			ELas += las.force(r_glob,timestep);
		}
		#if SHAPE == 1
			push_pic(r,p, charge, mass(k), mpw(k), mesh.fieldE, mesh.fieldB, mesh.fieldJ, ELas);
		#else
			push(r,p, charge, mass(k), mpw(k), mesh.fieldE, mesh.fieldB, mesh.fieldJ, ELas);
		#endif

		bound_resumption(particlesData(k),r,p);
		/*if( option.boundResumption == 1){

		  layer = layer_resumption_left(_world.region);
		  bool addFromLeft = _world.region.boundType[0].x() == OPEN 
		  				&& particlesData(k).coord.x() <= layer && r.x() > layer;
		  layer = layer_resumption_right(_world.region);
		  bool addFromRight = _world.region.boundType[1].x() == OPEN 
		  				&& particlesData(k).coord.x() >= layer && r.x() < layer;
		  
		  if (addFromLeft || addFromRight){
		    
		    particle = particlesData(k);
			particle.coord.x() = (addFromLeft ? r.x()-Dx : r.x()+Dx);
			particle.coord.y() = r.y();
			particle.coord.z() = r.z();
			particle.pulse = double3(p.x(),Gauss(temperature), Gauss(temperature) );
			particlesData.push_back(particle);		    
			}
		}*/		
		
		particlesData(k).coord = r;
		particlesData(k).pulse = p;
		k++;
	}
	
	k = 0;
	while (k < size()) {

		r = particlesData(k).coord;
		
		lostXLeft = (r.x() < Dx * _world.region.dampCells[0].x() );
		 
		lostXRight = (r.x() >= Dx*(_world.region.numCells.x() - _world.region.dampCells[1].x() ));
		
		lostY = (r.y() <= 2*Dy || r.y() >= Dy*_world.region.numCells.y() - 2*Dy);
		lostZ = (r.z() <= 2*Dz || r.z() >= Dz*_world.region.numCells.z() - 2*Dz);
		
		in_area = ! (lostXLeft || lostXRight || lostY || lostZ);
		
		if( in_area )  
		  ++k;
		else{
			if( !(lostZ || lostY) ){
				if(lostXLeft && _world.region.boundType[0].x() == NEIGHBOUR){
					CopyToBuf(particlesData(k),ParticlesBufLeft, _world.region);
				}
				if(lostXRight && _world.region.boundType[1].x() == NEIGHBOUR){
					CopyToBuf(particlesData(k),ParticlesBufRight, _world.region);
				}
			}
			particlesData.del(k);
		}
	}

	MPIExchangeParticles(ParticlesBufLeft, ParticlesBufRight,_world.MPIconf);

	addFromBuf(particlesData, ParticlesBufLeft, _world.region);
	addFromBuf(particlesData, ParticlesBufRight, _world.region);

}


void ParticlesArray::move_virt(Mesh& mesh,long timestep){
	bool is_virt_move;
	double x, y,z, xL,xR;
	double3 r;
	double3 p;
	if(charge == 0) return;
	if ( _world.region.boundType[0].x() != OPEN && _world.region.boundType[1].x() != OPEN) return;
	  
	long k = 0;
	

	while (k < particlesData.size() ) {
		is_virt_move = false;
		x = particlesData(k).coord.x();
		xL = Dx*(_world.region.dampCells[0].x()) ;
		xR = Dx*(_world.region.numCells.x() - _world.region.dampCells[1].x());

		if(x < xL + 2.*Dx &&  _world.region.boundType[0].x() == OPEN){
			x = x - 2*Dx;
			is_virt_move = true;
		}
		if(x > xR - 2.*Dx  &&  _world.region.boundType[1].x() == OPEN){
			x = x + 2*Dx;
			is_virt_move = true;
		}
		
		if( is_virt_move ){
		
			y = particlesData(k).coord.y();
			z = particlesData(k).coord.z();
			p = particlesData(k).pulse;

			r = double3(x,y,z);
			#if SHAPE == 1	
				push_pic(r, p, charge, mass(k), mpw(k), mesh.fieldE, mesh.fieldB, mesh.fieldJ);		
			#else
				push(r, p, charge, mass(k), mpw(k), mesh.fieldE, mesh.fieldB, mesh.fieldJ);		
			#endif
		}
		k++;
	}
}
