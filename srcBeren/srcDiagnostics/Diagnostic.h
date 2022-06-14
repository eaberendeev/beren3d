#ifndef DIAGNOSTIC_H_
#define DIAGNOSTIC_H_
#define DIAGNOSTIC_H_
#include "World.h"
#include "Mesh.h"
#include "Particles.h"
#include "Timer.h"
void write_array2D(const Array2D<double>& data, long size1, long size2, const char* filename, const MPI_Topology& MPIconf);
void write_array3D(const Array3D<double>& data, long size1, long size2, long size3, const char* filename, const MPI_Topology& MPIconf);




struct DiagDataParams{
    std::vector<double> radiationDiagRadiuses;
    std::vector<double> zondCoordsCircleX;
    std::vector<double> sliceFieldsPlaneX, sliceFieldsPlaneY, sliceFieldsPlaneZ;
    std::vector<double> sliceRadiationPlaneX, sliceRadiationPlaneY, sliceRadiationPlaneZ;
    std::vector<double3> zondCoords;
    std::vector<double3> zondCoordsLineX,zondCoordsLineY,zondCoordsLineZ;
    std::vector<long> outTime3D;

};

struct RadialDiagData{
    RadialDiagData(const Region &region, double radius): powerRadCircleLine(region.numNodes.x()),
                            powerRadCircle2D(region.numNodes.x(),720), 
                            circleE(region.numNodes.x(),720), 
                            circleB(region.numNodes.x(),720), 
                            powerRadCircle{0.}, radiationDiagRadius{radius} {

    }

    Array1D<double> powerRadCircleLine;
    Array2D<double> powerRadCircle2D;
    Array2D<double3> circleE;
    Array2D<double3> circleB;
    double powerRadCircle;
    double radiationDiagRadius;
    
    void calc_radiation_pointing_circle2D(const Mesh &mesh);


};

struct DiagData{

    
    DiagData(const Region &region){
        //powerRadLine(region ),
          //                  powerRad1D(region.numNodes.x() ){

        std::vector< std::vector<std::string> > vecStringParams;
        read_params_to_string("Diagnostics","./Diagnostics.cfg",vecStringParams);
        for (const auto& line: vecStringParams[0]){
            set_params_from_string(line);
        }
        for( const auto& radius : params.radiationDiagRadiuses){
            std::cout << "RadialDiagData created for radius "<< radius <<"\n";
            radialDiag.emplace_back( RadialDiagData(region, radius) );
        }
        for( const auto& elem : params.sliceRadiationPlaneX){
            std::cout  << "sliceRadiationPlaneX created for coord " << elem <<"\n";
            powerRadPlaneX.emplace_back( Array2D<double>(region.numNodes.y(), region.numNodes.z() ) );
        }
        for( const auto& elem : params.sliceRadiationPlaneY){
            std::cout  << "sliceRadiationPlaneY created for coord " << elem <<"\n";
            powerRadPlaneY.emplace_back( Array2D<double>(region.numNodes.x(), region.numNodes.z() ) );
        }
        for( const auto& elem : params.sliceRadiationPlaneZ){
            std::cout  << "sliceRadiationPlaneZ created for coord " << elem <<"\n";
            powerRadPlaneZ.emplace_back( Array2D<double>(region.numNodes.x(), region.numNodes.y() ) );
        }        
       Reset() ;

    };
    void set_params_from_string(const std::string& line);

    void calc_energy(const Mesh &mesh,const std::vector<ParticlesArray> &species);
    //void calc_radiation_pointing_circle_2D(const Mesh &mesh);
    //void calc_radiation_pointing_planeZ(const Mesh &mesh, const Region &region);
    void calc_radiation_pointing_planeX(Array2D<double>& dataPlaneX, double coordX, const Mesh &mesh, const Region &region);
    void calc_radiation_pointing_planeY(Array2D<double>& dataPlaneY, double coordY, const Mesh &mesh, const Region &region);
    void calc_radiation_pointing_planeZ(Array2D<double>& dataPlaneZ, double coordZ, const Mesh &mesh, const Region &region);

    void Reset(){
        
        for(auto& data : powerRadPlaneX){
            data.clear();
        }
        for(auto& data : powerRadPlaneY){
            data.clear();
        }
        for(auto& data : powerRadPlaneZ){
            data.clear();
        }        

        for( auto& data : radialDiag){
            data.powerRadCircle = 0;
            data.powerRadCircleLine.clear(); 
            data.powerRadCircle2D.clear(); 
            data.circleE.clear(); 
            data.circleB.clear(); 
        }
    };

    //BoundData<double> powerRad;
    //BoundData<double> powerRadAvg;
    //BoundData<Array1D<double> > powerRadLine;

    std::vector< Array2D<double> > powerRadPlaneX;
    std::vector< Array2D<double> > powerRadPlaneY;
    std::vector< Array2D<double> > powerRadPlaneZ;
    std::vector< RadialDiagData > radialDiag;

    std::map<std::string,double> energyParticlesKinetic;
    std::map<std::string,double> energyParticlesInject;
    
    double energyFieldE, energyFieldB;
    DiagDataParams params;

};




struct Writer{
protected:
    const World &_world;
    const Mesh &_mesh;
    const std::vector<ParticlesArray> &_species;
public:
    FILE *fDiagEnergies;
    DiagData diagData;

    Writer(const World &world, const Mesh &mesh,std::vector<ParticlesArray> &species);

    void output( long timestep);
    ~Writer(){
        if( !_world.MPIconf.is_master() ) return;
            fclose(fDiagEnergies);  
    } 

    void write_particles2D(long timestep);
    void write_particles3D(long timestep);
    void write_energies(long timestep);
    void write_radiation_line(long timestep);
    void write_radiation_avg_line(long timestep);
    
    void write_radiation_circle2D(long timestep);
    void write_radiation_planes(long timestep);


    void diag_zond(long timestep);
    //void diag_zond_lineX_bin(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB,long timestep);
    void write_fields_lineX(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, long timestep);
    void write_fields_lineY(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, long timestep);
    void write_fields_lineZ(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, long timestep);
    void write_fields3D(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB,  const long& timestep);
    void write_fields2D_planeX(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, double coord, const long& timestep);
    void write_fields2D_planeY(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, double coord, const long& timestep);
    void write_fields2D_planeZ(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, double coord, const long& timestep);
    void write_fields_circle( long timestep);
    void write_fields2D_circle(const Array2D<double3>& fieldE, const Array2D<double3>& fieldB, long series, const long& timestep);

//    void write_fields2D(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB,const std::string& axes, long index,  const long& timestep);
    void write_fields2D(const Array2D<double3>& fieldE, const Array2D<double3>& fieldB, const long& timestep);

};
#endif 	
