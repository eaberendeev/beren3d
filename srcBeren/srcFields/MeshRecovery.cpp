#include "Mesh.h"

void read_array_from_recovery(const std::string& dataName, Array3D<double3>& field, const MPI_Topology &MPIconf){

  FILE *fFields;
  char filename[50];
  long info1,info2,info3;
  auto size = field.size();
  //auto size_r = field.size_d2();
      
  if( MPIconf.is_master_depth() ){

    sprintf(filename, (dataName + "%03ld").c_str(),MPIconf.rank_line() );  
     
    fFields = fopen(filename, "rb");
        
    size_t result = fread(&info1, sizeof(long), 1, fFields);
    result = fread(&info2, sizeof(long), 1, fFields);
    result = fread(&info3, sizeof(long), 1, fFields);
    
    //assert(result == sizeof(long) && "Number of reading byte is incorrect");

    assert(info1 ==  size.x() && info2 ==  size.y() && info3 == size.z() && "Invalid Field size from read");
      
    result = fread(&field.data(0), 3*size.x()*size.y()*size.z()*sizeof(double), 1, fFields);    
    assert(result == 3*size.x()*size.y()*size.z()*sizeof(double) && "Number of reading byte is incorrect");

    fclose(fFields);

  }

  MPI_Bcast( &field.data(0), 3*size.x()*size.y()*size.z(), MPI_DOUBLE, 0, MPIconf.comm_depth() );
  
}


void Mesh::read_from_recovery(const MPI_Topology &MPIconf){
	read_array_from_recovery(".//Recovery//Fields//FieldsE",fieldE,MPIconf);	
	read_array_from_recovery(".//Recovery//Fields//FieldsB",fieldB,MPIconf);	
	
	if(DAMP_FIELDS == PML){
	//	read_array_from_recovery(".//Recovery//Fields//FieldsEPML",fieldEp,MPIconf);	
//		read_array_from_recovery(".//Recovery//Fields//FieldsBPML",fieldBp,MPIconf);
  }	
} 

void write_array_to_recovery(const std::string& dataName, const Array3D<double3>& field,const MPI_Topology &MPIconf){
    FILE *fFields;
    char filename[50];
    long info;
    auto size = field.size();
    //auto size_r = field.size_d2();

    sprintf(filename, (dataName + "%03ld").c_str(),MPIconf.rank_line() );

    fFields = fopen(filename, "wb");

    info = size.x();
    fwrite(&info, sizeof(long), 1, fFields);
    info = size.y();
    fwrite(&info, sizeof(long), 1, fFields);
    info = size.z();
    fwrite(&info, sizeof(long), 1, fFields);

    fwrite(&field.data(0), 3*size.x()*size.y()*size.z()*sizeof(double), 1, fFields);
    
  fclose(fFields);
  
}

void Mesh::write_recovery(const MPI_Topology &MPIconf) const {
  
  write_array_to_recovery(".//Recovery//Fields//FieldsE", fieldE,MPIconf);
  write_array_to_recovery(".//Recovery//Fields//FieldsB", fieldB,MPIconf);
  
  if(DAMP_FIELDS == PML){
    //write_array_to_recovery(".//Recovery//Fields//FieldsEPML", fieldEp,MPIconf);
   // write_array_to_recovery(".//Recovery//Fields//FieldsBPML", fieldBp,MPIconf);
  }
}
