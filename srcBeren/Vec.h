#ifndef VEC_H_
#define VEC_H_
#include "defines.h"
#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>
using  ulog = unsigned long;
// *** vec2 ****
template <typename T>
struct vec2 {
    vec2 (const T u, const T v) : d{u,v} {}
    vec2 (const T a[2]) : d{a[0],a[1]} {}
    vec2 (): d{0,0} {}
    T& operator()(int i) {return d[i];}
    const T& operator()(int i) const {return d[i];}
    vec2<T>& operator=(double s) {d[0]=s;d[1]=s;return (*this);}
    vec2<T>& operator+=(vec2<T> o) {d[0]+=o(0);d[1]+=o(1);return(*this);}
    vec2<T>& operator-=(vec2<T> o) {d[0]-=o(0);d[1]-=o(1);return(*this);}
    vec2<T> operator/(double s) {vec2<T>o; o(0)=d[0]/s;o(1)=d[1]/s;return o;}
    vec2<T> operator/=(double s) {d[0]/=s;d[1]/=s;return (*this);}

    //dot product of two vectors
    friend T dot(const vec2<T> &v1, const vec2<T> &v2) {
        T s=0;  for (int i=0;i<2;i++) s+=v1(i)*v2(i);
        return s;   }

    //vector magnitude
    friend T mag(const vec2<T> &v) {return sqrt(dot(v,v));}

    //unit vector
    friend vec2<T> unit(const vec2<T> &v) {return vec2(v)/mag(v);}

    T& x() {return d[0];}
    T x() const {return d[0];}
    T& y() {return d[1];}
    T y() const {return d[1];}
protected:
    T d[2];
};

//vec2-vec2 operations
template<typename T>    //addition of two vec3s
vec2<T> operator+(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T> (a(0)+b(0),a(1)+b(1));   }
template<typename T>    //subtraction of two vec2s
vec2<T> operator-(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T> (a(0)-b(0),a(1)-b(1));   }
template<typename T>    //element-wise multiplication of two vec2s
vec2<T> operator*(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T> (a(0)*b(0),a(1)*b(1));   }
template<typename T>    //element wise division of two vec3s
vec2<T> operator/(const vec2<T>& a, const vec2<T>& b) {
    return vec2<T> (a(0)/b(0),a(1)/b(1));   }

//vec2 - scalar operations
template<typename T>        //scalar multiplication
vec2<T> operator*(const vec2<T> &a, T s) {
    return vec2<T>(a(0)*s, a(1)*s);}
template<typename T>        //scalar multiplication 2
vec2<T> operator*(T s,const vec2<T> &a) {
    return vec2<T>(a(0)*s, a(1)*s);}

//output
template<typename T>    //ostream output
std::ostream& operator<<(std::ostream &out, vec2<T>& v) {
    out<<v(0)<<" "<<v(1)<<" 0"; //paraview does not support 2-component arrays
    return out;
}

using double2 = vec2<double>;
using long2 = vec2<long>;
using int2 = vec2<int>;

template <typename T>
struct vec3 {
    vec3 (const T u, const T v, const T w) : d{u,v,w} {}
    vec3 (const T a[3]) : d{a[0],a[1],a[2]} {}
    vec3 (): d{0,0,0} {}

    T& operator()(int i) {return d[i];}
    T operator()(int i) const {return d[i];}
    vec3<T>& operator=(double s) {d[0]=s;d[1]=s;d[2]=s;return (*this);}
    vec3<T>& operator+=(vec3<T> v) {d[0]+=v(0);d[1]+=v(1);d[2]+=v(2);return(*this);}
    vec3<T>& operator-=(vec3<T> v) {d[0]-=v(0);d[1]-=v(1);d[2]-=v(2);return(*this);}
    vec3<T> operator/(double s) {vec3<T> v; v(0) = d[0] / s; v(1) = d[1] / s; v(2) = d[2] / s; return v;}
    vec3<T> operator/=(double s) {d[0]/=s;d[1]/=s;d[2]/=s;return (*this);}

    //dot product of two vectors
    friend T dot(const vec3<T> &v1, const vec3<T> &v2) {
        T s=0;  for (int i=0;i<3;i++) s+=v1(i)*v2(i);
        return s;   }

    //vector magnitude
    friend T mag(const vec3<T> &v) {return sqrt(dot(v,v));}

    //unit vector
    friend vec3<T> unit(const vec3<T> &v) {if (mag(v)>0) return vec3(v)/mag(v); else return {0,0,0};}

    //cross product
    friend vec3<T> cross(const vec3<T> &a, const vec3<T> &b) {
        return {a(1)*b(2)-a(2)*b(1), a(2)*b(0)-a(0)*b(2), a(0)*b(1)-a(1)*b(0)};
    }

    T& x() {return d[0];}
    T x() const {return d[0];}
    T& y() {return d[1];}
    T y() const {return d[1];}
    T& z() {return d[2];}
    T z() const {return d[2];}
    T& r() {return d[1];}
    T r() const {return d[1];}
    T& p() {return d[2];}
    T p() const {return d[2];}
protected:
    T d[3];
};

//vec3-vec3 operations
template<typename T>    //addition of two vec3s
vec3<T> operator+(const vec3<T>& a, const vec3<T>& b) {
    return vec3<T> (a(0)+b(0),a(1)+b(1),a(2)+b(2)); }
template<typename T>    //subtraction of two vec3s
vec3<T> operator-(const vec3<T>& a, const vec3<T>& b) {
    return vec3<T> (a(0)-b(0),a(1)-b(1),a(2)-b(2)); }
template<typename T>    //element-wise multiplication of two vec3s
vec3<T> operator*(const vec3<T>& a, const vec3<T>& b) {
    return vec3<T> (a(0)*b(0),a(1)*b(1),a(2)*b(2)); }
template<typename T>    //element wise division of two vec3s
vec3<T> operator/(const vec3<T>& a, const vec3<T>& b) {
    return vec3<T> (a(0)/b(0),a(1)/b(1),a(2)/b(2)); }

//vec3 - scalar operations
template<typename T>        //scalar multiplication
vec3<T> operator*(const vec3<T> &a, T s) {
    return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}
template<typename T>        //scalar multiplication 2
vec3<T> operator*(T s,const vec3<T> &a) {
    return vec3<T>(a(0)*s, a(1)*s, a(2)*s);}

//output
template<typename T>    //ostream output
std::ostream& operator<<(std::ostream &out, const vec3<T>& v) {
    out<<v(0)<<" "<<v(1)<<" "<<v(2);
    return out;
}

using double3 = vec3<double>;
using int3 = vec3<int>;
using long3 = vec3<long>;
template <typename T> 
struct Array3D{

    Array3D(long n1, long n2, long n3){
       allocate(n1,n2,n3);
    }
    Array3D(const long3& nn){
       allocate(nn(0),nn(1),nn(2));
    }
    Array3D(Array3D &&other): _data{other._data}, 
                  _size1{other._size1},_size2{other._size2},_size3{other._size3} {
        other._data = nullptr;
        other._size1 = other._size2 = other._size3 =  0; 
    }
    Array3D(const Array3D &other): _size1{other._size1},_size2{other._size2},_size3{other._size3} { 
       	allocate( _size1,_size2,_size3 );
    	for(auto i = 0; i < capacity(); ++i){
            _data[i] = other._data[i] ;
        }
    }
    Array3D(){
        _data = nullptr;
        _size1 = _size2 =_size3 = 0.;
    }
    
    void allocate(long n1, long n2, long n3){
        _data = new T[ n1 * n2 * n3];
        _size1 = n1;
        _size2 = n2;
        _size3 = n3;
    }
    
    void clear(){
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = 0.;
        }
    }

    void free(){
        if (_data != nullptr)
            delete[] _data;
        _size1 = _size2 =_size3 = 0.;
    }

    ~Array3D(){
       free();
    }
    
    Array3D<T>& operator=(T s) {
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = s ;
        }
        return (*this);
    }
    
    Array3D<T>& operator=(const Array3D<T>& array3D) {
    	// Self-assignment check
    	if (this == &array3D)
        	return *this;
        // Checking the conformity of dimensions
       	assert( capacity() == array3D.capacity() );
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = array3D._data[i] ;
        }
        return (*this);
    }   
    Array3D<T>& operator=(Array3D<T> &&array3D) {
    	// Self-assignment check
    	if (this == &array3D)
        	return *this;
        // Checking the conformity of dimensions
       	assert( capacity() == array3D.capacity() );
		delete[] _data;
       	_data = array3D._data;
        _size1 = array3D._size1;
        _size2 = array3D._size2;
        _size3 = array3D._size3;
        array3D._data = nullptr;
        array3D._size1 = array3D._size2 = array3D._size3 =  0; 
        return (*this);
    }   

    Array3D<T>& operator-=(const Array3D<T>& array3D) {
        assert( capacity() == array3D.capacity() );
        for(auto i = 0; i < capacity(); ++i){
            _data[i] -= array3D._data[i] ;
        }
        return (*this);
    }
    
    Array3D<T>& operator+=(const Array3D<T>& array3D) {
        assert( capacity() == array3D.capacity() );
        for(auto i = 0; i < capacity(); ++i){
            _data[i] += array3D._data[i] ;
        }
        return (*this);
    }          
    
    T& operator() (long i, long j, long k) {
        #if DEBUG > 0 
            if ( capacity() <= i * _size2 * _size3 + j * _size3 + k || i * _size2 * _size3 + j * _size3 + k < 0 ){
                std::string msg = "IndexError3";
                std::cout << msg << ": "<< long3(i,j,k)<< "\n";
                throw msg;
            };
        #endif
    	return _data[i * _size2 * _size3 + j * _size3 + k ];
    }

    const T& operator() (long i, long j, long k) const{ 
        #if DEBUG > 0 
            if ( capacity() <= i * _size2 * _size3 + j * _size3 + k || i * _size2 * _size3 + j * _size3 + k < 0 ){
                std::string msg = "IndexError3";
                std::cout << msg << ": "<< long3(i,j,k) << "\n";
                throw msg;
            };
        #endif
    	return _data[i * _size2 * _size3 + j * _size3 + k];
    }
    T& data(long i) {
        #if DEBUG > 0 
            if ( capacity() <= i || i < 0 ){
                std::string msg = "IndexError3";
                std::cout << msg << ": " << i << "\n";
                throw msg;
            };
        #endif
        return _data[i];
    }
    const T& data(long i) const {
        #if DEBUG > 0 
            if ( capacity() <= i || i < 0){
                std::string msg = "IndexError3";
                std::cout << msg << ": " << i << "\n";
                throw msg;
            };
        #endif
        return _data[i];
    }
  
    long3 size() const{
        return long3(_size1,_size2,_size3);
    }
    long capacity() const{
        return _size1*_size2*_size3;
    }

protected:
    T* _data;
    long _size1, _size2, _size3;
};

template <typename T> 
struct Array2D{
    Array2D(long n1, long n2){
       allocate(n1,n2);
    }
    Array2D(long2 nn){
       allocate(nn(0),nn(1));
    }
    Array2D(Array2D &&other): _data{other._data}, _size1{other._size1}, _size2{other._size2}{
        other._data = nullptr;
        other._size1 = other._size2 = 0; 
    }
    Array2D(){
        _data = nullptr;
        _size1 = _size2 = 0; 
    }
    
    void allocate(long n1, long n2){
        _data = new T[ n1 * n2];
        _size1 = n1;
        _size2 = n2;
    }
    
    void clear(){
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = 0.;
        }
    }
    void free(){
        if (_data != nullptr)
            delete[] _data;
        _size1 = _size2 = 0.;
    }
    
    ~Array2D(){
       free();
    }
    Array2D<T>& operator=(T s) {
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = s ;
        }
        return (*this);
    }

    Array2D<T>& operator=(const Array2D<T>& array2D) {
        assert( capacity() == array2D.capacity() && "Array 2D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = array2D._data[i] ;
        }
        return (*this);
    }
    Array2D<T>& operator=(const Array2D<T>&& other) {
        //assert( capacity() == other.capacity() );
        
        if(&other == this) return (*this);

        free();
        _size1 = other._size1;
        _size2 = other._size2;
        _data = other._data;
        other._data = nullptr;
        other._size1 = other._size2 = 0;
        
        return (*this);
    }

    Array2D<T>& operator-=(const Array2D<T>& array2D) {
        assert( capacity() == array2D.capacity()  && "Array 2D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] -= array2D._data[i] ;
        }
        return (*this);
    }
    
    Array2D<T>& operator+=(const Array2D<T>& array2D) {
        assert( capacity() == array2D.capacity()  && "Array 2D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] += array2D._data[i] ;
        }
        return (*this);
    }

    T& operator() (long i, long j) {
        #if DEBUG > 0 
            if ( capacity() <= i * _size2  + j ||  i * _size2  + j < 0 ){
                std::string msg = "IndexError2";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif
	return _data[i * _size2 + j];
    }

    const T& operator() (long i, long j) const{	
        #if DEBUG > 0 
            if ( capacity() <= i * _size2  + j ||  i * _size2  + j < 0 ){
                std::string msg = "IndexError2";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif	
        return _data[i * _size2 + j];
    }
    
    T& data(long i) {
        #if DEBUG > 0 
            if ( capacity() <= i || i < 0 ){
                std::string msg = "IndexError2";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif
        return _data[i];
    }
    const T& data(long i) const {
        #if DEBUG > 0 
            if ( capacity() <= i || i < 0){
                std::string msg = "IndexError2";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif
        return _data[i];
    }
  
    long2 size() const{
        return long2(_size1,_size2);
    }
    long capacity() const{
        return _size1*_size2;
    }
    T sum_d2(long i) const{
        T sum = 0;
        for( long t =  0; t < _size2; ++t){
          sum += _data[i * _size2 + t];
        }
        return sum;
    }

protected:
    T* _data;
    long _size1,_size2;
};

template <typename T> 
struct Array1D{
    
    Array1D(long size){
       allocate(size);
    }
    Array1D(Array1D &&other): _data{other._data}, _size{other._size} {
        other._data = nullptr;
        other._size = 0; 
    }
    Array1D(){
        _size = 0.;
    }
    
    void allocate(long sizeDim){
        _data = new T[ sizeDim];
        _size =  sizeDim;
    }
    
    void clear(){
        for(auto i = 0; i <_size; ++i){
            _data[i] = 0.;
        }
    }
    void free(){
        if (_data != nullptr)
            delete[] _data;
        _size = 0.;
    }
    
    Array1D<T>& operator=(T s) {
        for(auto i = 0; i < _size; ++i){
            _data[i] = s ;
        }
        return (*this);
    }

    Array1D<T>& operator=(const Array1D<T>& other) {
        assert( capacity() == other.capacity() && "Array 1D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] = other._data[i] ;
        }
        return (*this);
    }

    Array1D<T>& operator-=(const Array1D<T>& other) {
        assert( capacity() == other.capacity()  && "Array 1D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] -= other._data[i] ;
        }
        return (*this);
    }
    
    Array1D<T>& operator+=(const Array1D<T>& other) {
        assert( capacity() == other.capacity()  && "Array 1D sizes do not match!\n");
        for(auto i = 0; i < capacity(); ++i){
            _data[i] += other._data[i] ;
        }
        return (*this);
    }
    Array1D<T>& operator=(const Array1D<T>&& other) {
        
        if( &other == this) return (*this);

        free();
        _size = other._size;
        _data = other._data;
        other._data = nullptr;
        other._size = 0;
        return (*this);
    }
    ~Array1D(){
        free();
    }
    
    T& operator() (long i) {
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexError";
                std::cout << msg << ": index = " << i <<". capacity = "<< capacity() << "\n";
                throw msg;
            };
        #endif      
        return _data[i];
    }

    const T& operator() (long i) const{ 
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexError";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif          
        return _data[i];
    }
    long size() const{
        return _size;
    }
    long capacity() const{
        return _size;
    }
    T& data(long i) {
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexError";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif  
        return _data[i];
    }
    const T& data(long i) const {
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexError";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif  
        return _data[i];
    }  

protected:
    T *_data;
    long _size;  
};
/*
template <typename T> 
struct Array{

    Array(long maxSize){
        allocate(maxSize);
        _size = 0;
    }
    
    Array(){
        _data = nullptr;
        _maxSize = 0;
        _size = 0;
    };
    Array(Array &&arr): _data(arr._data), _size(arr._size), _maxSize(arr._maxSize) {
        arr._data = nullptr;
        arr._size = 0;
        arr._maxSize = 0;
    }
    ~Array(){
        free();
    }
        
    void allocate(long maxSize){
        _data = new T[maxSize];
        _maxSize = maxSize;
    }

    void clear(){
        _size = 0;
    }   
    
    void free(){
        if( _data!= nullptr)
            delete[] _data;
        _size = 0;
    }   
    

    void push_back(const T& elem){
        _data[_size] = elem;
        ++_size;
        #if DEBUG > 0 
            if ( _size >= capacity() ){
                std::string msg = "AddError";
                std::cout << msg << ": " << _size << " " << _maxSize << "\n";
                throw msg;
            };
        #endif  
    }
    
    void del(long k){
        #if DEBUG > 0 
            if ( _size < k ||  k < 0 ){
                std::string msg = "DelError";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif  
        --_size;
        _data[k] = _data[_size];
    }
    long size() const{
        return _size;
    }
    void resize(long newSize){
        _size = newSize;
    }
    long capacity() const{
        return _maxSize;
    }
    T& operator() (long i) {
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexErrorP";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif          
        return _data[i];
    }

    const T& operator() (long i) const{
        #if DEBUG > 0 
            if ( capacity() <= i ||  i < 0 ){
                std::string msg = "IndexErrorP";
                std::cout << msg << "\n";
                throw msg;
            };
        #endif  
        return _data[i];
    }
    T back(){
        return _data[_size-1];
    }
    void pop_back(){
        --_size;
    }
        T* _data;
        long _size, _maxSize;
};
*/
template <class T>
struct Array : std::vector<T> {
    using std::vector<T>::reserve;
    //using std::vector<T>::size;
    using std::vector<T>::pop_back;
    using std::vector<T>::erase;
    Array(long maxSize){
        reserve(maxSize);
        //_size = 0;
    }    
    T& operator() (long i) {
        return (*this)[i];
    }

    const T& operator() (long i) const{
        return (*this)[i];
    }
    void del(long k){
        (*this)[k] = (*this)[ size()-1];
        pop_back();
    }
    long size() const{
        return long(std::vector<T>::size() );
    }

};
#endif 
