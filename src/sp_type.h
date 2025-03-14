#ifndef SP_TYPE
#define SP_TYPE
#include "sp_misc.h"

template <class T, size_t m_dat>
class Array1{
	T dat[m_dat];
 public:
	Array1(){}
	Array1(T v){for(int m=0; m<m_dat; m++) dat[m] = v;}
	Array1(const Array1 &from){for(int m=0; m<m_dat; m++) dat[m] = from.dat[m];}
	Array1(const T *xt){for(int m=0; m<m_dat; m++) dat[m] = xt[m];}
	~Array1(){}
//
	void setx(T vs[]){for(int m=0; m<m_dat; m++) dat[m] = vs[m];}
	void fill(T v){for(int m=0; m<m_dat; m++) dat[m] = v;}
//	operator T*(){return dat;}
	T *getdat(){return dat;}
	T *getx(){return dat;}
	T &operator[](size_t it){return dat[it];}
	Array1 operator+(const Array1 &v){
		T xt[m_dat];
		for(int m=0; m<m_dat; m++) xt[m] = dat[m] + v.dat[m];
		return Array1(xt);
	}
	Array1 operator-(const Array1 &v){
		T xt[m_dat];
		for(int m=0; m<m_dat; m++) xt[m] = dat[m] - v.dat[m];
		return Array1(xt);
	}
	Array1 operator*(const T v){
		T xt[m_dat];
		for(int m=0; m<m_dat; m++) xt[m] = dat[m] * v;
		return Array1(xt);
	}
	Array1 operator/(const T v){
		T xt[m_dat];
		for(int m=0; m<m_dat; m++) xt[m] = dat[m] / v;
		return Array1(xt);
	}
	T normalize(){
		double s = 0;
		for(int m=0; m<m_dat; m++) s += dat[m] * dat[m];
		s = max(1e-99, sqrt(s));
		for(int m=0; m<m_dat; m++) dat[m] /= s;
		return s;
	}
	Array1 &operator=(const Array1 &from){
		for(int m=0; m<m_dat; m++) dat[m] = from.dat[m]; return(*this);}
	Array1 &operator+=(const Array1 &from){
		for(int m=0; m<m_dat; m++) dat[m] += from.dat[m]; return(*this);}
	Array1 &operator-=(const Array1 &from) {
		for(int m=0; m<m_dat; m++) dat[m] -= from.dat[m]; return(*this);}
	Array1 &operator*=(const T v) {
		for(int m=0; m<m_dat; m++) dat[m] *= v; return (*this);}
	Array1 &operator/=(const T v) {
		for(int m=0; m<m_dat; m++) dat[m] /= v; return (*this);}
	T dot_product(const Array1 &v){
		T s = 0;
		for(int m=0; m<m_dat; m++) s += dat[m] * v.dat[m];
		return s;
	}
	T distance2(Array1 const &v){
		T x_, s_ = 0;
		//#pragma omp simd lastprivate(s_) //this is not a good optimization 
		for(int m=0; m<m_dat; m++){
			x_ = dat[m] - v.dat[m]; s_ += x_ * x_;
		}
		return s_;
	}
//
	friend inline ostream &operator<<(ostream &os, Array1 const &v){
		for(int m=0; m<m_dat; m++) os<<' '<<v.dat[m];
		os<<endl;
		return os;
	}
	friend inline int operator<(Array1 const &v1, Array1 const &v2){
		for(int m=0; m<m_dat; m++){
			if(v1.dat[m] < v2.dat[m]) return 1;
			else if(v1.dat[m] > v2.dat[m]) return 0;
		}
		return 0;
	}
};

typedef Array1<double, 3> Xvec;

template <typename T>
class TwoDVector {
public:
  TwoDVector() : rows_(0), cols_(0) {}

  void resize(size_t rows, size_t cols) {
    rows_ = rows;
    cols_ = cols;
    data_.resize(rows_ * cols_);
    index_.resize(rows_);
    for (size_t i = 0; i < rows_; ++i) {
      index_[i] = &data_[i * cols_];
    }
  }

  void reset() {
    std::fill(data_.begin(), data_.end(), T());
  }

  void fill(T value) {
	std::fill(data_.begin(), data_.end(), value);
  }

  T* operator[](size_t i) {
    return index_[i];
  }

  const T* operator[](size_t i) const {
    return index_[i];
  }

private:
  size_t rows_;
  size_t cols_;
  std::vector<T> data_;
  std::vector<T*> index_;
};


#endif
