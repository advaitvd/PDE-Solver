#pragma once
#include <iostream>
#include <cmath>
using namespace std;

template<class T>
class Matrix{
private:
    T **M;
    int num_rows, num_cols;
public:
    Matrix();
    Matrix(int, int);
    int get_ncols();
    int get_nrows();
    T** get_M();
    void set_element(T,int,int);
    T get_element(int,int);
    void round_zeros();
    Matrix<T> operator+(Matrix<T>);
    Matrix<T> operator-(Matrix<T>);
    Matrix<T> operator*(Matrix<T>);
    Matrix<T> operator+(T);
    Matrix<T> operator-(T);
    Matrix<T> operator*(T);
    Matrix<T>& operator=(const Matrix<T>&);
    
    template<class U>
    friend Matrix<T> operator+(T,Matrix<T>);
    
    template<class U>
    friend Matrix<T> operator-(T,Matrix<T>);
    
    template<class U>
    friend Matrix<T> operator*(T,Matrix<T>);
    
    template<class U>
    friend ostream& operator<<(ostream&, Matrix<T>&);
    
    template<class U>
    friend istream& operator>>(istream&, Matrix<T>&);
    // Define other operators as per requirement.
};

template<class T>
Matrix<T>::Matrix():M(NULL),num_rows(0),num_cols(0){
}

template<class T>
Matrix<T>::Matrix(int rows,int cols):num_rows(rows),num_cols(cols){
	this->M=new T*[rows];
	for(int i=0;i<rows;i++){
		this->M[i]=new T[cols];
	}
}

template<class T>
int Matrix<T>::get_ncols(){
	return num_cols;
}

template<class T>
void Matrix<T>::round_zeros(){
	for(int i=0;i<num_rows;i++){
		for(int j=0;j<num_cols;j++){
			if(fabs(M[i][j])<1e-7){
				M[i][j]=0.0;
			}
		}
	}
}

template<class T>
int Matrix<T>::get_nrows(){
	return num_rows;
}

template<class T>
T** Matrix<T>::get_M(){
	return this->M;
}

template<class T>
void Matrix<T>::set_element(T val,int i,int j){
	this->M[i][j]=val;
}

template<class T>
T Matrix<T>::get_element(int i,int j){
	return this->M[i][j];
}

template<class T>
Matrix<T> Matrix<T>::operator+(Matrix<T> mat){
	if(mat.num_rows!=this->num_rows || mat.num_cols!=this->num_cols){
		cout<<"Invalid"<<endl;
		exit(0);
	}else{
		Matrix<T> res(mat.num_rows,mat.num_cols);
		for(int i=0;i<mat.num_rows;i++){
			for(int j=0;j<mat.num_cols;j++){
				res.M[i][j]=mat.M[i][j]+this->M[i][j];
			}
		}
		return res;
	}
}

template<class T>
Matrix<T> Matrix<T>::operator-(Matrix<T> mat){
	if(mat.num_rows!=this->num_rows || mat.num_cols!=this->num_cols){
		cout<<"Invalid"<<endl;
		exit(0);
	}else{
		Matrix<T> res(mat.num_rows,mat.num_cols);
		for(int i=0;i<mat.num_rows;i++){
			for(int j=0;j<mat.num_cols;j++){
				res.M[i][j]=this->M[i][j]-mat.M[i][j];
			}
		}
		return res;
	}
}

template<class T>
Matrix<T> Matrix<T>::operator*(Matrix<T> mat){
	if(this->num_cols!=mat.num_rows){
		cout<<"Invalid"<<endl;
		exit(0);
	}else{
		Matrix<T> res(this->num_rows,mat.num_cols);
		for(int i=0;i<this->num_rows;i++){
			for(int k=0;k<mat.num_cols;k++){
				T temp=0;
				for(int j=0;j<this->num_cols;j++){
					temp=temp+this->M[i][j]*mat.M[j][k];
				}
				res.M[i][k]=temp;
			}
		}
		return res;
	}
}

template<class T>
Matrix<T> Matrix<T>::operator+(T num){
	Matrix<T> res=Matrix<T>(this->num_rows,this->num_cols);
	for(int i=0;i<this->num_rows;i++){
		for(int j=0;j<this->num_cols;j++){
			res.M[i][j]=this->M[i][j]+num;
		}
	}
	return res;
}

template<class T>
Matrix<T> Matrix<T>::operator-(T num){
	Matrix<T> res=Matrix<T>(this->num_rows,this->num_cols);
	for(int i=0;i<this->num_rows;i++){
		for(int j=0;j<this->num_cols;j++){
			res.M[i][j]=this->M[i][j]-num;
		}
	}
	return res;
}

template<class T>
Matrix<T> Matrix<T>::operator*(T num){
	Matrix<T> res=Matrix<T>(this->num_rows,this->num_cols);
	for(int i=0;i<this->num_rows;i++){
		for(int j=0;j<this->num_cols;j++){
			res.M[i][j]=this->M[i][j]*num;
		}
	}
	return res;
}

template<class T>
Matrix<T> operator+(T num,Matrix<T> mat){
	Matrix<T> res=Matrix<T>(mat.num_rows,mat.num_cols);
	for(int i=0;i<mat.num_rows;i++){
		for(int j=0;j<mat.num_cols;j++){
			res.M[i][j]=num+mat.M[i][j];
		}
	}
	return res;
}

template<class T>
Matrix<T> operator-(T num,Matrix<T> mat){
	Matrix<T> res=Matrix<T>(mat.get_nrows(),mat.get_ncols());
	for(int i=0;i<mat.get_nrows();i++){
		for(int j=0;j<mat.get_ncols();j++){
			res.set_element(num-mat.get_element(i,j),i,j);
		}
	}
	return res;
}

template<class T>
Matrix<T> operator*(T num,Matrix<T> mat){
	Matrix<T> res=Matrix<T>(mat.num_rows,mat.num_cols);
	for(int i=0;i<mat.num_rows;i++){
		for(int j=0;j<mat.num_cols;j++){
			res.M[i][j]=num*(mat.M[i][j]);
		}
	}
	return res;
}

template<class T>
istream& operator>>(istream& in, Matrix<T>& mat){
	for(int i=0;i<mat.num_rows;i++){
		for(int j=0;j<mat.num_cols;j++){
			in>>mat.M[i][j];
		}
	}
	return in;
}

template<class T>
ostream& operator<<(ostream& out, Matrix<T>& mat){
	T** m=mat.get_M();
	for(int i=0;i<mat.get_nrows();i++){
		for(int j=0;j<mat.get_ncols();j++){
			out<<m[i][j]<<" ";
		}
		out<<endl;
	}
	return out;
}

template<class T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& mat){
	this->num_rows=mat.num_rows;
	this->num_cols=mat.num_cols;
	this->M=new T*[mat.num_rows];
	for(int i=0;i<mat.num_rows;i++){
		this->M[i]=new T[mat.num_cols];
		for(int j=0;j<mat.num_cols;j++){
			this->M[i][j]=mat.M[i][j];
		}
	}
	return *this;
}
