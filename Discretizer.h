#pragma once
#include "Matrix.h"
#include "AD.h"
#include "charList.h"
#include<cmath>
#include<unordered_map>
#include "NonLinearSolver.h"
using namespace std;

template<class T>
class Discretizer
{
	private:
		int nx,ny,iters;
		T a,b,c,d;
		T f1,f2,g1,g2;
		string A,B,C,D,E;
		T dx,dy;
		vector<AD<T>> X,Y;
		Matrix<AD<T>>U;
		string linear_solver,non_linear_solver;
	public:
		Discretizer(T,T,T,T,T,T,T,T,T,string,string,string,string,string,int,int,int,string,string);
		Matrix<T> solve();
};

template<class T>
Discretizer<T>::Discretizer(T initial_guess,T a,T b,T c,T d,T f1,T f2,T g1,T g2,string A,string B,string C,string D,string E,int Nx,int Ny,int iters,string non_linear_solver,string linear_solver)
	:nx(Nx),ny(Ny),iters(iters),a(a),b(b),c(c),d(d),f1(f1),f2(f2),g1(g1),g2(g2),A(A),B(B),C(C),D(D),E(E),linear_solver(linear_solver),non_linear_solver(non_linear_solver){
	dx=(b-a)/(T)(Nx-1);
	dy=(d-c)/(T)(Ny-1);
	X=vector<AD<T>>(Nx);
	Y=vector<AD<T>>(Ny);
	
	for(int i=0;i<Nx;i++){
		X[i]=AD<T>(a+i*dx,0);
		X[i].setIndVar(3);
	}
	for(int i=0;i<Ny;i++){
		Y[i]=AD<T>(c+i*dy,1);
		Y[i].setIndVar(3);
	}
	
	U=Matrix<AD<T>>(Nx,Ny);
	int varCount=(Nx-2)*(Ny-2);
	int id=0;
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			AD<T> u;
			if(i!=0 && i!=Nx-1 && j!=0 && j!=Ny-1){
				AD<T>temp(initial_guess,id);
				u=temp;
				id++;
				u.setIndVar(varCount);
			}else if(i==0){
				u.set_f(f1);
			}else if(j==0){
				u.set_f(g1);
			}else if(i==Nx-1){
				u.set_f(f2);
			}else if(j==Ny-1){
				u.set_f(g2);
			}
			U.set_element(u,i,j);
		}
	}
	cout<<"Successfully Discretized!"<<endl;
}

template<class T>
Matrix<T> Discretizer<T>::solve(){
	NonLinearSolver<T>* ptr=nullptr;
	if(this->non_linear_solver=="Newton"){
		ptr=new Newton<T>();
	}else if(this->non_linear_solver=="Broyden"){
		ptr=new Broyden<T>();
	}else{
		cout<<"Enter a valid Non-linear solver!"<<endl;
		exit(0);
	}
	Matrix<T> res=ptr->solve(nx,ny,dx,dy,iters,linear_solver,U,X,Y,A,B,C,D,E);
	
	delete ptr;
	return res;
}

