#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include "Matrix.h"
#include "AD.h"
#include "LinearSolver.h"
#include "charList.h"
#include<unordered_map>
using namespace std;

template<class T>
class NonLinearSolver{
	public:
		virtual Matrix<T> solve(int nx,int ny,T dx,T dy,int iters,string linear_solver,Matrix<AD<T>> &U,vector<AD<T>> &X,vector<AD<T>> &Y,string A,string B,string C,string D,string E)=0;
		virtual ~NonLinearSolver(){}
};

template<class T>
class Newton: public NonLinearSolver<T>{
	public:
		virtual Matrix<T> solve(int nx,int ny,T dx,T dy,int iters,string linear_solver,Matrix<AD<T>> &U,vector<AD<T>> &X,vector<AD<T>> &Y,string A,string B,string C,string D,string E);
};

template<class T>
class Broyden: public NonLinearSolver<T>{
	public:
		virtual Matrix<T> solve(int nx,int ny,T dx,T dy,int iters,string linear_solver,Matrix<AD<T>> &U,vector<AD<T>> &X,vector<AD<T>> &Y,string A,string B,string C,string D,string E);
};

template<class T>
Matrix<T> Broyden<T>::solve(int nx,int ny,T dx,T dy,int iters,string linear_solver,Matrix<AD<T>> &U,vector<AD<T>> &X,vector<AD<T>> &Y,string A,string B,string C,string D,string E){
	Matrix<T> res(U.get_nrows(),U.get_ncols());
	Matrix<T> J;
	AD<T> F[(nx-2)*(ny-2)];
	T A_M,B_M,C_M,D_M,E_M;
	unordered_map<string,AD<T>> vars;
	vector<T> XX((nx-2)*(ny-2)),XX_prev((nx-2)*(ny-2)),bb((nx-2)*(ny-2)),bb_prev((nx-2)*(ny-2));
	
	for(int tt=0;tt<iters;tt++){
		for(int i=1;i<nx-1;i++){
			for(int j=1;j<ny-1;j++){
				charList<T> A_expr(A),B_expr(B),C_expr(C),D_expr(D),E_expr(E);
				vars["x"]=X[i];
				vars["y"]=Y[j];
				vars["u"]=AD<T>(U.get_element(i,j).getf(),2);
				XX[(i-1)*(ny-2)+j-1]=(U.get_element(i,j).getf());
				vars["u"].setIndVar(3);
				
				A_expr.infix_to_postfix(vars);
				B_expr.infix_to_postfix(vars);
				C_expr.infix_to_postfix(vars);
				D_expr.infix_to_postfix(vars);
				E_expr.infix_to_postfix(vars);
				
				A_M=A_expr.evaluate(3).getf();
				B_M=B_expr.evaluate(3).getf();
				C_M=C_expr.evaluate(3).getf();
				D_M=D_expr.evaluate(3).getf();
				E_M=E_expr.evaluate(3).getf();
				
				AD<T> tmp1,tmp2;
				if(i+1==nx-1){
					tmp1=A_M*(U.get_element(i+1,j).getf()-2.0*U.get_element(i,j)+U.get_element(i-1,j))/(2.0*pow(dx,2.0));
					tmp1=tmp1+(U.get_element(i+1,j).getf()-U.get_element(i-1,j))*C_M/(2.0*dx);
				}else if(i-1==0){
					tmp1=A_M*(U.get_element(i+1,j)-2.0*U.get_element(i,j)+U.get_element(i-1,j).getf())/(2.0*pow(dx,2.0));
					tmp1=tmp1+(U.get_element(i+1,j)-U.get_element(i-1,j).getf())*C_M/(2.0*dx);
				}else{
					tmp1=A_M*(U.get_element(i+1,j)-2.0*U.get_element(i,j)+U.get_element(i-1,j))/(2.0*pow(dx,2.0));
					tmp1=tmp1+(U.get_element(i+1,j)-U.get_element(i-1,j))*C_M/(2.0*dx);
				}
				
				if(j+1==ny-1){
					tmp2=(U.get_element(i,j+1).getf()-2.0*U.get_element(i,j)+U.get_element(i,j-1))*B_M/(2.0*pow(dy,2.0));
					tmp2=tmp2+(U.get_element(i,j+1).getf()-U.get_element(i,j-1))*D_M/(2.0*dy);
				}else if(j-1==0){
					tmp2=(U.get_element(i,j+1)-2.0*U.get_element(i,j)+U.get_element(i,j-1).getf())*B_M/(2.0*pow(dy,2.0));
					tmp2=tmp2+(U.get_element(i,j+1)-U.get_element(i,j-1).getf())*D_M/(2.0*dy);
				}else{
					tmp2=(U.get_element(i,j+1)-2.0*U.get_element(i,j)+U.get_element(i,j-1))*B_M/(2.0*pow(dy,2.0));
					tmp2=tmp2+(U.get_element(i,j+1)-U.get_element(i,j-1))*D_M/(2.0*dy);
				}
				
				F[(i-1)*(ny-2)+j-1]=tmp1+tmp2-E_M;
				bb[(i-1)*(ny-2)+j-1]=(F[(i-1)*(ny-2)+j-1].getf());
			}
		}
		if(tt==0){
			
			//find the jacobian only in the first iteration
			J=getJacobian(F,(nx-2)*(ny-2),(nx-2)*(ny-2));
			J=J.inverse();
			Matrix<T> FF(bb.size(),1);
			for(int i=0;i<bb.size();i++){
				FF.set_element(bb[i],i,0);
			}
			
			Matrix<T> yy=J*FF;
			
			for(int i=0;i<XX_prev.size();i++){
				XX_prev[i]=XX[i];
				XX[i]=XX[i]-yy.get_element(i,0);
			}
			bb_prev=vector<T>(bb);
		}else{
			Matrix<T> y(bb_prev.size(),1),s(XX_prev.size(),1);
			for(int i=0;i<bb_prev.size();i++){
				y.set_element(bb[i]-bb_prev[i],i,0);
				s.set_element(XX[i]-XX_prev[i],i,0);
			}
			
			T val=(s.transpose()*s).get_element(0,0);
			cout<<val<<endl;
//			cout<<J.get_nrows()<<"x"<<J.get_ncols()<<endl;
//			cout<<s.get_nrows()<<"x"<<s.get_ncols()<<endl;
			J=J+(y-J*s)*s.transpose()/val;
			
			Matrix<T> FF(bb.size(),1);
			for(int i=0;i<bb.size();i++){
				FF.set_element(bb[i],i,0);
			}
			Matrix<T> yy=J*FF;
			XX=vector<T>(XX_prev.size());
			for(int i=0;i<XX_prev.size();i++){
				XX_prev[i]=XX[i];
				XX[i]=XX[i]-yy.get_element(i,0);
			}
			bb_prev=vector<T>(bb);
		}
		for(int i=1;i<nx-1;i++){
			for(int j=1;j<ny-1;j++){
				AD<T> elem=U.get_element(i,j);
				elem.set_f(XX[(i-1)*(ny-2)+j-1]);
				U.set_element(elem,i,j);
			}
		}
	}
	
	for(int i=0;i<U.get_nrows();i++){
		for(int j=0;j<U.get_ncols();j++){
			res.set_element(U.get_element(i,j).getf(),i,j);
		}
	}

	return res;
}

template<class T>
Matrix<T> Newton<T>::solve(int nx,int ny,T dx,T dy,int iters,string linear_solver,Matrix<AD<T>> &U,vector<AD<T>> &X,vector<AD<T>> &Y,string A,string B,string C,string D,string E){
	LinearSolver<T>* linsol=nullptr;
	if(linear_solver=="LU_Decomposition"){
		linsol=new LU_Decomposition<T>();
	}else if(linear_solver=="TriDiagonal"){
		linsol=new TriDiagonal<T>();
	}else if(linear_solver=="Gauss_Jacobi"){
		linsol=new Gauss_Jacobi<T>();
	}else if(linear_solver=="Gauss_Seidal"){
		linsol=new Gauss_Seidal<T>();
	}else if(linear_solver=="SOR"){
		linsol=new SOR<T>();
	}else if(linear_solver=="GaussElimination"){
		linsol=new GaussElimination<T>();
	}else{
		cout<<"Enter a valid linear solver!"<<endl;
		exit(0);
	}
	Matrix<T> res(U.get_nrows(),U.get_ncols());
	AD<T> F[(nx-2)*(ny-2)];
	T A_M,B_M,C_M,D_M,E_M;
	unordered_map<string,AD<T>> vars;
//	bool flag=true;
	for(int tt=0;tt<iters;tt++){
		vector<T> XX,XX_prev,bb,xx;
		for(int i=1;i<nx-1;i++){
			for(int j=1;j<ny-1;j++){
//				cout<<U.get_element(i,j).id<<" "<<U.get_element(i,j).df.size()<<endl;
				charList<T> A_expr(A),B_expr(B),C_expr(C),D_expr(D),E_expr(E);
				vars["x"]=X[i];
				vars["y"]=Y[j];
				vars["u"]=AD<T>(U.get_element(i,j).getf(),2);
				XX_prev.push_back(U.get_element(i,j).getf());
				vars["u"].setIndVar(3);
				A_expr.infix_to_postfix(vars);
				B_expr.infix_to_postfix(vars);
				C_expr.infix_to_postfix(vars);
				D_expr.infix_to_postfix(vars);
				E_expr.infix_to_postfix(vars);
				
				A_M=A_expr.evaluate(3).getf();
				B_M=B_expr.evaluate(3).getf();
				C_M=C_expr.evaluate(3).getf();
				D_M=D_expr.evaluate(3).getf();
				E_M=E_expr.evaluate(3).getf();
				
//				if(flag){
//					flag=false;
//					cout<<A_M<<" "<<B_M<<" "<<C_M<<" "<<D_M<<" "<<E_M<<endl;
//				}
				
				AD<T> tmp1,tmp2;
				if(i+1==nx-1){
					tmp1=A_M*(U.get_element(i+1,j).getf()-2.0*U.get_element(i,j)+U.get_element(i-1,j))/(2.0*pow(dx,2.0));
					tmp1=tmp1+(U.get_element(i+1,j).getf()-U.get_element(i-1,j))*C_M/(2.0*dx);
				}else if(i-1==0){
					tmp1=A_M*(U.get_element(i+1,j)-2.0*U.get_element(i,j)+U.get_element(i-1,j).getf())/(2.0*pow(dx,2.0));
					tmp1=tmp1+(U.get_element(i+1,j)-U.get_element(i-1,j).getf())*C_M/(2.0*dx);
				}else{
					tmp1=A_M*(U.get_element(i+1,j)-2.0*U.get_element(i,j)+U.get_element(i-1,j))/(2.0*pow(dx,2.0));
					tmp1=tmp1+(U.get_element(i+1,j)-U.get_element(i-1,j))*C_M/(2.0*dx);
				}
				
				if(j+1==ny-1){
					tmp2=(U.get_element(i,j+1).getf()-2.0*U.get_element(i,j)+U.get_element(i,j-1))*B_M/(2.0*pow(dy,2.0));
					tmp2=tmp2+(U.get_element(i,j+1).getf()-U.get_element(i,j-1))*D_M/(2.0*dy);
				}else if(j-1==0){
					tmp2=(U.get_element(i,j+1)-2.0*U.get_element(i,j)+U.get_element(i,j-1).getf())*B_M/(2.0*pow(dy,2.0));
					tmp2=tmp2+(U.get_element(i,j+1)-U.get_element(i,j-1).getf())*D_M/(2.0*dy);
				}else{
					tmp2=(U.get_element(i,j+1)-2.0*U.get_element(i,j)+U.get_element(i,j-1))*B_M/(2.0*pow(dy,2.0));
					tmp2=tmp2+(U.get_element(i,j+1)-U.get_element(i,j-1))*D_M/(2.0*dy);
				}
				
				F[(i-1)*(ny-2)+j-1]=tmp1+tmp2-E_M;
				bb.push_back(F[(i-1)*(ny-2)+j-1].getf());
			}
		}
		Matrix<T> A=getJacobian(F,(nx-2)*(ny-2),(nx-2)*(ny-2));
		//cout<<"Jacobian found successfully!"<<endl;
		//solve the system XX=XX_prev-A_inv*bb;
		//Input for the Nonlinearsolver:
		//A,XX,XX_prev
		vector<T> y=linsol->solve(A,bb);
		XX=vector<T>(y.size());
		for(int q=0;q<y.size();q++){
			XX[q]=XX_prev[q]-y[q];
		}
		for(int i=1;i<nx-1;i++){
			for(int j=1;j<ny-1;j++){
				AD<T> elem=U.get_element(i,j);
				elem.set_f(XX[(i-1)*(ny-2)+j-1]);
				U.set_element(elem,i,j);
			}
		}
	}
	for(int i=0;i<U.get_nrows();i++){
		for(int j=0;j<U.get_ncols();j++){
			res.set_element(U.get_element(i,j).getf(),i,j);
		}
	}
	delete linsol;
	return res;
}

