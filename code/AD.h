#pragma once
#include "Matrix.h"
#include <cmath>
#include <vector>
#include <assert.h>
using namespace std;

template<class T>
class AD
{
	public:
		T f;
		vector<T> df;
		int id;
	public:
		AD();
		AD(T,int);
		AD<T>& operator=(const AD<T>&);
		void setIndVar(int);
		void set_f(T val);
		void set_df_const(int n);
		T getf();
		T getDf(int);
		vector<T> getGradient();
		
		template<class U>
		friend Matrix<T> getJacobian(AD<T>*,int n,int m);
		
		AD<T> operator*(AD<T>);//
		AD<T> operator+(AD<T>);//
		AD<T> operator-(AD<T>);//
		AD<T> operator/(AD<T>);//
		AD<T> operator^(AD<T>);//
		
		AD<T> operator*(T);//
		AD<T> operator+(T);//
		AD<T> operator-(T);//
		AD<T> operator/(T);//
		AD<T> operator^(T);//
		
		template<class U>
		friend AD<T> operator*(T,AD<T>);//
		template<class U>
		friend AD<T> operator+(T,AD<T>);//
		template<class U>
		friend AD<T> operator-(T,AD<T>);//
		template<class U>
		friend AD<T> operator/(T,AD<T>);//
		template<class U>
		friend AD<T> operator^(T,AD<T>);
		
		template<class U>
		friend AD<T> sin(AD<T>);
		
		template<class U>
		friend AD<T> cos(AD<T>);
		
		template<class U>
		friend AD<T> tan(AD<T>);
		
		template<class U>
		friend AD<T> cosec(AD<T>);
		
		template<class U>
		friend AD<T> sec(AD<T>);
		
		template<class U>
		friend AD<T> cot(AD<T>);
		
		template<class U>
		friend AD<T> arcsin(AD<T>);
		
		template<class U>
		friend AD<T> arccos(AD<T>);
		
		template<class U>
		friend AD<T> arctan(AD<T>);
		
		template<class U>
		friend AD<T> sinh(AD<T>);
		
		template<class U>
		friend AD<T> cosh(AD<T>);
		
		template<class U>
		friend AD<T> tanh(AD<T>);
		
		template<class U>
		friend AD<T> log(AD<T>);
		
		template<class U>
		friend AD<T> abs(AD<T>);
		
		template<class U>
		friend AD<T> exp(AD<T>);
};

template<class T>
AD<T>::AD(){
	f=0;
}

template<class T>
AD<T>::AD(T value,int id){
	this->f=value;
	this->id=id;
}

template<class T>
AD<T>& AD<T>::operator=(const AD<T>& rhs){
	this->f=rhs.f;
	this->df=vector<T>(rhs.df);
	this->id=rhs.id;
	return *this;
}

template<class T>
void AD<T>::setIndVar(int varCount){
	this->df=vector<T>(varCount,0.0);
	this->df[this->id]=1;
}

template<class T>
void AD<T>::set_f(T val){
	this->f=val;
}

template<class T>
void AD<T>::set_df_const(int n){
	this->df=vector<T>(n,0.0);
}

template<class T>
T AD<T>::getf(){
	return this->f;
}

template<class T>
T AD<T>::getDf(int index){
	return df[index];
}

template<class T>
vector<T> AD<T>::getGradient(){
	vector<T> gradient(this->df);
	return gradient;
}

template<class T>
Matrix<T> getJacobian(AD<T> funList[],int n,int m){
	Matrix<T> Mat=Matrix<T>(n,m);
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++){
			Mat.set_element(funList[i].getDf(j),i,j);
		}
	}
	return Mat;
}

template<class T>
AD<T> AD<T>::operator*(AD<T> g){
	AD<T> h;
	h.f=this->f*g.f;
	int n=this->df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=this->f*g.df[i]+g.f*this->df[i];
	}
	return h;
}

template<class T>
AD<T> AD<T>::operator+(AD<T> g){
	AD<T> h;
	h.f=this->f+g.f;
	int n=this->df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=this->df[i]+g.df[i];
	}
	return h;
}

template<class T>
AD<T> AD<T>::operator-(AD<T> g){
	AD<T> h;
	h.f=this->f-g.f;
	int n=this->df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=this->df[i]-g.df[i];
	}
	return h;
}

template<class T>
AD<T> AD<T>::operator/(AD<T> g){
	assert(g.f!=0);
	AD<T> h;
	h.f=this->f/g.f;
	int n=this->df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=(this->df[i]*g.f-this->f*g.df[i])/(g.f*g.f);
	}
	return h;
}

template<class T>
AD<T> AD<T>::operator^(AD<T> g){
	AD<T> h;
	assert(this->f!=0);
	h.f=pow(this->f,g.f);
	int n=this->df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=h.f*(g.f*(this->df[i])/this->f+g.df[i]*log(this->f));
	}
	return h;
}


template<class T>
AD<T> AD<T>::operator*(T s){
	AD<T> h;
	h.f=s*this->f;
	int n=this->df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=s*this->df[i];
	}
	return h;
}

template<class T>
AD<T> AD<T>::operator+(T s){
	AD<T> h;
	h.f=this->f+s;
	int n=this->df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=this->df[i];
	}
	return h;
}

template<class T>
AD<T> AD<T>::operator-(T s){
	AD h;
	h.f=this->f-s;
	int n=this->df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=this->df[i];
	}
	return h;
}

template<class T>
AD<T> AD<T>::operator/(T s){
	AD<T> h;
	assert(s!=0);
	h.f=this->f/s;
	int n=this->df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=this->df[i]/s;
	}
	return h;
}

template<class T>
AD<T> AD<T>::operator^(T s){
	AD<T> h;
	h.f=pow(this->f,s);
	int n=this->df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=s*pow(this->f,s-1)*this->df[i];
	}
	return h;
}

template<class T>
AD<T> operator*(T s,AD<T> g){
	AD<T> h;
	h.f=s*g.f;
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=s*g.df[i];
	}
	return h;
}

template<class T>
AD<T> operator+(T s,AD<T> g){
	AD<T> h;
	h.f=s+g.f;
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=g.df[i];
	}
	return h;
}

template<class T>
AD<T> operator-(T s,AD<T> g){
	AD<T> h;
	h.f=s-g.f;
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=-g.df[i];
	}
	return h;
}

template<class T>
AD<T> operator/(T s,AD<T> g){
	AD<T> h;
	assert(g.f!=0);
	h.f=s/g.f;
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=s*(-1/(g.f*g.f))*g.df[i];
	}
	return h;
}

template<class T>
AD<T> operator^(T s,AD<T> g){
	AD<T> h;
	h.f=pow(s,g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=pow(s,g.f)*log(s)*g.df[i];
	}
	return h;
}

template<class T>
AD<T> sin(AD<T> g){
	AD<T> h;
	h.f=sin(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=cos(g.f)*g.df[i];
	}
	return h;
}

template<class T>
AD<T> cos(AD<T> g){
	AD<T> h;
	h.f=cos(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=-sin(g.f)*g.df[i];
	}
	return h;
}

template<class T>
AD<T> tan(AD<T> g){
	AD<T> h;
	h.f=tan(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	assert(cos(g.f)!=0);
	for(int i=0;i<n;i++){
		h.df[i]=g.df[i]/(cos(g.f)*cos(g.f));
	}
	return h;
}

template<class T>
AD<T> cosec(AD<T> g){
	AD<T> h;
	assert(sin(g.f)!=0);
	h.f=1/sin(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	assert(tan(g.f)!=0&&sin(g.f)!=0);
	for(int i=0;i<n;i++){
		h.df[i]=-g.df[i]/(tan(g.f)*sin(g.f));
	}
	return h;
}

template<class T>
AD<T> sec(AD<T> g){
	AD<T> h;
	assert(cos(g.f)!=0);
	h.f=1/cos(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=(g.df[i]*tan(g.f))/cos(g.f);
	}
	return h;
}

template<class T>
AD<T> cot(AD<T> g){
	AD<T> h;
	assert(tan(g.f)!=0);
	h.f=1/tan(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	assert(sin(g.f)!=0);
	for(int i=0;i<n;i++){
		h.df[i]=-g.df[i]/(sin(g.f)*sin(g.f));
	}
	return h;
}

template<class T>
AD<T> arcsin(AD<T> g){
	assert(g.f>=-1 && g.f<=1);
	AD<T> h;
	h.f=asin(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=g.df[i]/sqrt(1-g.f*g.f);
	}
	return h;
}

template<class T>
AD<T> arccos(AD<T> g){
	assert(g.f>=-1 && g.f<=1);
	AD<T> h;
	h.f=acos(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=-g.df[i]/sqrt(1-g.f*g.f);
	}
	return h;
}

template<class T>
AD<T> arctan(AD<T> g){
	AD<T> h;
	h.f=atan(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=-g.df[i]/(1+g.f*g.f);
	}
	return h;
}

template<class T>
AD<T> sinh(AD<T> g){
	AD<T> h;
	h.f=sinh(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=g.df[i]*cosh(g.f);
	}
	return h;
}

template<class T>
AD<T> cosh(AD<T> g){
	AD<T> h;
	h.f=cosh(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=g.df[i]*sinh(g.f);
	}
	return h;
}

template<class T>
AD<T> tanh(AD<T> g){
	assert(cosh(g.f)!=0);
	AD<T> h;
	h.f=tanh(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=g.df[i]*(pow(cosh(g.f),2)-pow(sinh(g.f),2))/pow(cosh(g.f),2);
	}
	return h;
}

template<class T>
AD<T> log(AD<T> g){
	assert(g.f>0);
	AD<T> h;
	h.f=log(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=g.df[i]/g.f;
	}
	return h;
}

template<class T>
AD<T> abs(AD<T> g){
	AD<T> h;
	h.f=abs(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		if(g.f<0){
			h.df[i]=-g.df[i];
		}else{
			h.df[i]=g.df[i];
		}
	}
	return h;
}

template<class T>
AD<T> exp(AD<T> g){
	AD<T> h;
	h.f=exp(g.f);
	int n=g.df.size();
	h.df=vector<T>(n);
	for(int i=0;i<n;i++){
		h.df[i]=g.df[i]*exp(g.f);
	}
	return h;
}

