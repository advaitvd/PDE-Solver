#pragma once
#include<vector>
#include<iostream>
#include "Matrix.h"
using namespace std;

template<class T>
class LinearSolver
{
	private:
	public:
		virtual vector<T> solve(Matrix<T>,vector<T>)=0;
		virtual ~LinearSolver(){}
};

template<class T>
class GaussElimination: public LinearSolver<T>{
	private:
	public:
		virtual vector<T> solve(Matrix<T> M, vector<T> u);
		virtual ~GaussElimination(){}
};

template<class T>
class LU_Decomposition: public LinearSolver<T>{
	private:
	public:
		virtual vector<T> solve(Matrix<T> M, vector<T> u);
		virtual ~LU_Decomposition(){}
};

template<class T>
class TriDiagonal: public LinearSolver<T>{
	private:
	public:
		virtual vector<T> solve(Matrix<T> M, vector<T> u);
		virtual ~TriDiagonal(){}
};

template<class T>
class Gauss_Jacobi: public LinearSolver<T>{
	private:
	public:
		virtual vector<T> solve(Matrix<T> M, vector<T> u);
		virtual ~Gauss_Jacobi(){}
};

template<class T>
class Gauss_Seidal: public LinearSolver<T>{
	private:
	public:
		virtual vector<T> solve(Matrix<T> M, vector<T> u);
		virtual ~Gauss_Seidal(){}
};

template<class T>
class SOR: public LinearSolver<T>{
	private:
	public:
		virtual vector<T> solve(Matrix<T> M, vector<T> u);
		virtual ~SOR(){}
};

template<class T>
vector<T> SOR<T>::solve(Matrix<T> M, vector<T> u){
    vector<T> x,xn;
    int i,j,flag;
    T sum,eps=0.001,w=1;
    
    x=vector<T>(u.size());
    xn=vector<T>(u.size());
    
    
    for(i=0;i<u.size();i++)
        x[i]=0;
    do
    {
        for(i=0;i<u.size();i++)
        {
            sum=u[i]*w+M.get_element(i,i)*x[i];
            for(j=0;j<u.size();j++)
            {
                if(j<i)
                    sum-=M.get_element(i,j)*xn[j]*w;
                else if(j>=i)
                    sum-=M.get_element(i,j)*x[j]*w;
                xn[i]=sum/M.get_element(i,i);
            }
        }
        flag=0;
        for(i=0;i<u.size();i++)
            if(fabs(x[i]-xn[i])>eps)
                flag=1;
        if(flag==1)
            for(i=0;i<u.size();i++)
                x[i]=xn[i];
    }while(flag==1);
    return xn;
}

template<class T>
vector<T> Gauss_Seidal<T>::solve(Matrix<T> M, vector<T> u)
{
    vector<T> x,xn;
    int test,i,j;
    T g,sum,eps=0.0001;
    //cout<<"SATYA";
//    for(i=0;i<u.size();i++){
//       for(j=0;j<u.size();j++){
//         cout<<M.get_element(i,j)<<" ";
//       }
//       cout<<"\n";
//    }
      

    x=vector<T>(u.size());
    xn=vector<T>(u.size());
    for(i=0;i<u.size();i++)
    {
        x[i] = 0;
    }
    
    for(i=0;i<u.size();i++)
        x[i]=0;

        
    test=0;
    for(i=0;i<u.size();i++)
    {
        sum=0;
        for(j=0;j<u.size();j++)
            if(i!=j)
                sum+=fabs(M.get_element(i,j));
        if(sum>fabs(M.get_element(i,i)))
            test=1;
    }
 // Diagonally dominance verification by columns.
    if(test==1)
    {
        test=0;
        for(j=0;j<u.size();j++)
        {
            sum=0;
            for(i=0;i<u.size();i++)
                if(i!=j)
                    sum+=fabs(M.get_element(i,j));
            if(sum>fabs(M.get_element(j,j)))
                test=1;
        }
    }

    if(test==1)
    {
        cout<<"The co-efficient matrix is not diagonally dominant\n";
        cout<<"The Gauss-Seidel method doesn't converge surely\n";
        exit(0);
    }
    do
    {
        for(i=0;i<u.size();i++)
        {
            sum=u[i];
            for(j=0;j<u.size();j++)
            {
                if(j<i)
                    sum-=M.get_element(i,j)*xn[j];
                else if(j>i)
                    sum-=M.get_element(i,j)*x[j];
            }
            xn[i]=sum/M.get_element(i,i);
        }
        test=0;
        for(i=0;i<u.size();i++)
            if(fabs(x[i]-xn[i])>eps)
                test=1;
        if(test==1)
            for(i=0;i<u.size();i++)
                x[i] = xn[i];
    }while(test==1);
    return xn;
}

template<class T>
vector<T> Gauss_Jacobi<T>::solve(Matrix<T> M, vector<T> u)
{
    vector<T> x,xn;
    int i,j,flag;
    T sum,eps;
    
    x=vector<T>(u.size());
    xn=vector<T>(u.size());
    cout<<"Enter the accuracy u want\n";
    cin>>eps;
    for(i=0;i<u.size();i++)
        x[i]=0;
    flag=0;
    for(i=0;i<u.size();i++)
    {
        sum=0;
        for(j=0;j<u.size();j++)
            if(i!=j)
                sum+=fabs(M.get_element(i,j));
        if(sum>fabs(M.get_element(i,i)))
            flag=1;
    }
    if(flag==1)
    {
        flag=0;
        for(j=0;j<u.size();j++)
        {
            sum=0;
            for(i=0;i<u.size();i++)
                if(i!=j)
                    sum+=fabs(M.get_element(i,j));
            if(sum>fabs(M.get_element(j,j)))
                flag=1;
        }
    }
    if(flag==1)
    {
        cout<<"The co-efficient matrix is not diagonally dominant\n";
        cout<<"The Gauss-Jacobi method doesn't converge surely\n";
        exit(0);
    }
    do
    {
        for(i=0;i<u.size();i++)
        {
            sum=u[i];
            for(j=0;j<u.size();j++)
                if(j!=i)
                    sum-=M.get_element(i,j)*x[j];
            xn[i]=sum/M.get_element(i,i);
        }
        flag=0;
        for(i=0;i<u.size();i++)
            if(fabs(x[i]-xn[i])>eps)
                flag=1;
        if(flag==1)
            for(i=0;i<u.size();i++)
                x[i]=xn[i];
    }while(flag==1);
    return xn;
}

template<class T>
vector<T> TriDiagonal<T>::solve(Matrix<T> M, vector<T> u)
{
    vector<T> x(u.size());
    T A[u.size()-1],B[u.size()],C[u.size()-1],D[u.size()];
    T C_star[u.size()-1],D_star[u.size()];
    int i,j;
    
    for(i=0;i<u.size()-2;i++)
        for(j=i;j<u.size()-2;j++)
    {
        if(M.get_element(i,j+2)!=0)
        {
            cout<<"Method can't be applied";
            exit(0);
        }

    }
    for(i=u.size()-1;i>1;i--)
        for(j=i;j>1;j--)
    {
        if(M.get_element(i,j-2)!=0)
        {
            cout<<"Method can't be applied";
            exit(0);
        }
    }
    for(i=0;i<u.size();i++)
    {
        B[i]=M.get_element(i,i);
        D[i]=u[i];
    }
    for(i=1;i<u.size();i++)
    {
       A[i]=M.get_element(i,i-1);
    }
    for(i=0;i<u.size()-1;i++)
    {
       C[i]=M.get_element(i,i+1);
    }

    C_star[0]=C[0]/B[0];
    D_star[0]=D[0]/B[0];
    for(i=1;i<u.size();i++)
    {
        C_star[i]=C[i]/(B[i]-A[i]*C_star[i-1]);
        D_star[i]=(D[i]-A[i]*D_star[i-1])/(B[i]-A[i]*C_star[i-1]);
    }
    x[u.size()-1]=D_star[u.size()-1];
    for(i=u.size()-2;i>=0;i--)
    {
        x[i]=D_star[i]-C_star[i]*x[i+1];
    }
      return x;
}


template<class T>
vector<T> LU_Decomposition<T>::solve(Matrix<T> M, vector<T> u)
{
    vector<T> x(u.size()), z(u.size());
    Matrix<T> l(u.size(),u.size()), u1(u.size(),u.size());
    int i,j,k;
    
    for(i=0;i<u.size();i++)
    	l.set_element(M.get_element(i,0),i,0);
    for(j=1;j<u.size();j++)
    	u1.set_element(M.get_element(0,j)/l.get_element(0,0),0,j);
    for(i=0;i<u.size();i++)
    	u1.set_element(1,i,i);
    for(i=1;i<u.size();i++)
        for(j=1;j<u.size();j++)
            if(i>=j)
            {
            	l.set_element(M.get_element(i,j),i,j);
                for(k=0;k<=j-1;k++){
                	l.set_element(l.get_element(i,j)-l.get_element(i,j)*u1.get_element(k,j),i,j);
                }
            }
            else
            {
                u1.set_element(M.get_element(i,j),i,j);
                for(k=0;k<=i-1;k++)
                	u1.set_element(u1.get_element(i,j)-l.get_element(i,k)*u1.get_element(k,j),i,j);
                u1.set_element(u1.get_element(i,j)/l.get_element(i,i),i,j);
            }
//    cout<<"The lower triangular matrix L:"<<endl;
//    for(i=0;i<u.size();i++)
//    {
//        for(j=0;j<=i;j++)
//            cout<<"\t"<<l.get_element(i,j);
//        cout<<endl;
//    }
//    cout<<"\nThe upper triangular matrix U:"<<endl;
//    for(i=0;i<u.size();i++)
//    {
//        for(j=0;j<i;j++)
//            cout<<"\t";
//        for(j=i;j<u.size();j++)
//            cout<<"\t"<<u1.get_element(i,j);
//        cout<<endl;
//    }
    z[0]=u[0]/l.get_element(0,0);
    for(i=1;i<u.size();i++)
    {
        z[i]=u[i];
        for(j=0;j<=i-1;j++)
            z[i]-=l.get_element(i,j)*z[j];
        z[i]/=l.get_element(i,i);
    }
    x[u.size()-1]=z[u.size()-1];
    for(i=u.size()-2;i>=0;i--)
    {
        x[i]=z[i];
        for(j=i+1;j<u.size();j++)
            x[i]-=u1.get_element(i,j)*x[j];
    }
    return x;
}

template<class T>
vector<T> GaussElimination<T>::solve(Matrix<T> M, vector<T> u)
{
    vector<T> x(u.size());
    T temp;
    int i,j,k;
    for(i=0;i<u.size();i++)
    {
        temp = fabs(M.get_element(i,i));
        k = i;
        for (j=i+1;j<u.size();j++)
            if(temp<fabs(M.get_element(j,i)))
            {
                temp = fabs(M.get_element(j,i));
                k = j;
            }
        if (fabs(M.get_element(k,i))<0.00001)
        {
            cout << "The matrix is singular: The system has either no solution or infinitely many solution";
            exit(0);
        }
        if(k!=i)
        {
            for(j=0;j<u.size();j++)
            {
                temp = M.get_element(k,j);
                M.set_element(M.get_element(i,j),k,j);
                M.set_element(temp,i,j);
            }
            temp = u[k];
            u[k] = u[i];
            u[i] = temp;
        }
        for(j=i+1;j<u.size();j++)
        {
            temp = M.get_element(j,i)/M.get_element(i,i);
            for(k=0;k<u.size();k++)
            	M.set_element(M.get_element(j,k) - temp*M.get_element(i,k),j,k);
            u[j] = u[j] - temp*u[i];
        }
    }
    x[u.size()-1] = u[u.size()-1] / M.get_element(u.size()-1, u.size()-1);
    for(i=u.size()-2;i>=0;i--)
    {
        x[i] = u[i];
        for(j=i+1;j<u.size();j++)
            x[i] = x[i] - M.get_element(i,j)*x[j];
        x[i] = x[i] / M.get_element(i,i);
    }
    return x;
}


