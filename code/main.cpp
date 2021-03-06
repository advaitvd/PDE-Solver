#include<iostream>
#include<fstream>

#include "AD.h"
#include "Discretizer.h"
using namespace std;
/*
Please Enter a test case as described below:
a b c d f1 f2 g1 g2
A:string
B:string
C:string
D:string
E:string
The string describing expressions for A,B,C,D,E must be in terms of x,y and u and mathematical operators.

Common functions like sin,cos,tan,cosec,cot,sec,abs,log,etc are also allowed.

*/

/*
Testcases given in the problem pdf. Use these as input in the terminal.
testcase1:
0 1 0 1 0 0 0 0
1
1
0
0
-2

testcase2:
0 1 0 1 0 0 0 0
1
1
0
0
-x*y

testcase3:
0 1 0 2 0 0 0 0
1
1
u
1
-(u^2-x^4-2*x^2-x+2)

*/

int main(){
	
	double a,b,c,d,f1,f2,g1,g2;
	cin>>a>>b>>c>>d>>f1>>f2>>g1>>g2;
	string A,B,C,D,E;
	cin>>A>>B>>C>>D>>E;
	if(A[0]=='-'){
		A="0"+A;
	}
	if(B[0]=='-'){
		B="0"+B;
	}
	if(C[0]=='-'){
		C="0"+C;
	}
	if(D[0]=='-'){
		D="0"+D;
	}
	if(E[0]=='-'){
		E="0"+E;
	}
	
	int Nx=16;
	int Ny=16;
	double initial_guess=0.1;
	string non_linear_solver="Newton";
	string linear_solver="Gauss_Seidal";
	
	Discretizer<double> M(initial_guess,a,b,c,d,f1,f2,g1,g2,A,B,C,D,E,Nx,Ny,1,non_linear_solver,linear_solver);
	
	Matrix<double> temp=M.solve();
	temp.round_zeros();
	
	cout<<temp;
	
	ofstream fout("output.txt");
	fout<<temp;
	fout.close();


//The following is the matlab  code for plotting  of the result
/*
%%
clc;clear;

%%
z=importdata('output.txt');
[n,m]=size(z);
x=linspace(0,1,n);
y=linspace(0,1,m);

surf(x,y,z);
xlabel('x');
ylabel('y');
zlabel('u');
*/	
	return 0;
}




