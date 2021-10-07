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