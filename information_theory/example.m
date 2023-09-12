clear all;
clc; 
close all;

ref_x=[1 2 3 4;1.2 2.2 3.2 4.2;1.3 2.3 3.3 4.3];

ref_y=[1 2 3 4;1.2 2.2 3.2 4.2;1.3 2.3 3.3 4.3]+2;
x=ref_x(:);
y=ref_y(:);

plot(x,y,'o')


xy=[x y]

