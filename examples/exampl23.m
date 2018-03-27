% Exampl23     Calculate and plot the 2-dimensional Normal (Gaussian)
%              probability density function.
%
% Chapter 2 of the lecture notes ae4-304.
%
% Program revised August 1992, February 2004 [MM]

clc
clf
clear

disp('   Example 2.3');
disp(' ')
disp('   Calculate and plot the 2-dimensional Normal (Gaussian)');
disp('   probability density function.');
disp('   ');
disp('   This program can produce Figures 2-15, 2-16 and 2-17 of');
disp('   the lecture notes ae4-304.');

x = -3:0.1:3; y = x;

% Definition of distribution parameters
mx  = input('   Average value of stochastic variable x           : ');
sx  = input('   Standard deviation of x                          : ');
my  = input('   Average value of stochastic variable y           : ');
sy  = input('   Standard deviation of y                          : ');
Kxy = input('   Correlation coefficient between x and y, 0<Kxy<1 : ');

if abs(Kxy) > 0.99 | abs(Kxy) < 0
    break
end    

fmax=1/(2*pi*sx*sy*sqrt(1-Kxy^2));
for i=1:length(x)
    for j=1:length(y)
        G=((x(i)-mx)^2/sx^2-2*Kxy*(x(i)-mx)*(y(j)-my)/...
          (sx*sy)+(y(j)-my)^2/sy^2);
        G=G/(1-Kxy^2);
        fxy(j,i)=fmax*exp(-G/2);
    end
end

clf
mesh(x,y,fxy);
title('The 2-dimensional Normal p.d.f.')
pause

contour(x,y,fxy);
title('The 2-dimensional Normal p.d.f. viewed from above')

% EOF