% Filename : exampl61.m
%
% Calculation of correlation coefficients of turbulence
% velocities according to Batchelor using Dryden correlation
% functions at two points A and B.

clc, clf, clear

disp('   Example 6.1');
disp('   Calculates the correlation coefficient between velocity');
disp('   vectors at two points as a function of the longitudinal');
disp('   scale length L according to Batchelor using Dryden');
disp('   correlation functions f(.) and g(.).');
disp('   ');
disp('   This program produces Figure 6-17 of the lecture notes:');
disp('   Aircraft Responses to Atmospheric Turbulence.');

x1 = input('   Give x-separation [m] between the two points (-40) : ');
x2 = input('   Give y-separation [m] between the two points ( 20) : ');
x3 = input('   Give z-separation [m] between the two points (-10) : ');
xi = [x1,x2,x3];

disp('   ');
disp('   Velocity directions:');
disp('   1 --- Longitudinal');
disp('   2 --- Lateral');
disp('   3 --- Normal');
disp('   ');
uA = input('   Give velocity direction 1st point              (3) : ');
uB = input('   Give velocity direction 2nd point              (3) : ');

% Turbulence velocity directions for calculation of correlation;
% 1: longitudinal;
% 2: lateral;
% 3: normal;

if uA == uB
   delta = 1; 
else
   delta = 0; 
end;

Lg = logspace(1,4,50);		% running Lg from 10-10000m

for i = 1:50
    L = Lg(i);
    % Specific functions f(xi) and g(xi) according to Dryden;
    f = exp(-norm(xi)/L);
    g = f*(1-norm(xi)/(2*L));

    % Correlation according to Batchelor;
    K(i) = ((f-g)/norm(xi)^2)*xi(uA)*xi(uB)+delta*g;
end

% PLOTTING RESULTS
semilogx(Lg,K);
axis('square'); axis([10^1 10^4 -.2 1]);
xlabel('scale length of turbulence Lg [m]');
ylabel('correlation coefficient K');
grid