% Filename : exampl62.m
%
% Simulation of atmospheric turbulence using Dryden model

clf, clc, clear

disp('   Example 6.2');
disp('   Simulation of vertical gust velocity using Dryden model.');
disp('   ');
disp('   This example produces Figure 6-19 of the lecture notes:');
disp('   Aircraft Responses to Atmospheric Turbulence.');

sigma = input('   Enter turbulence intensity sigma [m/s]  (0.282) : ');
Lg1   = input('   Enter turbulence scale length Lg1 [m]   (  150) : ');
Lg2   = input('   Enter turbulence scale length Lg2 [m]   ( 1500) : ');
V     = input('   Enter airspeed V [m/s]                  (   35) : ');

% Define time basis
dt = 0.1; T = 120;
t  = [0:dt:T];
N = length(t);

% White noise input
w = randn(1,N)/sqrt(dt);  % note: divide by sqrt(dt), with dt the sample time
                          %       because of lsim characteristics

% Forming filter characteristics equation (6.41)
rat = V/Lg1;
A = [0 1;-rat^2 -2*rat];
B = sigma*[sqrt(3*rat);(1-2*sqrt(3))*sqrt((rat^3))];
C = [1 0];
D = [0];

% Output turbulence velocity
wg = lsim(A,B,C,D,w,t);

% Forming filter characteristics equation (6.41)
rat = V/Lg2;
A = [0 1;-rat^2 -2*rat];
B = sigma*[sqrt(3*rat);(1-2*sqrt(3))*sqrt((rat^3))];
C = [1 0];
D = [0];

% Output turbulence velocity
wgg = lsim(A,B,C,D,w,t);

% Plot the results
subplot(2,1,1)
plot(t,w);
xlabel('time [s]');
ylabel('w');
title('White Noise Filter Input');

subplot(2,1,2)
plot(t,wg,t,wgg,'--');
xlabel('time [s]');
ylabel('wg [m/s]');

title('Vertical Gust Velocity, -  : Lg = 150 m, -- : Lg = 1500 m');
