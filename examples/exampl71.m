% Filename : examp71.m
%
% Simulation of aircraft symmetric response to atmospheric turbulence.

clf, clc, clear

disp('   Example 7.1');
disp('   ');
disp('   Simulation of symmetric gust response of the');
disp('   Cessna Ce-500 "Citation" at an airspeed of 59.9 [m/s], ');
disp('   CRUISE condition.');
disp('   ');
disp('   This program produces Figure 7-6 of the lecture notes:');
disp('   Aircraft Responses to Amospheric Turbulence.');

cit2s       % loading A,B matrices

% TIME AXIS INPUT VECTOR DEFINITION
dt = input('   Give sampling time interval dt             (0.01) : ');
T  = 60; t = [0:dt:T]; N = length(t);

% INPUT VECTOR DEFINITION
nn = zeros(1,N);             % zero input elevator
w1 = randn(1,N)/sqrt(dt);    % scaled input hor. turbulence,
                             % note the sqrt(dt) because of lsim
w3 = randn(1,N)/sqrt(dt);    % scaled input vert. turbulence,
                             % note the sqrt(dt) because of lsim
u  = [nn' nn' w3'];          % input vector definition (vertical
                             % turbulence only, can be changed).

% SIMULATION OF MOTION VARIABLES
C = eye(7); D = zeros(7,3);
y = lsim(A,B,C,D,u,t);

% PLOTTING RESULTS
subplot(2,1,1);
plot(t,y(:,1)*180/pi)
xlabel('time [s]'); ylabel('u/V [-]'); title('airspeed deviation');

subplot(2,1,2);
plot(t,y(:,2)*180/pi)
xlabel('time [s]'); ylabel('alpha [deg]'); title('angle of attack');
pause

subplot(2,1,1);
plot(t,y(:,3)*180/pi)
xlabel('time [s]'); ylabel('theta [deg]'); title('pitch angle');

subplot(2,1,2);
plot(t,y(:,4)*180/pi)
xlabel('time [s]'); ylabel('qc/V [deg]'); title('pitch rate');
