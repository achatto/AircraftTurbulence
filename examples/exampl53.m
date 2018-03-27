% Exampl53    Calculates covariance matrix as function of time for a
%             second order mass-spring-damper system perturbed by white
%             noise using the impulse response method.
%
% Chapter 5 of lecture notes ae4-304
%
% Program revised August 1992, February 2004 [MM]

clc
clf
clear

disp('   Example 5.3                                                ');
disp('                                                              ');
disp('   Calculation of the growth in time of the covariance        ');
disp('   matrix of a second order mass-spring-damper system         ');
disp('   perturbed by white noise using the impulse response method:');
disp('                                                              ');
disp('               _ t                                            ');
disp('              |          T                                    ');
disp('     C  (k) = | h(v).h(v)  dv                                 ');
disp('      xx     _|                                               ');
disp('            0                                                 ');
disp('                                                              ');
disp('   This program can produce Figure 5.9 of the lecture notes:  ');
disp('   Aircraft Responses to Atmospheric Turbulence.              ');
disp('                                                              ');
disp('                                                              ');
disp('   2nd Order Model Definition:                                ');
disp('                                                              ');

% CT SYSTEM DYNAMICS
w0      = input('   Give undamped natural frequency [rad/s] : ');
zeta    = input('   Give damping ratio                      : ');

num=1;
den=[1/(w0)^2 2*zeta/w0 1];

m=1;k=w0^2;
c=zeta*2*m*w0;

A=[0 1;-k/m -c/m];  % state-space representation of a CT second
B=[0 1/m]';         % order system
C=[1 0];
D=[0];

% TIME AXIS
dt=.1; T=15; t=[0:dt:T]; N=length(t);

% COMPUTE CT SYSTEM IMPULSE RESPONSE
u=0*ones(1,N);              % zero input
x0=B;                       % initial condition
h=lsim(A,B,C,D,u,t,x0);     % calculation impulse response as a resonse
                            % to an initial condition.

% SQUARE IMPULSE RESPONSE
hsq=h.*h;                   % squared impulse response

% INTEGRATE THE SQUARED IMPULSE RESPONSE
% done rather crudely here
vary=0;
for i=1:N-1
  vary(i+1)=vary(i)+dt*hsq(i);	% integrated squared impulse response
end

% PLOT RESULTS
clf
subplot(2,2,1)
axis([0 15 -.5 1.5]);
plot(t,h,'-');
xlabel('time [s]');ylabel('h (t)');title('Impulse Response');

subplot(2,2,2)
plot(t,hsq,'-');
xlabel('time [s]');ylabel('h^2 (t)');title('Squared Impulse Response');
axis([0 15 -.5 1.5]);

subplot(2,2,3)
plot(t,vary,'-');xlabel('time [s]');ylabel('Cx1x1 (t)');
title('Variance');

% EOF