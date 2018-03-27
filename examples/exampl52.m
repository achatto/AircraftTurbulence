% Exampl52  Calculates covariance matrix as function of time for a
%           second order mass-spring-damper system perturbed by white 
%           noise, using the discrete-time recursive calculation method.
%
% Chapter 5 of lecture notes ae4-304
%
% Program revised August 1992, February 2004 [MM]

clc
clf
clear

disp('   Example 5.2                                                     ');
disp('                                                                   ');  
disp('   Calculation of the growth in time of the covariance             ');
disp('   matrix of a second order mass-spring-damper system              ');
disp('   perturbed by white noise. The covariance matrix C               ');
disp('   is calculated in time with:                                     ');
disp('                                                                   ');
disp('                            T                               T      ');
disp('   C  (k+1) = PHI(k) C  (k) PHI (k)   +   GAMMA(k) C  (k) GAMMA (k)');
disp('    xx                xx                            ww             ');
disp('                                                                   ');
disp('   This program can produce Figures 5.4 to 5.8 of the lecture      ');
disp('   notes: Aircraft Responses to Atmospheric Turbulence.            ');
disp('                                                                   ');
disp('                                                                   ');
disp('   2nd Order Model Definition:                                     ');
disp('                                                                   ');

% CT SYSTEM DYNAMICS
w0      = input('   Give undamped natural frequency [rad/s]     : ');
zeta    = input('   Give damping ratio                          : ');

num=1;
den=[1/(w0)^2 2*zeta/w0 1];

m=1;k=w0^2;
c=zeta*2*m*w0;

A=[0 1;-k/m -c/m];  % state-space representation of second
B=[0 1/m]';         % order system
C=[1 0];
D=[0];

% DEFINE TIME AXIS
dt= .1;	                    % sample time 0.1 seconds
t = [dt:dt:15];	            % time axis
N = 15/dt;                  % number of samples

% DISCRETIZE SYSTEM MATRICES
[Phi,Gamma]=c2d(A,B,dt);    % discretizing using MATLAB c2d command

% DEFINE WHITE NOISE CHARACTERISTICS
Wc      = input('   Enter CT white noise intensity              : ');
answ    = input('   Stepwise change of noise intensity ? (y/n)  :','s');

if answ=='y',
  Q     = input('   Enter Q (0<Q<1)                             : ');
  answ1 = input('   Noise intensity W=0 after t=T(1-Q) (y/n)    :','s');
  answ2 = input('   Noise intensity W=2*Wc after t=T(1-Q) (y/n) :','s');
  M=Q*N;
  for k=1:N-M;
    W(k)=Wc/dt;             % always apply equation (5.45)
  end
  for k=N-M+1:N;
    if answ1=='y',
     W(k)=0;
    end
    if answ2=='y',
     W(k)=2*Wc/dt;
    end
  end
end
if answ=='n';
  for k=1:N;
    W(k)=Wc/dt;
  end
end

% DEFINE INITIAL CONDITIONS
Cx1x1(1)= input('   Give initial value of Cyy(1,1)              : ');
Cx1x2(1)= input('   Give initial value of Cyy(1,2)=C(2,1)       : ');
Cx2x2(1)= input('   Give initial value of Cyy(2,2)              : ');
Cxx=[Cx1x1(1) Cx1x2(1);Cx1x2(1) Cx2x2(1)];

% DISCRETE SOLUTION Cxx(k+1)=Phi*Cxx(k)*Phi' + Gamma*W*Gamma';
for k=1:N-1;
  Cxx        = Phi*Cxx*Phi'+Gamma*W(k)*Gamma';
  Cx1x1(k+1) = Cxx(1,1); Cx1x2(k+1)=Cxx(2,1); Cx2x2(k+1)=Cxx(2,2);
end

ref=zeros(1,N);

% PLOT RESULTS
clf
subplot(221);
axis([0 15 -0.5 1.5]);
plot(t,Cx1x1,'-',t,ref,'-');
xlabel('time [s]');ylabel('Cx1x1 (t)');

subplot(222);
plot(t,Cx1x2,'-',t,ref,'-');
xlabel('time [s]');ylabel('Cx1x2 (t)');

subplot(224);
plot(t,Cx2x2,'-',t,ref,'-');
xlabel('time [s]');ylabel('Cx2x2 (t)');

% EOF