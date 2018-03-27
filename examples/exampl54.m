% Exampl54  Calculates covariance matrix of a second order mass-spring-damper
%           system perturbed by white noise, using the Monte Carlo
%           method (ensemble average = average over all realizations)
%
% Chapter 5 of lecture notes ae4-304
%
% Program revised August 1992, February 2004 [MM]

clc
clf
clear

disp('   Example 5.4                                                ');
disp('                                                              ');
disp('   Calculation of the Covariance matrix of a second order     ');
disp('   mass-spring-damper system perturbed by white noise using   ');
disp('   Monte Carlo Method.                                        ');
disp('   The covariance matrix C is calculated with:                ');
disp('                                                              ');
disp('             1   N                    T                       ');
disp('      C   = --- SUM (x(i)-m )(x(i)-m )                        ');
disp('       xx   N-1 i=1        x        x                         ');
disp('                                                              ');
disp('   This program can produce Figures 5.10-5.12 in the lecture  ');
disp('   notes: Aircraft Responses to Atmospheric Turbulence.       ');
disp('                                                              ');
disp('   2nd order dynamic model Definition:                        ');
disp('                                                              ');

% DEFINE CT SYSTEM DYNAMICS
w0   =  input('   Give undamped natural frequency [rad/s]     : ');
zeta =  input('   Give damping ratio                          : ');

m=1;k=w0^2; c=zeta*2*m*w0;

A=[0 1;-k/m -c/m];
B=[0 1/m]';

% SET TIME AXIS
dt   =  input('   Give sampling time interval dt              : ');
T    =  input('   Give total time interval T                  : ');
t=[0:dt:T-dt];
N=T/dt;

% COMPUTE DT EQUIVALENT
[Phi,Gamma]=c2d(A,B,dt);

% DEFINE CT WHITE NOISE CHARACTERISTICS
Wc   =  input('   Give CT white noise intensity               : ');

Wd = Wc/dt;     % NOTE: divide by sample time dt, this was not done 
                % correctly in the lecture notes ae4-404, May 1998
                % and earlier!! (Figures 5.8, 5.9 
                % and 5.10 are now correct)

% DEFINE # of EXPERIMENTS (MONTE CARLO METHOD)
NN   =  input('   Give number of experiments                  : ');

w = sqrt(Wd).*randn(N,NN);   % generate realizations of the random noise
                     
x1=zeros(N,NN);
x2=zeros(N,NN);
x=[x1(1,:);x2(1,:)];
w1=zeros(2,1);
xx=zeros(2,1);
for j=1:NN;
  for l=1:2;
    xx(l,1)=x(l,j);
  end
  for i=1:N-1;
    for k=1:2;
      w1(k,1)=Gamma(k,1)*w(i,j);
    end
    xx=Phi*xx+w1;
    x1(i+1,j)=xx(1,1);
    x2(i+1,j)=xx(2,1);
  end
end  
mean1=zeros(N,1);
mean2=zeros(N,1);
for i=1:N;
  for j=1:NN;
    mean1(i,1)=mean1(i,1)+x1(i,j);
    mean2(i,1)=mean2(i,1)+x2(i,j);
  end
end
mean1=(1/NN).*mean1;
mean2=(1/NN).*mean2;
Cx1x1=zeros(N,1);
Cx1x2=zeros(N,1);
Cx2x2=zeros(N,1);
for i=1:N;
  for j=1:NN;
    Cx1x1(i,1)=Cx1x1(i,1)+(x1(i,j)-mean1(i,1))*(x1(i,j)-mean1(i,1));
    Cx1x2(i,1)=Cx1x2(i,1)+(x1(i,j)-mean1(i,1))*(x2(i,j)-mean2(i,1));
    Cx2x2(i,1)=Cx2x2(i,1)+(x2(i,j)-mean2(i,1))*(x2(i,j)-mean2(i,1));
  end
  Cx1x1(i,1)=Cx1x1(i,1)/(NN-1);
  Cx1x2(i,1)=Cx1x2(i,1)/(NN-1);
  Cx2x2(i,1)=Cx2x2(i,1)/(NN-1);
end

% PLOT RESULTS
clf
subplot(2,1,1);
plot(t,mean1);xlabel('time (s)');ylabel('mean1(t)');

subplot(2,1,2);
plot(t,mean2);xlabel('time (s)');ylabel('mean2(t)');
pause

% ANALYTICAL CALCULATIONS
Cy1y1(1)=Cx1x1(1,1);Cy1y2(1)=Cx1x2(1,1);Cy2y2(1)=Cx2x2(1,1);
Cyy=[0 0;
     0 0];
for i=1:N-1;
    Cyy=Phi*Cyy*Phi'+Gamma*Wd*Gamma';
    Cy1y1(i+1)=Cyy(1,1);
    Cy1y2(i+1)=Cyy(1,2);
    Cy2y2(i+1)=Cyy(2,2);
end

% SHOW ANALYTICAL CALCULATIONS AND "MONTE-CARLO" DATA
clf
subplot(2,2,1);
plot(t,Cx1x1,'--',t,Cy1y1); xlabel('time (s)'); ylabel('Cx1x1(t)');

subplot(2,2,2);
plot(t,Cx1x2,'--',t,Cy1y2); xlabel('time (s)'); ylabel('Cx1x2(t)');

subplot(2,2,4);
plot(t,Cx2x2,'--',t,Cy2y2); xlabel('time (s)'); ylabel('Cx2x2(t)');

% EOF