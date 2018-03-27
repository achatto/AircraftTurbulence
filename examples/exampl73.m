% Filename : exampl73.m
%
% Calculates the power spectral density of the normal acceleration.

clf, clc, clear;

disp('   Example 7.3');
disp('   ');
disp('   Calculation of the power spectral density of the normal');
disp('   acceleration due to longitudinal and vertical turbulence.');
disp('   Also, the effect of a lagfree autopilot will be investigated.');
disp('   ');
disp('   This program produces Figures 7-12 and 7-13 of the lecture');
disp('   notes: Aircraft Responses to Atmospheric Turbulence.');

% GET AIRCRAFT DYNAMICS
cit2s

% DEFINE FREQUENCY AXES
N = 1000;   
omega = logspace(-3,3,N);     % frequency axis

% MISCELLANEOUS
D = [0 0 0];
g = 9.80665;

% HORIZONTAL TURBULENCE
C = [0 1 0 0 0 0 0];
[numa,den] = ss2tf(A,B,C,D,2);            % transfer function alpha-w1

C = [ 0 0 1 0 0 0 0];
[numt,den] = ss2tf(A,B,C,D,2);            % transfer function theta-w1

numn = numt-numa;                         % transfer function an-w1
numn(length(numn)+1) = 0;
numn = (V/g)*numn;

% COMPUTE FREQUENCY RESPONSE FUNCTION
[mag,phase] = bode(numn,den,omega);
Snn = mag.*mag;

% VERTICAL TURBULENCE
C = [0 1 0 0 0 0 0];
[numa,den] = ss2tf(A,B,C,D,3);            % transfer function alpha-w3

C = [0 0 1 0 0 0 0];
[numt,den] = ss2tf(A,B,C,D,3);            % transfer function theta-w3

numn = numt-numa;                         % transfer function an-w3
numn(length(numn)+1) = 0;
numn = (V/g)*numn;

% COMPUTE FREQUENCY RESPONSE FUNCTION
[mag,phase] = bode(numn,den,omega);
Snn1 = mag.*mag;

% VERTICAL TURBULENCE AIRCRAFT WITH PITCH ATTITUDE HOLD SYSTEM
C = [0 1 0 0 0 0 0];
[numa,den] = ss2tf(At,B,C,D,3);           % transfer function alpha-w3

C = [0 0 1 0 0 0 0];
[numt,den] = ss2tf(At,B,C,D,3);           % transfer function theta-w3

numn = numt-numa;                         % transfer function an-w3
numn(length(numn)+1) = 0;
numn = (V/g)*numn;

% COMPUTE FREQUENCY RESPONSE FUNCTION
[mag,phase] = bode(numn,den,omega);
Snnt1 = mag.*mag;

% PLOT POWER SPECTRAL DENSITIES
clf;
axis('square')
loglog(omega,Snnt1,'--',omega,Snn1);
xlabel('omega [rad/s]'); ylabel('Snn');
text(0.5,0.85,'--- pitch attitude hold','sc')
text(0.5,0.80,'-   elevator fixed','sc')
title('Power Spectral Density of Normal Acceleration');
pause

clf;
axis('square')
loglog(omega,Snn,'--',omega,Snn1);
xlabel('omega [rad/s]'); ylabel('Snn');
text(0.48,0.85,'--- horizontal turbulence','sc')
text(0.48,0.80,'-   vertical turbulence','sc')
title('Power Spectral Density of Normal Acceleration');
pause

% CALCULATION OF VARIANCES USING VERY CRUDE INTEGRATION
dw = diff(omega);
dw(length(dw)+1)=0;   % make vector length equal to N again

disp(' ');
disp(' Variance of n due to horizontal turbulence')
varn    = sum(Snn'.*dw)/pi
pause
disp(' Variance of az due to horizontal turbulence')
varaz   = sum(Snn'.*dw)*g^2/pi
pause
disp(' Variance of n due to vertical turbulence')
varn1   = sum(Snn1'.*dw)/pi
pause
disp(' Variance of az due to vertical turbulence')
varaz1  = sum(Snn1'.*dw)*g^2/pi
pause
disp(' Variance of n due to vertical turbulence for aircraft')
disp(' with pitch attitude hold system')
varnt1  = sum(Snnt1'.*dw)/pi
pause
disp(' Variance of az due to vertical turbulence for aircraft')
disp(' with pitch attitude hold system')
varazt1 = sum(Snnt1'.*dw)*g^2/pi
pause
