% Filename : exampl82.m
% 

clc, clf, clear

disp('   Example 8.2');
disp('   ');
disp('   Calculation of the power spectral density of the roll angle');
disp('   due to asymmetric vertical turbulence.');
disp('   ');
disp('   The effect of an autopilot which keeps the roll-angle constant');
disp('   is also investigated.');
disp('   ');
disp('   This program produces Figures 8-18 and 8-19 of the lecture');
disp('   notes: Aircraft Responses to Atmospheric Turbulence.');

% COMPUTE SYSTEM DYNAMICS
% NOTE: A  represents the UNCONTROLLED SYSTEM
%       A2 represents the CONTROLLED SYSTEM with gain as in lecture notes
dynamics;

% FREQUENCY VECTOR
omega = logspace(-2,2,200);

% DEFINE C AND D MATRICES
%C = [0 1 0 0 0 0 0 0 0 0];
%D = [0 0 0 0 0];

% CALCULATION OF THE POWER SPECTRAL DENSITY FUNCTION OF PHI DUE TO
% VERTICAL TURBULENCE

% Uncontrolled aircraft
[numphi,den] = ss2tf(A,B,C,D,4);
[mag,phase]  = bode(numphi,den,omega);
Sphi         = mag.*mag;

% Controlled aircraft
[numphi,den] = ss2tf(A2,B,C,D,4);
[mag,phase]  = bode(numphi,den,omega);
Sphit        = mag.*mag;

clf
axis('square')
loglog(omega,Sphi,'-',omega,Sphit,'--');
xlabel('omega [rad/s]'); ylabel('Sphi [rad2]');
%print -depsc2 -r1200 fig8_18
pause

% CALCULATION OF THE COVARIANCE MATRIX
dt   = 0.05; T = 300; t = [0:dt:T]; N = length(t);

Wdis = 1/dt;   % discrete-time intensity of white noise

% Calculation of discrete time system matrices for the 
% uncontrolled and controlled aircraft
[Phi,Gamma]   = c2d(A,B(:,4),dt);
[Phit,Gammat] = c2d(A2,B(:,4),dt);

% Initial covariance matrices
Cxx  = zeros(10,10);
Cxxt = zeros(10,10);

% Store only variance of phi
for k=2:N
   Cxx       = Phi*Cxx*Phi'    + Gamma*Wdis*Gamma';
   Cxxt      = Phit*Cxxt*Phit' + Gammat*Wdis*Gammat';
   Cx2x2(k)  = Cxx(2,2);
   Cx2x2t(k) = Cxxt(2,2);
end

% PLOT RESULTS
clf
subplot(2,1,1)
plot(t,Cx2x2);  xlabel('time [s]'); ylabel('variance phi');
title('uncontrolled aircraft');

subplot(2,1,2)
plot(t,Cx2x2t); xlabel('time [s]'); ylabel('variance phi');
title('controlled aircraft');

%print -depsc2 -r1200 fig8_19
pause

% CHECK: CALCULATION OF THE VARIANCE BY NUMERICALLY INTEGRATING THE 
%        POWER SPECTRAL DENSITY FUNCTION OF THE ROLL ANGLE PHI.
dw      = diff(omega);
varphi  = sum(Sphi(1:length(dw) )'.*dw)/pi
varphit = sum(Sphit(1:length(dw))'.*dw)/pi
