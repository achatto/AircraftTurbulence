% Filename : exampl83.m
%
% Computation of analytical PSDs and time-simulation of aircraft 
% asymmetric response to atmospheric turbulence.

clc, clf, clear

disp('   Example 8.3');
disp('   ');
disp('   This example compares analytically obtained auto power');
disp('   spectral densities with experimentally obtained');
disp('   periodograms of the asymmetric motion variables.');
disp('   ');
disp('   This program produces Figures 8-20 and 8-21 of the lecture');
disp('   notes: Aircraft Responses to Atmospheric Turbulence.');
disp('   ');

% GET AIRCRAFT DYNAMICS
cit2a

% NOTE (see also cit2a.m) 
%
% THE CESSNA CITATION CE-500 IS NOT STABLE IN SPIRAL MODE (FOR THE cit2a.m 
% FLIGHT CONDITION), HENCE THE FEEDBACK CONTROLLER FOR PHI IS USED AS IN : 
%
%   delta_a = K_phi*phi (K_phi for THIS flight condition)
%
% THEREFORE, CONTROLLED AIRCRAFT SYSTEM MATRICES WILL BE USED FOR RESULTS;
%
%      A = A2
   
% DEFINE C and D MATRICES
C = [1 0 0 0 0 0 0 0 0 0   % beta
     0 1 0 0 0 0 0 0 0 0   % phi
     0 0 1 0 0 0 0 0 0 0   % pb/2V
     0 0 0 1 0 0 0 0 0 0   % rb/2V
     0 0 0 0 0 0 0 0 1 0]; % betag

D = [0 0 0 0 0
     0 0 0 0 0
     0 0 0 0 0
     0 0 0 0 0
     0 0 0 0 0];

% DEFINE FREQUENCY VECTOR
w = logspace(-2,2,300);

% COMPUTE ANALYTIC POWER SPECTRAL DENSITIES
% RESPONSE TO HORIZONTAL LATERAL TURBULENCE
temp = bode(A2,B,C(1,:),D(1,:),5,w); Sbeta  = temp.*temp;
temp = bode(A2,B,C(2,:),D(2,:),5,w); Sphi   = temp.*temp;
temp = bode(A2,B,C(3,:),D(3,:),5,w); Spp    = temp.*temp;
temp = bode(A2,B,C(4,:),D(4,:),5,w); Srr    = temp.*temp;
temp = bode(A2,B,C(5,:),D(5,:),5,w); Sbetag = temp.*temp;

Sxx  = [Sbeta Sphi Spp Srr Sbetag];

% COMPUTE PSDS USING TIME DOMAIN DATA

% SET TIME AXIS
dt = 0.01; T = 60; t = [0:dt:T]; N = length(t);

% In this case responses to lateral gust vg are calculated (fifth input):
% no asymmetric vertical and longitudinal turbulence: u_g = w_g = 0.
v_g = randn(N,1)/sqrt(dt);    % sqrt(dt) because of lsim

nn = zeros(N,1);
u  = [nn nn nn nn v_g];

% COMPUTE SYSTEM RESPONSE
y     = lsim(A2,B,C,D,u,t);

beta  = y(:,1);
phi   = y(:,2);
pbtV  = y(:,3);
rbtV  = y(:,4);
betag = y(:,5);

% PLOT TIME RESPONSES
clf
subplot(2,2,1); plot(t,beta); xlabel('time, s'); ylabel('beta');
subplot(2,2,2); plot(t,phi);  xlabel('time, s'); ylabel('phi');
subplot(2,2,3); plot(t,pbtV); xlabel('time, s'); ylabel('pb/2V');
subplot(2,2,4); plot(t,rbtV); xlabel('time, s'); ylabel('rb/2V');
%print -depsc2 -r1200 fig8_20a
pause

clf
plot(t,betag); xlabel('time, s'); ylabel('betag');
%print -depsc2 -r1200 fig8_20b
pause

% COMPUTE PERIODOGRAM AND ESTIMATE PSD
% PERIODOGRAM
BETA  = dt*fft(beta);
PHI   = dt*fft(phi);
P     = dt*fft(pbtV);
R     = dt*fft(rbtV);
BETAg = dt*fft(betag);

% PSD ESTIMATE
Pbeta  = (1/T)*( BETA.*conj(BETA));
Pphi   = (1/T)*(  PHI.*conj(PHI));
Pp     = (1/T)*(    P.*conj(P));
Pr     = (1/T)*(    R.*conj(R));
Pbetag = (1/T)*(BETAg.*conj(BETAg));

% DEFINE FREQUENCY VECTOR
fs = 1/dt;     % sample frequency
omega = 2*pi*fs*(0:(N/2)-1)/N;

% PLOT ANALYTIC AND ESTIMATED PSDS IN ONE PLOT
clf
subplot(2,2,1); loglog(w,Sxx(:,1),'--',omega,Pbeta(1:N/2)); 
axis(10.^[-2 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Sbeta');
subplot(2,2,2); loglog(w,Sxx(:,2),'--',omega,Pphi(1:N/2));
axis(10.^[-2 2 -12 -2]); xlabel('omega [rad/s]'); ylabel('Sphi')
subplot(2,2,3); loglog(w,Sxx(:,3),'--',omega,Pp(1:N/2));
axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('Spp')
subplot(2,2,4); loglog(w,Sxx(:,4),'--',omega,Pr(1:N/2));
axis(10.^[-2 2 -14 -2]); xlabel('omega [rad/s]'); ylabel('Srr')
%print -depsc2 -r1200 fig8_21a
pause

clf
loglog(w,Sxx(:,5),'--',omega,Pbetag(1:N/2)); 
xlabel('omega [rad/s]'); ylabel('Sbetag [rad2]')
%print -depsc2 -r1200 fig8_21b
pause