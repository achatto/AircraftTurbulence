% Filename : exampl74b.m
%
% Calculation of the output spectral densities using the
% MATLAB FFT algorithm

clc, clf, clear

disp('   Example 7.4');
disp('   ')
disp('   This program calculates the periodograms of time-domain');
disp('   data using the MATLAB FFT algorithm. The experimentally');
disp('   obtained PSD will be compared with the analytically derived');
disp('   PSD of the motion variables.');
disp('   ');
disp('   This program produces Figures 7-14 to 7-17 of the lecture');
disp('   notes: Aircraft Responses to Atmospheric Turbulence.');
disp('   ')

% RUN EXAMPL74A 
exampl74a

% DEFINE TIME AXIS
dt = 0.05; fs = 1/dt;
T  = 200;  t = [0:dt:T]; N = length(t);

% CREATE ZERO INPUT SIGNAL
delta = zeros(1,N);

% CREATE NORMAL WHITE NOISE SIGNALS 
%  if Vertical turbulence   (u=3) : w1 = 0
%  if Horizontal turbulence (u=2) : w3 = 0
if u == 2
   w1 = randn(1,N)/sqrt(dt);   % sqrt(dt) because of lsim 
   w3 = zeros(1,N);
else
   w1 = zeros(1,N);
   w3 = randn(1,N)/sqrt(dt);   % sqrt(dt) because of lsim
end

inpsig  = [delta' w1' w3'];

% DEFINE C and D MATRICES
C = [1 0 0 0 0 0 0   % u/V
     0 1 0 0 0 0 0   % alpha
     0 0 1 0 0 0 0   % theta
     0 0 0 1 0 0 0   % qc/V
     0 0 0 0 1 0 0   % u_g/V
     0 0 0 0 0 1 0]; % alpha_g
D = zeros(6,3);

% COMPUTE TIME RESPONSE
y      = lsim(A,B,C,D,inpsig,t);

hatu   = y(:,1);
alpha  = y(:,2);
theta  = y(:,3);
qcV    = y(:,4);
hatug  = y(:,5);
alphag = y(:,6);

% Add a trailing zero for alpha array. Because we use the routine
% diff(w) which fills a vector of length(w)-1 with w(i+1)-w(i) for i=1 to
% length(w).

alphanz                    = alpha;
alphanz(length(alphanz)+1) = 0;

% Calculation of the normal load factor nz according to: 
nz = (V/g)*( (V/c)*qcV-diff(alphanz)/dt );
nz(length(nz))=nz(length(nz)-1);

% PLOT TIME RESPONSES
clf
subplot(3,2,1);
plot(t,hatu);  xlabel('Time [s]'); ylabel('u/V [rad]');
subplot(3,2,2);
plot(t,alpha); xlabel('Time [s]'); ylabel('alpha [rad]');
subplot(3,2,3);
plot(t,theta); xlabel('Time [s]'); ylabel('theta [rad]');
subplot(3,2,4);
plot(t,qcV);   xlabel('Time [s]'); ylabel('qc/V [rad]');
subplot(3,2,5);
plot(t,nz);    xlabel('Time [s]'); ylabel('nz');
subplot(3,2,6);
if u == 2
    plot(t,hatug); xlabel('Time [s]'); ylabel('ug/V [rad]');
    print -depsc2 -r1200 fig7_14
else
    plot(t,alphag); xlabel('Time [s]'); ylabel('alphag [rad]');
    print -depsc2 -r1200 fig7_16
end    
pause

% FFT ALL SIGNALS
U      = dt*fft(hatu);
ALPHA  = dt*fft(alpha);
THETA  = dt*fft(theta);
QCV    = dt*fft(qcV);
NZ     = dt*fft(nz);
Ug     = dt*fft(hatug);
ALPHAg = dt*fft(alphag);

% COMPUTE PSDs
Pu      = (1/T)*     U.*conj(U);
Palpha  = (1/T)* ALPHA.*conj(ALPHA);
Ptheta  = (1/T)* THETA.*conj(THETA);
PqcV    = (1/T)*   QCV.*conj(QCV);
Pnz     = (1/T)*    NZ.*conj(NZ);
Pug     = (1/T)*    Ug.*conj(Ug);
Palphag = (1/T)*ALPHAg.*conj(ALPHAg);

% DEFINE FREQUENCY VECTOR FOR PLOTTING
omega = 2*pi*fs*(0:(N/2)-1)/N;

% PLOT PSDs
clf
subplot(3,2,1), loglog(w,Sxx(:,1),'--',omega,Pu(1:N/2));
axis(10.^[-2,2,-15,0])
xlabel('omega [rad/s]'); ylabel('Suu [rad2]')

subplot(3,2,2), loglog(w,Sxx(:,2),'--',omega,Palpha(1:N/2));
axis(10.^[-2,2,-15,0])
xlabel('omega [rad/s]'); ylabel('Saa [rad2]')

subplot(3,2,3), loglog(w,Sxx(:,3),'--',omega,Ptheta(1:N/2));
axis(10.^[-2,2,-15,0])
xlabel('omega [rad/s]'); ylabel('Stt [rad2]')

subplot(3,2,4), loglog(w,Sxx(:,4),'--',omega,PqcV(1:N/2));
axis(10.^[-2,2,-20,-5])
xlabel('omega [rad/s]'); ylabel('Sqq [rad2]')

subplot(3,2,5), loglog(w,Sxx(:,5),'--',omega,Pnz(1:N/2));
axis(10.^[-2,2,-10,0])
xlabel('omega [rad/s]'); ylabel('Snznz [rad2]')

if u == 2
   subplot(3,2,6), loglog(w,Sxx(:,6),'--',omega,Pug(1:N/2));
   axis(10.^[-2,2,-15,0])
   xlabel('omega [rad/s]'); ylabel('Sugug [rad2]')
   print -depsc2 -r1200 fig7_15
else
   subplot(3,2,6), loglog(w,Sxx(:,6),'--',omega,Palphag(1:N/2));
   axis(10.^[-2,2,-15,0])
   xlabel('omega [rad/s]'); ylabel('Sagag [rad2]')
   print -depsc2 -r1200 fig7_17
end    
