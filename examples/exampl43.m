% Exampl43.m
%
% Chapter 4 of lecture notes ae4-304
%
% Calculates digitally the spectral densities of displacement and
% acceleration as a result of runway surface irregularities.
% The results are compared with the analytical solutions
% as found in Example 3.4.

% Program revised August 1992, February 2004 [MM]

clc; clg; clear

disp('   Example 4.3');
disp('   Compares the digitally calculated spectral densities');
disp('   of displacement and acceleration as a result of runway');
disp('   surface irregularities with the analytical solutions');
disp('   found in Example 3.4. The landing gear is modelled as a');
disp('   second order 1-DOF mass-spring-damper system, the input');
disp('   power spectral density is taken from AGARD-R-632.');
disp('   The output power spectral density is calculated with:');
disp('   ');
disp('           Syy(w) = |H(w)|^2 * Suu(w)       (analytical)');
disp('           Syy(w) = conj(Y[k]).Y[k]/N       (numerical)');
disp('   ');
disp('   while the power spectral density for the acceleration');
disp('   is found with');
disp('   ');
disp('          Saa(w) = w^4 Syy(w)');
disp('   ');
disp('   This program produces figures 4-11 and 4-12 of the lecture');
disp('   notes: Aircraft Responses to Atmospheric Turbulence.');

V  = input('   Enter ground speed V [m/s]          : ');
fs = input('   Enter sample frequency [Hz]         : ');
N  = input('   Enter number of samples             : ');
Wn = input('   Enter noise intensity               : ');
m  = input('   Enter aircraft mass [kg]            : ');

dt = 1/fs;  % sample time
T  = N*dt;  % time interval of realization
t=[0:dt:T-dt]; % time axis

tau1=0.4/V; tau2=7/V;

% LANDING GEAR PARAMETERS (Can be changed)
%
% aircraft mass
c = 20000;       % damping constant
k = 183000;      % spring constant

% DIGITAL CALCULATION OF INPUT SIGNAL
% the shaping filter ('rumble' filter definition)
num1= (6.3e-4)/V;     
den = [tau1*tau2 tau1+tau2 1];
sys = tf(sqrt(num1),den);

% create CT white noise
w = sqrt(Wn).*randn(1,N);
% filter it before you use it as an input to Matlab's lsim!
[B,A]=butter(3,0.99); 
w = filter(B,A,w);

% compute response of the rumble filter to "white noise",
% yields the colored noise signal u
u = lsim(sys,w,t,'zoh');

% DIGITAL CALCULATION OF OUTPUT SIGNAL
% the suspension system (with known transfer function H(s))
num=[c k]; den=[m c k];
sys = tf(num,den);

% compute response of the suspension system to the
% colored noise u
y = lsim(sys,u,t,'zoh');

% plot the results
clf
subplot(3,1,1)
plot(t,w); 
vv=axis; vv(1,1)=0; vv(1,2)=T; axis(vv);
xlabel('time [sec]'); ylabel('w(t)'); title('White Noise w(t)');

subplot(3,1,2)
plot(t,u); 
vv=axis; vv(1,1)=0; vv(1,2)=T; axis(vv);
xlabel('time [sec]'); ylabel('u(t)'); title('Surface irregularity u(t)');

subplot(3,1,3)
plot(t,y);
vv=axis; vv(1,1)=0; vv(1,2)=T; axis(vv);
xlabel('time [sec]');ylabel('y(t)');
title('Suspension system output y(t)');
pause

% TRANSFORM ALL TO THE FREQUENCY DOMAIN
U = fft(u,N); 
W = fft(w,N);
Y = fft(y,N);

Suucalc = U.*conj(U)/N;
Swwcalc = W.*conj(W)/N;
Syycalc = Y.*conj(Y)/N;

Suycalc = conj(U(1:N)).*Y(1:N)/N;

omega = 2*pi*fs*(1:N/2)/N; omega=omega';

% ANALYTICAL EXPRESSIONS FOR INPUT POWER SPECTRAL DENSITIES
for i=1:N/2
    Swwanal(i,1) = Wn;
    Suuanal(i,1) = num1/((1+tau1^2*omega(i)^2)*(1+tau2^2*omega(i)^2))*Swwanal(i);
end
% compute the accelerations
Saacalc = (omega.^4).*Syycalc(1:N/2);

% ANALYTICAL EXPRESSIONS FOR OUTPUT POWER SPECTRAL DENSITIES
h = freqs(num,den,omega);
maganal   = abs(h); 
phaseanal = 180*angle(h)/pi;

Suyanal = h.*Suuanal;
Syyanal = (maganal.^2).*Suuanal;
Saaanal = (omega.^4).*Syyanal;

% ESTIMATION OF FREQUENCY RESPONSE FUNCTION FROM
% POWER SPECTRAL DENSITY RATIO
hcalc = Suycalc./Suucalc;

magcalc   = abs(hcalc);  
phasecalc = 180*angle(hcalc)/pi;

% PLOT THE RESULTS
clf
subplot(2,2,1);
loglog(omega/(2*pi),Swwanal,omega/(2*pi),Swwcalc(1:N/2),'--');
vv=axis; vv=[0.001 10 0.001 10]; axis(vv);
xlabel('frequency [Hz]');ylabel('Sww');
title('PSD White Noise');

subplot(2,2,2);
loglog(omega/(2*pi),Suuanal,omega/(2*pi),Suucalc(1:N/2),'--');
vv=axis; vv=[0.001 10 10^(-10) 0.1]; axis(vv);
xlabel('frequency [Hz]');ylabel('Suu');
title('PSD Forming Filter Output')

subplot(2,2,3);
loglog(omega/(2*pi),real(Suyanal),omega/(2*pi),real(Suycalc(1:N/2)),'--');
vv=axis; vv=[0.001 10 10^(-11) 0.01]; axis(vv);
xlabel('frequency [Hz]');ylabel('Re (Suy)');
title('Cross PSD w-u');

subplot(2,2,4);
semilogx(omega/(2*pi),imag(Suyanal),omega/(2*pi),imag(Suycalc(1:N/2)),'--');
xlabel('frequency [Hz]');ylabel('Im (Suy)');
pause

clf
subplot(2,2,1);
loglog(omega/(2*pi),Syyanal,omega/(2*pi),Syycalc(1:N/2),'--');
xlabel('frequency [Hz]');ylabel('Syy');
title('PSD Model Output');

subplot(2,2,2);
loglog(omega/(2*pi),Saaanal,omega/(2*pi),Saacalc(1:N/2),'--');
xlabel('frequency [Hz]');ylabel('Saa');
title('PSD Normal Acceleration');

subplot(2,2,3);
loglog(omega/(2*pi),maganal,omega/(2*pi),magcalc(1:N/2),'--');
xlabel('frequency [Hz]');ylabel('gain (H)');
title('Frequency Response Function');

subplot(2,2,4);
semilogx(omega/(2*pi),phaseanal,omega/(2*pi),phasecalc(1:N/2),'--');
xlabel('frequency [Hz]');ylabel('phase angle (H)');

% EOF