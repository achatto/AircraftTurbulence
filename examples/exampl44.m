% Exampl44.m
%
% Chapter 4 of the lecture notes ae4-304
%
% Calculates digitally the spectral densities of displacement
% and acceleration as a result of runway surface irregularities.
% The results are compared with the analytical solutions as
% found in example 3.4.

% Program revised August 1992, February 2004 [MM]

clf
clear

disp('   Example 4.4');
disp('   Compares the digitally calculated spectral densities of');
disp('   displacement and acceleration as a result of runway surface');
disp('   irregularities with the analytical solutions found in');
disp('   Example 3.4. The landing gear is modelled as a second order');
disp('   1-DOF mass-spring-damper system, the input power spectral');
disp('   density is taken from AGARD-R-632.');
disp('   The output power spectral density is calculated with:');
disp('   ');
disp('          Syy(w) = |H(w)|^2 * Suu(w)       (analytical)');
disp('          Syy(w) = spectrum(Y,Ns);         (numerical)');
disp('   ');
disp('   while the power spectral density of the acceleration');
disp('   is found with');
disp('   ');
disp('          Saa(w) = w^4 Syy(w)');
disp('   ');
disp('   This program produces figure 4-13 of the lecture notes:');
disp('   Aircraft Responses to Atmospheric Turbulence.');

V  = 10;%    input('   Enter ground speed V [m/s]           : ');
fs = 10;%    input('   Enter sample frequency [Hz]          : ');
N  = 1024;%    input('   Enter number of samples              : ');
Ns =  64;%  input('   Enter number of samples in spectrum  : ');
Wn =   1;%  input('   Enter noise intensity                : ');

dt=1/fs; T=N*dt; t=[0:dt:T-dt];

w=sqrt(Wn)*randn(1,N);

tau1=0.4/V; tau2=7/V;

% LANDING GEAR PARAMETERS (Can be changed)
% aircraft mass. Can be changed.
%
m  =  2290;%   input('   Enter aircraft mass m [kg]           : ');

c=20000;         % damping constant. Can be changed.
k=183000;        % spring constant. Can be changed.

% DIGITAL CALCULATION OF INPUT POWER SPECTRAL DENSITIES USING FFT

num=(6.3e-4)/V;	 % 'rumble' filter definition;
den=[tau1*tau2 tau1+tau2 1];
sys=tf(sqrt(num),den);

u=lsim(sys,w,t);

Swwspct=spectrum(w,Ns);

omega=2*pi*fs*(1:N/Ns:N/2)/N;

% ANALYTICAL EXPRESSIONS FOR INPUT POWER SPECTRAL DENSITIES
for i=1:Ns/2
  Swwanal(i,1) = Wn;
  Suuanal(i,1) = num/((1+tau1^2*omega(i)^2)*(1+tau2^2*omega(i)^2))*Swwanal(i);
end

% DIGITAL CALCULATION OF CROSS- AND OUTPUT SPECTRA USING FFT
num = [c k]; den=[m c k];		% transfer function H(s)
sys = tf(num,den);
y   = lsim(sys,u,t);		% time history y using u

Suyspct = spectrum(u,y,Ns);
Saaspct = (omega.^4).*Suyspct(:,2)';

% ANALYTICAL EXPRESSIONS FOR OUTPUT POWER SPECTRAL DENSITIES
h = freqs(num,den,omega);

maganal   = abs(h); 
phaseanal = 180*angle(h)/pi;

Suyanal = h.*Suuanal;
Syyanal = (maganal.^2).*Suuanal;
Saaanal = (omega.^4).*Syyanal;

% ESTIMATION OF FREQUENCY RESPONSE FUNCTION FROM POWER SPECTRAL DENSITY RATIO
hspct =Suyspct(:,4);	%Suyspct(:,3)./Suyspct(:,1);

magspct   = abs(hspct);  
phasespct = 180*angle(hspct)/pi;

% PLOTTING THE RESULTS
clf
subplot(2,2,1)

loglog(omega/(2*pi),Swwanal,omega/(2*pi),Swwspct(:,1),'--');
xlabel('frequency [Hz]');ylabel('Sww');
title('PSD White Noise');

subplot(2,2,2)
loglog(omega/(2*pi),Suuanal,omega/(2*pi),Suyspct(:,1),'--');
xlabel('frequency [Hz]');ylabel('Suu');
title('PSD Forming Filter Output')

subplot(2,2,3)
loglog(omega/(2*pi),real(Suyanal),omega/(2*pi),real(Suyspct(:,3)),'--');
xlabel('frequency [Hz]');ylabel('Re(Suy)');
title('Cross PSD w-u');

subplot(2,2,4)
semilogx(omega/(2*pi),imag(Suyanal),omega/(2*pi),imag(Suyspct(:,3)),'--');
xlabel('frequency [Hz]');ylabel('Im(Suy)');
pause;

clg
subplot(2,2,1);
loglog(omega/(2*pi),Syyanal,omega/(2*pi),Suyspct(:,2),'--');
xlabel('frequency [Hz]');ylabel('Syy');
title('PSD Model Output');

subplot(2,2,2);
loglog(omega/(2*pi),Saaanal,omega/(2*pi),Saaspct,'--');
xlabel('frequency [Hz]');ylabel('Saa');
title('PSD Normal Acceleration');

subplot(2,2,3);
loglog(omega/(2*pi),maganal,omega/(2*pi),magspct,'--');
xlabel('frequency [Hz]');ylabel('gain (H)');
title('Frequency Response Function');

subplot(2,2,4);
semilogx(omega/(2*pi),phaseanal,omega/(2*pi),phasespct,'--');
xlabel('frequency [Hz]');ylabel('phase angle (H)');
pause
