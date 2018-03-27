% Exampl42.m
%
% Chapter 4 of lecture notes ae4-304
%
% Digital calculation of covariance function and
% auto power spectral density.

% Program revised 1995, February 2004 [MM]

clf
clear

disp('   Example 4.2');
disp('   Digital calculation of covariance function and');
disp('   auto power spectral density of a stochastic signal x(t).');
disp('   ');
disp('   This program produces figures 4-X to 4-Y of the lecture');
disp('   notes: Aircraft Responses to Atmospheric Turbulence.');

fs    = input('   Enter sample frequency f (Hz)   : ');
T     = input('   Enter time T_end    (s)         : ');
dt=1/fs; N=T*fs; t=[dt:dt:T];

x     = input('   Enter function definition x(t)= : ');
int   = input('   Enter noise intensity           : ');

% generate DT white noise sequence
x  = x + sqrt(int)*randn(1,N);

% compute the auto covariance function
c  = xcov(x,'biased');         

% compute the Fast Fourier Transform
X  = fft(x);

% calculate the power spectral density
Sxx=X.*conj(X)/N;              

% define the frequency axis
f=fs*(0:N/2-1)/N;
% properly assign positive and negative frequencies
for i=1:N/2
    fg(i)     =  -f(N/2+1-i);
    Sx(i)     = Sxx(N/2+1-i); 
end
f   = [fg f];
Sxx = [Sx Sxx(1:N/2)];

% reference for white noise PSD
Refw= int*ones(1,N/2);

% PLOTTING THE RESULTS
clf
subplot(2,1,1);
plot(t,x);
xlabel('time [s]');
ylabel('x (t)');
title('SIGNAL');

subplot(2,1,2);
plot(t-T/2,c(2*N/4+1:6*N/4));
xlabel('tau [s]');
ylabel('Cxx (tau)');
title('Auto Covariance Function');
pause

clf
subplot(2,1,1);
plot(f,Sxx);
xlabel('frequency [Hz]'); 
ylabel('Sxx');
title('Auto Power Spectral Density (linear scales)');

subplot(2,1,2);
loglog(f(N/2:N-1),Sxx(N/2:N-1),f(N/2:N-1),Refw,'--');
xlabel('frequency [Hz]'); 
ylabel('Sxx');
title('Auto Power Spectral Density (loglog Scales)');

% EOF