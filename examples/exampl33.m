% Exampl33.m
%
% Chapter 3
%
% Calculates power spectral densities of displacement and
% acceleration as a result of runway surface irregularities.
% The landing gear is modelled as a second order 1-DOF
% mass-spring-damper system, the input power spectral
% density is taken from AGARD-R-632.
% The output power spectral density is calculated with
%
%          Syy(w) = |H(w)|^2 * Suu(w)
%
% while the power spectral density of the acceleration is found with
%
%          Saa(w) = w^4 Syy(w)
%
% Revised August 1992, February 2004 [MM]

clc, clf, clear

disp('   Example 3.3');
disp('   Calculates power spectral densities of displacement and');
disp('   acceleration as a result of runway surface irregularities.');
disp('   The landing gear is modelled as a second order 1-DOF');
disp('   mass-spring-damper system, the input power spectral');
disp('   density is taken from AGARD-R-632.');
disp('   The output power spectral density is calculated with');
disp('   ');
disp('            Syy(w) = |H(w)|^2 * Suu(w)');
disp('   ');
disp('   while the power spectral density of the acceleration is');
disp('   found with:');
disp('   ');
disp('            Saa(w) = w^4 Syy(w)');
disp('   ');
disp('   This program produces Figure 3.14 of the lecture notes:');
disp('   Aircraft Responses to Atmospheric Turbulence.')

omega=logspace(-2,2,100);	% frequency axis

% INPUT POWER SPECTRAL DENSITY OF 'RUNWAY RUMBLE'

% Ground speed. Can be changed.
V = input('   Enter ground speed V [m/s]  : ');

tau1=0.4/V; tau2=7/V;
for i=1:100
    Suu(i)=((6.3e-4)/V)/((1+(tau1*omega(i))^2)*(1+(tau2*omega(i))^2));
end

% LANDING GEAR DATA. TAKEN FROM DHC-2 'Beaver'
%
% aircraft mass. Can be changed.
m = input('   Enter aircraft mass m [kg]  : ');

c=20000;	% damping constant. Can be changed.
k=183000;	% spring constant. Can be changed.

% CALCULATION OF FREQUENCY RESPONSE USING TRANSFER FUNCTION H(s)
num=[c k]; den=[m c k];
h=freqs(num,den,omega);

% CALCULATION OF CROSS- AND AUTO POWER SPECTRAL DENSITIES
Suy = h.*Suu;
Syy = (abs(h).^2).*Suu;
Saa = (omega.^4).*Syy;

% PLOTTING RESULTS
clf
subplot(2,2,1);
loglog(omega,real(Suy));
xlabel('frequency, rad/s');
ylabel('Re (Suy)');
axis(10.^[-2 2 -12 0]); axis(axis);
grid
title('Cross P.S.D.')

subplot(2,2,2);
semilogx(omega,imag(Suy));
xlabel('frequency, rad/s');
ylabel('Im (Suy)');
grid
axis([10^(-2) 10^2 -2e-6 0]); axis(axis);

subplot(2,2,3);
loglog(omega,Syy);
xlabel('frequency, rad/s');
ylabel('Syy, m^2/(rad/s)');
grid
axis(10.^[-2,2,-12,0]); axis(axis);
title('Auto P.S.D. input')

subplot(2,2,4);
loglog(omega,Saa);
xlabel('frequency, rad/s');
ylabel('Saa, (m^2/s^4)/(rad/s)');
grid
axis(10.^[-2 2 -13 2]); axis(axis);
title('Auto P.S.D. output');

% EOF