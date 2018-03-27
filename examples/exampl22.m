% Exampl22   Calculates the signal, the auto product function, the
%            auto-covariance function and the auto correlation function
%            of an arbitrary signal.
%
% Chapter 2 of the lecture notes ae4-304.
% 
% Program revised August 1992, February 2004 [MM]

clc
clf
clear

disp('   Example 2.2                                          ');
disp('                                                        ');
disp('   Calculates signal, auto product-, -covariance- and   ');
disp('   -correlation function of a stochastic signal x(t).   ');
disp('                                                        ');
disp('   This program can produce Figures 2.13 and 2.14 of the');
disp('   lecture notes ae4-304.                               ');

% create time axis
T     = input('   Enter total time interval T       : ');
dt    = input('   Enter sampling time interval dt   : ');

t = [-T/2:dt:T/2];
N = length(t);

% create signal
x     = input('   Enter function definition x(t) =  : ');
int   = input('   Enter noise intensity             : ');

x = x + sqrt(int)*randn(1,N);

fprintf('\n\n'                                   );
fprintf('   Signal mean     = %f\n',mean(x)      );
fprintf('   Signal variance = %f\n',cov(x)       );
fprintf('   Signal std.dev. = %f\n',sqrt(cov(x)) );
pause

% compute variables of interest
r=xcorr(x,'unbiased');      % auto product function
c=xcov (x,'unbiased');      % auto covariance function
k=c/(std(x)^2);             % auto correlation function

% plot results
clf
subplot(2,1,1);
plot(t+T/2,x);
xlabel('time'); ylabel('x (t)'); title('signal');

subplot(2,1,2);
plot(t,r((T/2)/dt:3*(T/2)/dt));
xlabel('tau'); ylabel('Rxx (tau)'); title('auto product function');
pause;

clf
subplot(2,1,1)
plot(t,c((T/2)/dt:3*(T/2)/dt));
xlabel('tau'); ylabel('Cxx (tau)'); title('auto covariance function');

subplot(2,1,2)
plot(t,k((T/2)/dt:3*(T/2)/dt));
xlabel('tau'); ylabel('Kxx (tau)'); title('auto correlation function');

% EOF