% Filename : exampl72.m
%
% Calculation of covariance matrix of the (non-dimensional)
% motion variables.
%
% Steady-state (co)variances are calculated by solving the
% Lyapunov-equation AC+CA'+BWB'=0 or by integrating the power
% spectral density.
%
% Transient behaviour may be calculated with either the recursive
% relation C(k+1)=PHI C(k) PHI'+ GAMMA Wdis GAMMA' or using the
% impulse response method.

clc, clf, clear

disp('   Example 7.2');
disp('   ');
disp('   Calculation of covariance matrix of the (non-dimensional)');
disp('   motion variables. Steady-state (co)variances are calculated');
disp('   by solving the Lyapunov-equation AC+CA` +BWB` =0 or by');
disp('   integrating the power spectral density. Transient behaviour');
disp('   may be calculated with either the recursive relation');
disp('   C(k+1)=PHI C(k) PHI` + GAMMA Wdis GAMMA`  or using the');
disp('   impulse response method.');
disp('   ');
disp('   This program produces Figures 7-7 to 7-11 in the lecture');
disp('   notes: Aircraft Responses to Atmospheric Turbulence.');

% COMPUTE A/C DYNAMICS
cit2s            % loading A,B matrices

% GET INPUT PARAMETERS
Wc = input('   Give noise intensity                     (  1.0 ) : ');
dt = input('   Give sampling time interval dt           (  0.01) : ');
T  = input('   Give total time interval T_end           (150.0 ) : ');
Nf = input('   Give number of points in frequency axis  ( 200  ) : ');

% DEFINE NOISE INTENSITY
W  = Wc/dt;    % discrete time covariance, remember?

% DEFINE TIME AXIS
t  = [0:dt:T]; N = length(t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE THE STEADY-STATE SOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 1. SOLVING THE LYAPUNOV-EQUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bin = B(:,3);
L   = lyap(A,Bin*Wc*Bin');
% take only the part that belongs to the 4 states
L   = L(1:4,1:4);

disp('  ');
disp('  ');
disp('   Display solution Lyapunov equation');
disp(L)
pause
disp('   Note: the diagonal elements are the variances:');
disp(diag(L)')
pause

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 2. INTEGRATING THE ANALYTICAL POWER SPECTRAL DENSITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The motion variables u/V, alpha, theta and qc/V are the first four
% elements of the seventh order state vector.

% Input vector u=[0 0 w3]'.
D = zeros(4,3);
C = [eye(4) zeros(4,3)];

% DEFINE FREQUENCY AXIS
omega = logspace(-2,2,Nf);

% COMPUTE FREQUENCY RESPONSE
mag = bode(A,B,C,D,3,omega);

% COMPUTE POWER SPECTRA
Suu = mag(:,1).^2;
Saa = mag(:,2).^2;
Stt = mag(:,3).^2;
Sqq = mag(:,4).^2;

% PLOT POWER SPECTRA
clf
subplot(2,2,1); 
loglog(omega,Suu); xlabel('omega [rad/sec]'); ylabel('Suu [rad^2]');
subplot(2,2,2); 
loglog(omega,Saa); xlabel('omega [rad/sec]'); ylabel('Saa [rad^2]');
subplot(2,2,3); 
loglog(omega,Stt); xlabel('omega [rad/sec]'); ylabel('Stt [rad^2]');
subplot(2,2,4); 
loglog(omega,Sqq); xlabel('omega [rad/sec]'); ylabel('Sqq [rad^2]');
pause

% NUMERICAL INTEGRATION OF PSD's
do = diff(omega)';   % compute "difference vector" in omega
                     % i.e., omega(k+1)-omega(k);
% then perform (very crude) integration
var(1) = sum(do.*Suu(1:Nf-1));
var(2) = sum(do.*Saa(1:Nf-1));
var(3) = sum(do.*Stt(1:Nf-1));
var(4) = sum(do.*Sqq(1:Nf-1));

disp('   Numerical integration of PSD yields the variances:')
var = var/pi
disp('   Note that when we have more frequency points the')
disp('   integration will become more accurate and the')
disp('   variances will approximate the Lyapunov solution.')
pause

hh=ones(1,N);
var1 = var(1)*hh; var2 = var(2)*hh; var3 = var(3)*hh; var4 = var(4)*hh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALCULATION OF GROWTH IN TIME OF COVARIANCE MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 1. BY RECURSIVE CALCULATION 
%
%    Cxx(k+1) = PHI Cxx(k) PHI' + GAMMA Wdis GAMMA'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Discretize the system matrices
[Phi,Gamma] = c2d(A,B(:,3),dt);

% Initial conditions
Cxx = zeros(7,7);

% Discrete Solution of response equation
disp('   Computing response by recursive calculation')
hh=2;
for k=2:N
    Cxx=Phi*Cxx*Phi' + Gamma*W*Gamma';
    Cx1x1(k)=Cxx(1,1);Cx1x2(k)=Cxx(1,2);Cx1x3(k)=Cxx(1,3);Cx1x4(k)=Cxx(1,4);
                      Cx2x2(k)=Cxx(2,2);Cx2x3(k)=Cxx(2,3);Cx2x4(k)=Cxx(2,4);
                                        Cx3x3(k)=Cxx(3,3);Cx3x4(k)=Cxx(3,4);
                                                          Cx4x4(k)=Cxx(4,4);
    if hh > 100, hh=1; fprintf('step %g\n',k); end; hh=hh+1; 
end
fprintf('ready')

I=ones(1,N);
Css11=L(1,1)*I;Css12=L(1,2)*I;Css13=L(1,3)*I;Css14=L(1,4)*I;
               Css22=L(2,2)*I;Css23=L(2,3)*I;Css24=L(2,4)*I;
                              Css33=L(3,3)*I;Css34=L(3,4)*I;
                                             Css44=L(4,4)*I;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Method 2. USING THE IMPULSE RESPONSE METHOD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ZERO INPUT, INITIAL CONDITION EQUALS B (for input 3 in this case)
u = zeros(3,N); x0=B(:,3);

% CALCULATION OF IMPULSE RESPONSES
h = lsim(A,B,C,D,u,t,x0);	

% PLOT IMPULSE RESPONSE
clf
subplot(2,1,1)
plot(t,h(:,1)); xlabel('Time [sec]'); ylabel('h u w3(t)');
subplot(2,1,2)
plot(t,h(:,2)); xlabel('Time [sec]'); ylabel('h alpha w3(t)');
pause

clf
subplot(2,1,1)
plot(t,h(:,3)); xlabel('Time [sec]'); ylabel('h theta w3(t)');
subplot(2,1,2)
plot(t,h(:,4)); xlabel('Time [sec]'); ylabel('h qc/V w3(t)');
pause

% CALCULATION OF PRODUCT MATRIX OF IMPULSE RESPONSES
h11=h(:,1).*h(:,1);h12=h(:,1).*h(:,2);h13=h(:,1).*h(:,3);h14=h(:,1).*h(:,4);
                   h22=h(:,2).*h(:,2);h23=h(:,2).*h(:,3);h24=h(:,2).*h(:,4);
                                      h33=h(:,3).*h(:,3);h34=h(:,3).*h(:,4);
                                                         h44=h(:,4).*h(:,4);

% PLOT (CROSS) PRODUCTS OF IMPULSE RESPONSES
clf
subplot(2,2,1)
plot(t,h11); xlabel('Time [sec]'); ylabel('h1*h1(t)');
subplot(2,2,2)
plot(t,h12); xlabel('Time [sec]'); ylabel('h1*h2(t)');
subplot(2,2,4)
plot(t,h22); xlabel('Time [sec]'); ylabel('h2*h2(t)');
pause

clf
subplot(2,2,1)
plot(t,h13); xlabel('Time [sec]'); ylabel('h1*h3(t)');
subplot(2,2,2)
plot(t,h14); xlabel('Time [sec]'); ylabel('h1*h4(t)');
subplot(2,2,3)
plot(t,h23); xlabel('Time [sec]'); ylabel('h2*h3(t)');
subplot(2,2,4)
plot(t,h24); xlabel('Time [sec]'); ylabel('h2*h4(t)');
pause

clf
subplot(2,2,1)
plot(t,h33); xlabel('Time [sec]'); ylabel('h3*h3(t)');
subplot(2,2,2)
plot(t,h34); xlabel('Time [sec]'); ylabel('h3*h4(t)');
subplot(2,2,4)
plot(t,h44); xlabel('Time [sec]'); ylabel('h4*h4(t)');
pause

% INTEGRATION OF PRODUCT MATRIX OF IMPULSE RESPONSES
var11(1)=0; var12(1)=0; var13(1)=0; var14(1)=0;
            var22(1)=0; var23(1)=0; var24(1)=0;
                        var33(1)=0; var34(1)=0;
                                    var44(1)=0;
                                    
dth11 = dt*h11; dth12 = dt*h12; dth13 = dt*h13; dth14 = dt*h14;
                dth22 = dt*h22; dth23 = dt*h23; dth24 = dt*h24;
                                dth33 = dt*h33; dth34 = dt*h34;
                                                dth44 = dt*h44;
for i=1:N-1
    var11(i+1) = var11(i) + dth11(i);
    var12(i+1) = var12(i) + dth12(i);
    var13(i+1) = var13(i) + dth13(i);
    var14(i+1) = var14(i) + dth14(i);
    var22(i+1) = var22(i) + dth22(i);
    var23(i+1) = var23(i) + dth23(i);
    var24(i+1) = var24(i) + dth24(i);
    var33(i+1) = var33(i) + dth33(i);
    var34(i+1) = var34(i) + dth34(i);
    var44(i+1) = var44(i) + dth44(i);
end

% PLOT VARIANCES FROM IMPULSE RESPONSE METHOD
clf
subplot(2,2,1)
plot(t,Css11,'--',t,var11); xlabel('time [s]'); ylabel('Cx1x1');
subplot(2,2,2)
plot(t,Css12,'--',t,var12); xlabel('time [s]'); ylabel('Cx1x2');
subplot(2,2,4)
plot(t,Css22,'--',t,var22); xlabel('time [s]'); ylabel('Cx2x2');
pause

clf
subplot(2,2,1)
plot(t,Css13,'--',t,var13); xlabel('time [s]'); ylabel('Cx1x3');
subplot(2,2,2)
plot(t,Css14,'--',t,var14); xlabel('time [s]'); ylabel('Cx1x4');
subplot(2,2,3)
plot(t,Css23,'--',t,var23); xlabel('time [s]'); ylabel('Cx2x3');
subplot(2,2,4)
plot(t,Css24,'--',t,var24); xlabel('time [s]'); ylabel('Cx2x4');
pause

clf
subplot(2,2,1)
plot(t,Css33,'--',t,var33); xlabel('time [s]'); ylabel('Cx3x3');
subplot(2,2,2)
plot(t,Css34,'--',t,var34); xlabel('time [s]'); ylabel('Cx3x4');
subplot(2,2,4)
plot(t,Css44,'--',t,var44); xlabel('time [s]'); ylabel('Cx4x4');
pause

% PLOT RESULTS FROM RECURSIVE EQUATION WITH STEADY-STATE DASHED
clf
subplot(2,2,1)
plot(t,Cx1x1,t,Css11,'--'); xlabel('time [s]'); ylabel('Cx1x1');
subplot(2,2,2)
plot(t,Cx1x2,t,Css12,'--'); xlabel('time [s]'); ylabel('Cx1x2');
subplot(2,2,4)
plot(t,Cx2x2,t,Css22,'--'); xlabel('time [s]'); ylabel('Cx2x2')
pause

clf
subplot(2,2,1)
plot(t,Cx1x3,t,Css13,'--'); xlabel('time [s]'); ylabel('Cx1x3');
subplot(2,2,2)
plot(t,Cx1x4,t,Css14,'--'); xlabel('time [s]'); ylabel('Cx1x4');
subplot(2,2,3)
plot(t,Cx2x3,t,Css23,'--'); xlabel('time [s]'); ylabel('Cx2x3');
subplot(2,2,4)
plot(t,Cx2x4,t,Css24,'--'); xlabel('time [s]'); ylabel('Cx2x4');
pause

clf
subplot(2,2,1)
plot(t,Cx3x3,t,Css33,'--'); xlabel('time [s]'); ylabel('Cx3x3');
subplot(2,2,2)
plot(t,Cx3x4,t,Css34,'--'); xlabel('time [s]'); ylabel('Cx3x4');
subplot(2,2,4)
plot(t,Cx4x4,t,Css44,'--'); xlabel('time [s]'); ylabel('Cx4x4');
pause
