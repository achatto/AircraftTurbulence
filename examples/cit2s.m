% Filename:  cit2s.m
%
% Calculation of state matrix and input matrix for calculation
% of symmetric aircraft response to atmospheric turbulence.
% The system model is in the form
%
%       .
%       x = Ax + Bu
%       -    -    -
%
% with 
%     x = [u/V alpha theta qc/V u_g/V alpha_g alpha_g*]' 
% and 
%     u = [delta_e w_1 w_3]'.
%
% The turbulence filters are according to Dryden.

%
%  Cessna Citation Ce-500, landing (1)
%

% INPUT TURBULENCE- AND AIRCRAFT PARAMETERS

% AIRCRAFT FLIGHT CONDITION 'LANDING'.
V     = 59.9;
m     = 4547.8;
twmuc = 2*102.7;
KY2   = 0.980;
c     = 2.022;
S     = 24.2;
lh    = 5.5;

% TURBULENCE PARAMETERS
disp(' ');
sigma = input('   Enter turbulence intensity sigma [m/s]     (  1) : ');
Lg    = input('   Enter turbulence scale length Lg [m]       (150) : ');

sigmaug_V = sigma/V;
sigmaag   = sigma/V;

% AIRCRAFT SYMMETRIC AERODYNAMIC DERIVATIVES : 
CX0 = 0.0000;     CZ0  =-1.1360;     Cm0  =  0.0000;
CXu =-0.2199;     CZu  =-2.2720;     Cmu  =  0.0000;
CXa = 0.4653;     CZa  =-5.1600;     Cma  = -0.4300;
CXq = 0.0000;     CZq  =-3.8600;     Cmq  = -7.0400;
CXd = 0.0000;     CZd  =-0.6238;     Cmd  = -1.5530;
CXfa= 0.0000;     CZfa =-1.4300;     Cmfa = -3.7000;
                  CZfug= 0.0000;     Cmfug= -Cm0*lh/c;
                  CZfag= CZfa-CZq;   Cmfag=  Cmfa-Cmq;

% CALCULATION OF AIRCRAFT SYMMETRIC STABILITY DERIVATIVES
xu   = (V/c)*(CXu/twmuc);
xa   = (V/c)*(CXa/twmuc);
xt   = (V/c)*(CZ0/twmuc);
xq   = 0;
xd   = (V/c)*(CXd/twmuc);
xug  = xu;
xfug = 0;
xag  = xa;
xfag = 0;

zu   = (V/c)*( CZu/(twmuc-CZfa));
za   = (V/c)*( CZa/(twmuc-CZfa));
zt   = (V/c)*(-CX0/(twmuc-CZfa));
zq   = (V/c)*((CZq+twmuc)/(twmuc-CZfa));
zd   = (V/c)*( CZd/(twmuc-CZfa));
zug  = zu;
zfug = (V/c)*( CZfug/(twmuc-CZfa));
zag  = za;
zfag = (V/c)*( CZfag/(twmuc-CZfa));

mu   = (V/c)*(( Cmu+CZu*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
ma   = (V/c)*(( Cma+CZa*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mt   = (V/c)*((-CX0*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mq   = (V/c)*(Cmq+Cmfa*(twmuc+CZq)/(twmuc-CZfa))/(twmuc*KY2);
md   = (V/c)*((Cmd+CZd*Cmfa/(twmuc-CZfa))/(twmuc*KY2));
mug  = mu;
mfug = (V/c)*(Cmfug+CZfug*Cmfa/(twmuc-CZfa))/(twmuc*KY2);
mag  = ma;
mfag = (V/c)*(Cmfag+CZfag*Cmfa/(twmuc-CZfa))/(twmuc*KY2);

% STATE- AND INPUT MATRICES
A=[xu xa xt 0    xug                  xag       0;
   zu za zt zq   zug-zfug*V/Lg*(c/V)  zag       zfag*(c/V);
   0  0  0  V/c  0                    0         0;
   mu ma mt mq   mug-mfug*V/Lg*(c/V)  mag       mfag*(c/V);
   0  0  0  0   -V/Lg                 0         0;
   0  0  0  0    0                    0         1;
   0  0  0  0    0                   -(V/Lg)^2 -2*V/Lg];

B=...
 [xd 0                                 0;
  zd zfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) zfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  0                                 0;
  md mfug*(c/V)*sigmaug_V*sqrt(2*V/Lg) mfag*(c/V)*sigmaag*sqrt(3*V/Lg);
  0  sigmaug_V*sqrt(2*V/Lg)            0;
  0  0                                 sigmaag*sqrt(3*V/Lg);
  0  0                                 (1-2*sqrt(3))*sigmaag*sqrt((V/Lg)^3)];

% GAIN FACTORS FOR AUTOPILOT CHAPTER 7 AND DEFINITION OF MATRIX At
% (approximately "Dead Beat" damping)
Kt = -0.21; Kq = 0;      % gains on "theta" and "q"
K  = [0 0 Kt Kq 0 0 0];  % feedback matrix
At = A-B(:,1)*K;         % new A matrix = (A - BK) because of feedback

% SAVE ALL TO FILE, THEN CLEAR, THEN READ AGAIN
save dumfile A At B sigmaug_V sigmaag Lg V c
clear
load dumfile