% Filename : exampl74a.m
%
% Calculation of analytical power spectral densities and variance of
% motion variables

clf, clc, clear

% GET AIRCRAFT DYNAMICS
cit2s;

% DEFINE MISCELLANEOUS
g  = 9.80665;

% calculation of the frequency response of the normal acceleration factor
% horizontal turbulence : u=2, vertical turbulence u=3.

disp('   Example 7.3');
disp('   ');
disp('   Input 2 for horizontal turbulence and 3 for vertical turbulence');
disp('   excitation');
disp('   ');

u = input('   Input 2 or 3 : ');

% FREQUENCY AXIS
Nomega = 300; w = logspace(-2,2,Nomega);

% C and D MATRICES
Cu     = [1 0 0 0 0 0 0];
Calpha = [0 1 0 0 0 0 0];
Ctheta = [0 0 1 0 0 0 0];
Cq     = [0 0 0 1 0 0 0];
Cug    = [0 0 0 0 1 0 0];
Cag    = [0 0 0 0 0 1 0];
D      = [0 0 0];

Calphadot = A(2,:);
Dalphadot = B(2,:);

% COMPUTE FREQUENCY RESPONSE FUNCTION AND PSD
mag = bode(A,B,Cu,D,u,w);     Suu   = mag.*mag;
mag = bode(A,B,Calpha,D,u,w); Saa   = mag.*mag;
mag = bode(A,B,Ctheta,D,u,w); Stt   = mag.*mag;
mag = bode(A,B,Cq,D,u,w);     Sqq   = mag.*mag;
mag = bode(A,B,Cug,D,2,w);    Sugug = mag.*mag;
mag = bode(A,B,Cag,D,3,w);    Sagag = mag.*mag;

% COMPUTE FREQ. RESPONSE of NZ
sys    = ss(A,B(:,u),Calphadot,Dalphadot(:,u));
Hadotw = freqresp(sys,w);
sys    = ss(A,B(:,u),Cq,D(:,u));
Hqw    = freqresp(sys,w);

Hnz   = (V/g)*((V/c)*Hqw - Hadotw);
mag   = abs(squeeze(Hnz));    
Snznz = mag.*mag;

Sxx = [Suu Saa Stt Sqq Snznz Sugug Sagag];

% COMPUTE VARIANCE THROUGH CRUDE INTEGRATION OF PSDs
var = zeros(1,7);
for i=1:Nomega-1;
  for j=1:7
    var(j)=var(j)+(w(i+1)-w(i))*Sxx(i,j);
    vart(i,j)=var(j);
  end
end

var=var/pi;

disp('   ');
disp('   Variance of nz: ')
var(5)
disp('   Variance of az = g*nz: ')
var(5)*g*g
