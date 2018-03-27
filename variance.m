%Variance Calculation
clc, clf, clear

dynamics;
%omega space definition
Nomega = 300;
w = logspace(-2,2,Nomega);

%Full model
temp = bode(A2,B,C(1,:),D(1,:),5,w); Sbeta  = temp.*temp;
temp = bode(A2,B,C(2,:),D(2,:),5,w); Sphi   = temp.*temp;
temp = bode(A2,B,C(3,:),D(3,:),5,w); Spp    = temp.*temp;
temp = bode(A2,B,C(4,:),D(4,:),5,w); Srr    = temp.*temp;
temp = bode(A2,B,V*(A2(1,:)+[0 0 0 2*V/b 0    0    0    0    0  0]),V*B(1,:),5,w);Say = temp.*temp;

%Reduced model
temp = bode(Ar,Br,Cr(1,:),Dr(1,:),5,w); rSbeta  = temp.*temp;
temp = bode(Ar,Br,Cr(2,:),Dr(2,:),5,w); rSrr    = temp.*temp;
temp = bode(Ar,Br,V*(Ar(1,:)+[0 2*V/b 0    0    0    0    0  0]),V*Br(1,:),5,w);rSay = temp.*temp;


Sxx  = [Sbeta Sphi Spp Srr Say];
Sxxr = [rSbeta rSrr rSay];

% SET TIME AXIS
dt = 0.01; T = 60; t = [0:dt:T]; N = length(t);

v_g = randn(N,1)/sqrt(d t);    % sqrt(dt) because of lsim
w_g = randn(N,1)/sqrt(dt);
nn = zeros(N,1);
u  = [nn nn nn w_g v_g];

% COMPUTE SYSTEM RESPONSE
y     = lsim(A2,B,C,D,u,t);
yr    = lsim(Ar,Br,Cr,Dr,u,t);

%Calculate lateral acceleration ay

ay = V*(A2(1,:)*y'+ B(1,:)*u' + (2*V/b)*y(:,4)');    % ay= V*(beta_dot + r)
ay = ay';
%Reduced

ayr = V*(Ar(1,:)*yr' +Br(1,:)*u' +(2*V/b)*yr(:,2)');
ayr = ayr';


%Variance calculation through integration of analytic PSD
var_PSD = zeros(1,5);
for i=1:1:5
    for j=1:1:Nomega-1
        var_PSD(i)=var_PSD(i)+(w(j+1)-w(j))*Sxx(j,i);
    end
end
var_PSD = var_PSD/pi;

%For reduced model

var_PSD_r = zeros(1,3);
for i=1:1:3
    for j=1:1:Nomega-1
        var_PSD_r(i)=var_PSD_r(i)+(w(j+1)-w(j))*Sxxr(j,i);
    end
end
var_PSD_r = var_PSD_r/pi;

%Variance calculation through Lyapnov's equation

Bl=B(:,4:5);
L   = lyap(A2,Bl*1*Bl');
L   = L(1:4,1:4);
var_L = diag(L)';

%For reduced model

Bl=Br(:,4:5);
Lr   = lyap(Ar,Bl*1*Bl');
Lr   = Lr(1:2,1:2);
var_L_r = diag(Lr)';

%Variance calculation through matlab var function

varBeta  = var(y(:,1));
varPhi  = var(y(:,2));
varp  = var(y(:,3));
varr  = var(y(:,4));
varay  = var(ay);
    
var_3 = [varBeta varPhi varp varr varay];

%reduced model

varBeta_r  = var(yr(:,1));
varr_r  = var(yr(:,2));
varay_r  = var(ayr);
    
var_3_r = [varBeta_r varr_r varay_r];





