% Filename : exampl83.m
%
% Computation of analytical PSDs and time-simulation of aircraft 
% asymmetric response to atmospheric turbulence.

clc, clf, clear

disp('   Example 8.3');
disp('   ');
disp('   This example compares analytically obtained auto power');
disp('   spectral densities with experimentally obtained');
disp('   periodograms of the asymmetric motion variables.');
disp('   ');
disp('   This program produces Figures 8-20 and 8-21 of the lecture');
disp('   notes: Aircraft Responses to Atmospheric Turbulence.');
disp('   ');

% GET AIRCRAFT DYNAMICS
dynamics;

% NOTE (see also cit2a.m) 
%
% THE CESSNA CITATION CE-500 IS NOT STABLE IN SPIRAL MODE (FOR THE cit2a.m 
% FLIGHT CONDITION), HENCE THE FEEDBACK CONTROLLER FOR PHI IS USED AS IN : 
%
%   delta_a = K_phi*phi (K_phi for THIS flight condition)
%
% THEREFORE, CONTROLLED AIRCRAFT SYSTEM MATRICES WILL BE USED FOR RESULTS;
%
%      A = A2
   
% DEFINE C and D MATRICES
% C = [1 0 0 0 0 0 0 0 0 0   % beta
%      0 1 0 0 0 0 0 0 0 0   % phi
%      0 0 1 0 0 0 0 0 0 0   % pb/2V
%      0 0 0 1 0 0 0 0 0 0   % rb/2V
%      0 0 0 0 0 0 0 0 1 0]; % betag
% 
% D = [0 0 0 0 0
%      0 0 0 0 0
%      0 0 0 0 0
%      0 0 0 0 0
%      0 0 0 0 0];

% DEFINE FREQUENCY VECTOR
Nomega = 300;
w = logspace(-2,2,Nomega);

% COMPUTE ANALYTIC POWER SPECTRAL DENSITIES
% RESPONSE TO HORIZONTAL LATERAL TURBULENCE     ???????
% Full Model
temp = bode(A2,B,C(1,:),D(1,:),5,w); Sbeta  = temp.*temp;
temp = bode(A2,B,C(2,:),D(2,:),5,w); Sphi   = temp.*temp;
temp = bode(A2,B,C(3,:),D(3,:),5,w); Spp    = temp.*temp;
temp = bode(A2,B,C(4,:),D(4,:),5,w); Srr    = temp.*temp;
temp = bode(A2,B,V*(A2(1,:)+[0 0 0 2*V/b 0    0    0    0    0  0]),V*B(1,:),5,w);Say = temp.*temp;
%temp = bode(A2,B,C(5,:),D(5,:),5,w); Sbetag = temp.*temp;


temp = bode(Ar,Br,Cr(1,:),Dr(1,:),5,w); rSbeta  = temp.*temp;
%temp = bode(Ar,Br,Cr(2,:),Dr(2,:),5,w); rSphi   = temp.*temp;
%temp = bode(Ar,Br,Cr(3,:),Dr(3,:),5,w); rSpp    = temp.*temp;
temp = bode(Ar,Br,Cr(2,:),Dr(2,:),5,w); rSrr    = temp.*temp;
temp = bode(Ar,Br,V*(Ar(1,:)+[0 2*V/b 0    0    0    0    0  0]),V*Br(1,:),5,w);rSay = temp.*temp;
%temp = bode(Ar,Br,Cr(5,:),Dr(5,:),5,w); rSbetag = temp.*temp;

Sxx  = [Sbeta Sphi Spp Srr Say];
Sxxr = [rSbeta rSrr rSay];
% COMPUTE PSDS USING TIME DOMAIN DATA

% SET TIME AXIS
dt = 0.01; T = 60; t = [0:dt:T]; N = length(t);

% In this case responses to lateral gust vg are calculated (fifth input):
% no asymmetric vertical and longitudinal turbulence: u_g = w_g = 0.
v_g = randn(N,1)/sqrt(dt);    % sqrt(dt) because of lsim
w_g = randn(N,1)/sqrt(dt);
nn = zeros(N,1);
u  = [nn nn nn w_g v_g];

% COMPUTE SYSTEM RESPONSE
y     = lsim(A2,B,C,D,u,t);
yr    = lsim(Ar,Br,Cr,Dr,u,t);
beta  = y(:,1);
phi   = y(:,2);
pbtV  = y(:,3);
rbtV  = y(:,4);
%betag = y(:,5);
betar  = yr(:,1);
%phir   = yr(:,2);
%pbtVr  = yr(:,3);
rbtVr  = yr(:,2);

%Calculate lateral acceleration ay

ay = V*(A2(1,:)*y'+ B(1,:)*u' + (2*V/b)*y(:,4)');    % ay= V*(beta_dot + r)
ay = ay';
%Reduced

ayr = V*(Ar(1,:)*yr' +Br(1,:)*u' +(2*V/b)*yr(:,2)');
ayr = ayr';
% PLOT TIME RESPONSES
%clf
figure(1);
subplot(5,1,1); plot(t,beta); xlabel('time, s'); ylabel('beta (rad)');
subplot(5,1,2); plot(t,phi);  xlabel('time, s'); ylabel('phi (rad)');
subplot(5,1,3); plot(t,pbtV); xlabel('time, s'); ylabel('pb/2V (rad)');
subplot(5,1,4); plot(t,rbtV); xlabel('time, s'); ylabel('rb/2V (rad)');
subplot(5,1,5); plot(t,ay); xlabel('time, s'); ylabel('ay (m/s^2)');
suptitle('Time Responses (Full Aircraft Model)');



figure(2);
subplot(3,1,1); plot(t,betar); xlabel('time, s'); ylabel('beta (rad)');
%subplot(5,1,2); plot(t,phir);  xlabel('time, s'); ylabel('phi (rad)');
%subplot(5,1,3); plot(t,pbtVr); xlabel('time, s'); ylabel('pb/2V (rad)');
subplot(3,1,2); plot(t,rbtVr); xlabel('time, s'); ylabel('rb/2V (rad)');
subplot(3,1,3); plot(t,ayr); xlabel('time, s'); ylabel('ay (m/s^2)');
suptitle('Time Responses (Reduced Aircraft Model)')

%print -depsc2 -r1200 fig8_20a
%pause

% clf
% plot(t,betag); xlabel('time, s'); ylabel('betag');
%print -depsc2 -r1200 fig8_20b
%pause

% COMPUTE PERIODOGRAM AND ESTIMATE PSD
% CFT
BETA  = dt*fft(beta);
PHI   = dt*fft(phi);
P     = dt*fft(pbtV);
R     = dt*fft(rbtV);
AY    = dt*fft(ay);
%BETAg = dt*fft(betag);
BETAr  = dt*fft(betar);
%PHIr   = dt*fft(phir);
%Pr     = dt*fft(pbtVr);
Rr     = dt*fft(rbtVr);
AYr    = dt*fft(ayr);

% PSD ESTIMATE (Periodograms)
Pbeta  = (1/T)*( BETA.*conj(BETA));
Pphi   = (1/T)*(  PHI.*conj(PHI));
Pp     = (1/T)*(    P.*conj(P));
Pr     = (1/T)*(    R.*conj(R));
Pay    = (1/T)*(    AY.*conj(AY));
%Pbetag = (1/T)*(BETAg.*conj(BETAg));
Pbetar  = (1/T)*( BETAr.*conj(BETAr));
%Pphir   = (1/T)*(  PHIr.*conj(PHIr));
%Ppr     = (1/T)*(    Pr.*conj(Pr));
Prr     = (1/T)*(    Rr.*conj(Rr));
Payr    = (1/T)*(AYr.*conj(AYr));

Periodograms = [Pbeta Pphi Pp Pr Pay];
Periodograms_r = [Pbetar Prr Payr];
% DEFINE FREQUENCY VECTOR
fs = 1/dt;     % sample frequency
omega = 2*pi*fs*(0:(N/2)-1)/N;

% Smoothed PSD estimate

Pbeta_s  = zeros((length(Pbeta)-2),1);
Pphi_s   = zeros((length(Pphi)-2),1);
Pp_s     = zeros((length(Pp)-2),1);
Pr_s     = zeros((length(Pr)-2),1);
Pay_s  = zeros((length(Pay)-2),1);
Periodograms_ug = [Pbeta_s Pphi_s Pp_s Pr_s Pay_s];

for ii = 1:length(Periodograms(1,:))
    for jj = 2:length(Pbeta)-2
        Periodograms_ug(jj-1,ii) = 0.25*Periodograms(jj-1,ii)+0.5*Periodograms(jj,ii)+0.25*Periodograms(jj+1,ii);
    end
end

Pbeta_s_r  = zeros((length(Pbetar)-2),1);
Pr_s_r     = zeros((length(Prr)-2),1);
Pay_s_r  = zeros((length(Payr)-2),1);
Periodograms_ug_r = [Pbeta_s_r Pr_s_r Pay_s_r];

for ii = 1:length(Periodograms_r(1,:))
    for jj = 2:length(Pbetar)-2
        Periodograms_ug_r(jj-1,ii) = 0.25*Periodograms_r(jj-1,ii)+0.5*Periodograms_r(jj,ii)+0.25*Periodograms_r(jj+1,ii);
    end
end



%Plot PSD Estimates

figure(3);
subplot(2,2,1); loglog(w,Sxx(:,1),'g--',omega,Pbeta(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug(1:N/2-1,1),'-.','linewidth',1.5); 
axis(10.^[-2 2 -12 0]); xlabel('omega [rad/s]'); ylabel('Sbeta [rad^2/rad/s]')
subplot(2,2,2); loglog(w,Sxx(:,2),'g--',omega,Pphi(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug(1:N/2-1,2),'-.','linewidth',1.5);
axis(10.^[-2 2 -12 0]); xlabel('omega [rad/s]'); ylabel('Sphi [rad^2/rad/s]')
subplot(2,2,3); loglog(w,Sxx(:,3),'g--',omega,Pp(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug(1:N/2-1,3),'-.','linewidth',1.5);
axis(10.^[-2 2 -14 0]); xlabel('omega [rad/s]'); ylabel('Spp [rad^2/rad/s]')
subplot(2,2,4); loglog(w,Sxx(:,4),'g--',omega,Pr(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug(1:N/2-1,4),'-.','linewidth',1.5);
axis(10.^[-2 2 -14 0]); xlabel('omega [rad/s]'); ylabel('Srr [rad^2/rad/s]')
suptitle('State PSD Estimates (Full Aircraft Model)');

figure(4)
subplot(1,1,1); loglog(w,Sxx(:,5),'g--',omega,Pay(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug(1:N/2-1,5),'-.','linewidth',1.5);
axis(10.^[-2 2 -14 2]); xlabel('omega [rad/s]'); ylabel('Say [(m/s^2)^2/rad/s]');
suptitle('Lateral acceleration PSD estimate(Full Aircraft Model)');


figure(5);
subplot(2,2,1); loglog(w,Sxxr(:,1),'g--',omega,Pbetar(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug_r(1:N/2-1,1),'-.','linewidth',1.5); 
axis(10.^[-2 2 -12 0]); xlabel('omega [rad/s]'); ylabel('Sbeta [rad^2/rad/s]')
subplot(2,2,2); loglog(w,Sxxr(:,2),'g--',omega,Prr(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug_r(1:N/2-1,2),'-.','linewidth',1.5);
axis(10.^[-2 2 -14 0]); xlabel('omega [rad/s]'); ylabel('Srr [rad^2/rad/s]')
subplot(2,2,3); loglog(w,Sxxr(:,3),'g--',omega,Payr(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug_r(1:N/2-1,3),'-.','linewidth',1.5);
axis(10.^[-2 2 -14 2]); xlabel('omega [rad/s]'); ylabel('Say [(m/s^2)^2/rad/s]')
suptitle('State PSD Estimates (Full Aircraft Model)');

%Plot PSD comparison of the two aircraft models

figure(6);
subplot(2,2,1); loglog(w,Sxx(:,1),'r--',w,Sxxr(:,1)); 
axis(10.^[-2 2 -12 0]); xlabel('omega [rad/s]'); ylabel('Sbeta [rad^2/rad/s]')
subplot(2,2,2); loglog(w,Sxx(:,4),'r--',w,Sxxr(:,2));
axis(10.^[-2 2 -14 0]); xlabel('omega [rad/s]'); ylabel('Srr [rad^2/rad/s]')
subplot(2,2,3); loglog(w,Sxx(:,5),'r--',w,Sxxr(:,3));
axis(10.^[-2 2 -14 2]); xlabel('omega [rad/s]'); ylabel('Say [(m/s^2)^2/rad/s]')
suptitle('State PSD Estimates (Full (dashed red) vs Reduced (solid blue) Aircraft Models)');


%print -depsc2 -r1200 fig8_21a
% pause
% 
% clf
% loglog(w,Sxx(:,5),'--',omega,Pbetag(1:N/2)); 
% xlabel('omega [rad/s]'); ylabel('Sbetag [rad2]')
% %print -depsc2 -r1200 fig8_21b
% pause