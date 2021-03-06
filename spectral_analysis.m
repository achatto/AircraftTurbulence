% Filename : spectral_analysis.m
%
% Abhishek Chatterjee (4743075)
% Assignment : AE4304P Stochastic Aerospace Systems Practical
% Delft University of Technology
% Simulation of Aircraft Asymmetric Response to atmospheric Turbulence

% Calcultion of the analytical and experimental Power Spectral Density
% functions of the aircraft states. 

clc, clf, clear


% GET AIRCRAFT DYNAMICS
dynamics; %See dynamics.m
close all;

% DEFINE FREQUENCY VECTOR
Nomega = 300;
w = logspace(-2,2,Nomega);

% COMPUTE ANALYTIC POWER SPECTRAL DENSITIES

% Full Model
temp = bode(A2,B,C(1,:),D(1,:),4,w);
temp = temp + bode(A2,B,C(1,:),D(1,:),5,w); Sbeta  = temp.*temp;
temp = bode(A2,B,C(2,:),D(2,:),4,w);
temp = temp + bode(A2,B,C(2,:),D(2,:),5,w); Sphi   = temp.*temp;
temp = bode(A2,B,C(3,:),D(3,:),4,w);
temp = temp + bode(A2,B,C(3,:),D(3,:),5,w); Spp    = temp.*temp;
temp = bode(A2,B,C(4,:),D(4,:),4,w);
temp = temp + bode(A2,B,C(4,:),D(4,:),5,w); Srr    = temp.*temp;
temp = bode(A2,B,V*(A2(1,:)+[0 0 0 2*V/b 0    0    0    0    0  0]),V*B(1,:),4,w);
temp = temp + bode(A2,B,V*(A2(1,:)+[0 0 0 2*V/b 0    0    0    0    0  0]),V*B(1,:),5,w);Say = temp.*temp;

%Reduced Model
temp = bode(Ar,Br,Cr(1,:),Dr(1,:),4,w);
temp = temp + bode(Ar,Br,Cr(1,:),Dr(1,:),5,w); rSbeta  = temp.*temp;
temp = bode(Ar,Br,Cr(2,:),Dr(2,:),4,w);
temp = temp + bode(Ar,Br,Cr(2,:),Dr(2,:),5,w); rSrr    = temp.*temp;
temp = bode(Ar,Br,V*(Ar(1,:)+[0 2*V/b 0    0    0    0    0  0]),V*Br(1,:),4,w);
temp = temp + bode(Ar,Br,V*(Ar(1,:)+[0 2*V/b 0    0    0    0    0  0]),V*Br(1,:),5,w);rSay = temp.*temp;

Sxx  = [Sbeta Sphi Spp Srr Say];
Sxxr = [rSbeta rSrr rSay];

% COMPUTE EXPERIMENTAL PSDs

% COMPUTE PSDS USING TIME DOMAIN DATA

% SET TIME AXIS
dt = 0.01; T = 60; t = [0:dt:T]; N = length(t);

v_g = randn(N,1)/sqrt(dt);    % sqrt(dt) because of lsim
w_g = randn(N,1)/sqrt(dt);
nn = zeros(N,1);
u = [nn nn nn w_g v_g];

% COMPUTE SYSTEM RESPONSE
%Full Model
y     = lsim(A2,B,C,D,u,t);
beta  = y(:,1);
phi   = y(:,2);
pbtV  = y(:,3);
rbtV  = y(:,4);

%Reduced Model
yr    = lsim(Ar,Br,Cr,Dr,u,t);
betar  = yr(:,1);
rbtVr  = yr(:,2);

%Calculate lateral acceleration ay

%Full Model
ay = V*(A2(1,:)*y'+ B(1,:)*u' + (2*V/b)*y(:,4)');    % ay= V*(beta_dot + r)
ay = ay';
%Reduced Model
ayr = V*(Ar(1,:)*yr' +Br(1,:)*u' +(2*V/b)*yr(:,2)');
ayr = ayr';

% % PLOT TIME RESPONSES
% clf
% figure(1);
% subplot(5,1,1); plot(t,beta); xlabel('time, s'); ylabel('beta (rad)');
% subplot(5,1,2); plot(t,phi);  xlabel('time, s'); ylabel('phi (rad)');
% subplot(5,1,3); plot(t,pbtV); xlabel('time, s'); ylabel('pb/2V (rad)');
% subplot(5,1,4); plot(t,rbtV); xlabel('time, s'); ylabel('rb/2V (rad)');
% subplot(5,1,5); plot(t,ay); xlabel('time, s'); ylabel('ay (m/s^2)');
% suptitle('Time Responses (Full Aircraft Model)');
% 
% 
% 
% figure(2);
% subplot(3,1,1); plot(t,betar); xlabel('time, s'); ylabel('beta (rad)');
% %subplot(5,1,2); plot(t,phir);  xlabel('time, s'); ylabel('phi (rad)');
% %subplot(5,1,3); plot(t,pbtVr); xlabel('time, s'); ylabel('pb/2V (rad)');
% subplot(3,1,2); plot(t,rbtVr); xlabel('time, s'); ylabel('rb/2V (rad)');
% subplot(3,1,3); plot(t,ayr); xlabel('time, s'); ylabel('ay (m/s^2)');
% suptitle('Time Responses (Reduced Aircraft Model)')


% COMPUTE PERIODOGRAM AND ESTIMATE PSD
% CFT
%Full Model
BETA  = dt*fft(beta);
PHI   = dt*fft(phi);
P     = dt*fft(pbtV);
R     = dt*fft(rbtV);
AY    = dt*fft(ay);

%Reduced Model
BETAr  = dt*fft(betar);
Rr     = dt*fft(rbtVr);
AYr    = dt*fft(ayr);

% PSD ESTIMATE (Periodograms)

%Full Model
Pbeta  = (1/T)*( BETA.*conj(BETA));
Pphi   = (1/T)*(  PHI.*conj(PHI));
Pp     = (1/T)*(    P.*conj(P));
Pr     = (1/T)*(    R.*conj(R));
Pay    = (1/T)*(    AY.*conj(AY));

%Reduced Model
Pbetar  = (1/T)*( BETAr.*conj(BETAr));
Prr     = (1/T)*(    Rr.*conj(Rr));
Payr    = (1/T)*(AYr.*conj(AYr));

Periodograms = [Pbeta Pphi Pp Pr Pay];
Periodograms_r = [Pbetar Prr Payr];

% DEFINE FREQUENCY VECTOR
fs = 1/dt;     % sample frequency
omega = 2*pi*fs*(0:(N/2)-1)/N;

% SMOOTHED PERIODOGRAMS

%Full Model
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

%Reduced Model
Pbeta_s_r  = zeros((length(Pbetar)-2),1);
Pr_s_r     = zeros((length(Prr)-2),1);
Pay_s_r  = zeros((length(Payr)-2),1);
Periodograms_ug_r = [Pbeta_s_r Pr_s_r Pay_s_r];

for ii = 1:length(Periodograms_r(1,:))
    for jj = 2:length(Pbetar)-2
        Periodograms_ug_r(jj-1,ii) = 0.25*Periodograms_r(jj-1,ii)+0.5*Periodograms_r(jj,ii)+0.25*Periodograms_r(jj+1,ii);
    end
end



%PLOT PSD ESTIMATES

figure(3); %Full AIrcraft Model
subplot(2,2,1); loglog(w,Sxx(:,1),'g--',omega,Pbeta(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug(1:N/2-1,1),'-.','linewidth',1.5); 
axis(10.^[-2 2 -12 0]); xlabel('omega [rad/s]'); ylabel('PSD sideslip Sbeta [rad^2/rad/s]');set(gca,'fontsize',15);
subplot(2,2,2); loglog(w,Sxx(:,2),'g--',omega,Pphi(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug(1:N/2-1,2),'-.','linewidth',1.5);
axis(10.^[-2 2 -12 0]); xlabel('omega [rad/s]'); ylabel('PSD roll Sphi [rad^2/rad/s]');set(gca,'fontsize',15);
subplot(2,2,3); loglog(w,Sxx(:,3),'g--',omega,Pp(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug(1:N/2-1,3),'-.','linewidth',1.5);
axis(10.^[-2 2 -14 0]); xlabel('omega [rad/s]'); ylabel('PSD roll ratw Spp [rad^2/rad/s]');set(gca,'fontsize',15);
subplot(2,2,4); loglog(w,Sxx(:,4),'g--',omega,Pr(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug(1:N/2-1,4),'-.','linewidth',1.5);
axis(10.^[-2 2 -14 0]); xlabel('omega [rad/s]'); ylabel('PSD yaw rate Srr [rad^2/rad/s]');set(gca,'fontsize',15);
legend('Analytical PSD', 'Periodogram','Smoothed Periodogram');set(gca,'fontsize',10);
suptitle('State PSD Estimates (Full Aircraft Model)');set(gca,'fontsize',15);

figure(4); %Lateral Acceleration Full Aircraft Model
subplot(1,1,1); loglog(w,Sxx(:,5),'g--',omega,Pay(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug(1:N/2-1,5),'-.','linewidth',1.5);
axis(10.^[-2 2 -14 2]); xlabel('omega [rad/s]'); ylabel('PSD lateral acceleration Say [(m/s^2)^2/rad/s]');set(gca,'fontsize',15);
legend('Analytical PSD', 'Periodogram','Smoothed Periodogram');set(gca,'fontsize',10);
suptitle('Lateral acceleration PSD estimate(Full Aircraft Model)');set(gca,'fontsize',15);


figure(5); %Reduced Aircraft Model
subplot(2,1,1); loglog(w,Sxxr(:,1),'g--',omega,Pbetar(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug_r(1:N/2-1,1),'-.','linewidth',1.5); 
axis(10.^[-2 2 -12 0]); xlabel('omega [rad/s]'); ylabel('PSD sideslip Sbeta\n [rad^2/rad/s]');set(gca,'fontsize',15);
subplot(2,1,2); loglog(w,Sxxr(:,2),'g--',omega,Prr(1:N/2),'-',omega(1:length(omega)-1),Periodograms_ug_r(1:N/2-1,2),'-.','linewidth',1.5);
axis(10.^[-2 2 -14 0]); xlabel('omega [rad/s]'); ylabel('PSD yaw rate Srr\n [rad^2/rad/s]');set(gca,'fontsize',15);
legend('Analytical PSD', 'Periodogram','Smoothed Periodogram');set(gca,'fontsize',10);
suptitle('State PSD Estimates (Reduced Aircraft Model)');set(gca,'fontsize',15);


%Plot PSD comparison of the two aircraft models

figure(6);
subplot(2,1,1); loglog(w,Sxx(:,1),'r--',w,Sxxr(:,1)); 
axis(10.^[-2 2 -12 0]); xlabel('omega [rad/s]'); ylabel('PSD sideslip Sbeta [rad^2/rad/s]');set(gca,'fontsize',15);
subplot(2,1,2); loglog(w,Sxx(:,4),'r--',w,Sxxr(:,2));
axis(10.^[-2 2 -14 0]); xlabel('omega [rad/s]'); ylabel('PSD yaw rate Srr [rad^2/rad/s]');set(gca,'fontsize',15);
legend('Full Aircraft Model','Reduced Aircraft Model');set(gca,'fontsize',10);
suptitle('Aircraft state PSD Estimates (Full vs Reduced aircraft models)');set(gca,'fontsize',15);