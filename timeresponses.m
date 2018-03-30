% Filename : examp81c.m
%
% Simulation of aircraft asymmetric response to atmospheric turbulence.

clc, clf, clear, close all

% GET SYSTEM DYNAMICS
dynamics

% TIME AXIS AND INPUT VECTOR DEFINITION
dt = 0.005; T  = 60; t = [0:dt:T]; N = length(t);
nn = zeros(1,N);
% TURBULENCE INPUTS
u_g = randn(1,N)/sqrt(dt);    % sqrt(dt) because of lsim characteristics
v_g = randn(1,N)/sqrt(dt);
w_g = randn(1,N)/sqrt(dt);

% INPUT VECTORS
u1 = [nn' nn' nn' nn'  nn'];
%u1 = [nn' nn' nn' nn' nn'];
u2 = [nn' nn' nn'  nn' v_g'];
%u2 = [nn' nn' nn' nn' nn'];
u3 = [nn' nn' nn'  w_g'  nn'];
%u3 = [nn' nn' nn' nn' nn'];

% RESPONSE to u_g
y1 = lsim(A2,B,C,D,u1,t);
% RESPONSE to v_g
y2 = lsim(A2,B,C,D,u2,t);
% RESPONSE to w_g
y3 = lsim(A2,B,C,D,u3,t);
% RESPONSE to all together (linear system!)
yt = y1+y2+y3;

%%Alternate method of calculating ay. To calidate the calculated values
% C=[C; V*(A2(1,:)+[0 0 0 2*V/b 0    0    0    0    0  0])];
% D=[D; V*B(1,:)];
% 
% y1n = lsim(A2,B,C,D,u1,t);
% % RESPONSE to v_g
% y2n = lsim(A2,B,C,D,u2,t);
% % RESPONSE to w_g
% y3n = lsim(A2,B,C,D,u3,t);
% 
% ytn = y1n + y2n + y3n;

%%End alternate method

% For reduced model

% RESPONSE to u_g
y1r = lsim(Ar,Br,Cr,Dr,u1,t);
% RESPONSE to v_g
y2r = lsim(Ar,Br,Cr,Dr,u2,t);
% RESPONSE to w_g
y3r = lsim(Ar,Br,Cr,Dr,u3,t);
% RESPONSE to all together (linear system!)
ytr = y1r+y2r+y3r;

%Calculate lateral acceleration ay

ay = V*(A2(1,:)*yt'+ B(1,:)*u1' + B(1,:)*u2' + B(1,:)*u3' + (2*V/b)*yt(:,4)');    % ay= V*(beta_dot + r)

%Reduced

ayr = V*(Ar(1,:)*ytr' +Br(1,:)*u1' + Br(1,:)*u2' + Br(1,:)*u3'+ (2*V/b)*ytr(:,2)');

%Response plots

% beta_axis = [0 60 -0.07  0.07];
% phi_axis  = [0 60 -0.15  0.15];
% pb_axis   = [0 60 -1e-2  1e-2];
% rb_axis   = [0 60 -1e-2  1e-2];

% clf;
% figure(1);
% plot(t,yt(:,1));
% xlabel('time, s'); ylabel('beta');
% title('Sideslip angle (full model)');
% 
% figure(7);
% plot(t,yt(:,1));
% xlabel('time, s'); ylabel('beta');
% title('Sideslip angle (reduced model)');
% 
% figure(2);
% plot(t,yt(:,2));
% xlabel('time, s'); ylabel('phi');
% title('Roll angle (full model)');
% 
% figure(8);
% plot(t,yt(:,2));
% xlabel('time, s'); ylabel('beta');
% title('Sideslip angle (reduced model)');
% 
% figure(3);
% plot(t,yt(:,3));
% xlabel('time, s'); ylabel('pb/2V');
% title('Pitch rate (pb/2V) (full model)');
% 
% figure(9);
% plot(t,yt(:,3));
% xlabel('time, s'); ylabel('beta');
% title('Sideslip angle (reduced model)');
% 
% figure(4);
% plot(t,yt(:,4));
% xlabel('time, s'); ylabel('rb/2V');
% title('Yaw rate (rb/2V) (full model)');
% figure(10);
% plot(t,yt(:,1));
% xlabel('time, s'); ylabel('beta');
% title('Sideslip angle (reduced model)');
% 
% figure(5);
% plot(t,ay);
% xlabel('time, s'); ylabel('ay');
% title('Lateral acceleration (full model)');
% 
% figure(6);
% plot(t,ayr);
% xlabel('time, s'); ylabel('ay');
% title('Lateral acceleration (reduced model)');



% RESPONSE TO u_g
% clf;
% subplot(2,1,1); plot(t,y1(:,1)); axis(beta_axis); 
% ah=gca; set(ah,'Fontsize',14);
% xlabel('time, s'); ylabel('beta');
% subplot(2,1,2); plot(t,y1(:,2)); axis(phi_axis); 
% ah=gca; set(ah,'Fontsize',14);
% xlabel('time, s'); ylabel('phi');
% %print -depsc2 -r1200 fig8_16a1
% pause
% 
% subplot(2,1,1); plot(t,y1(:,3)); axis(pb_axis); 
% xlabel('time, s'); ylabel('pb/2V');
% subplot(2,1,2); plot(t,y1(:,4)); axis(rb_axis);
% xlabel('time, s'); ylabel('rb/2V');
% %print -depsc2 -r1200 fig8_16b1
% pause
% 
% % RESPONSE TO v_g
% subplot(2,1,1); plot(t,y2(:,1)); axis(beta_axis); 
% xlabel('time, s'); ylabel('beta');
% subplot(2,1,2); plot(t,y2(:,2)); axis(phi_axis);
% xlabel('time, s'); ylabel('phi');
% %print -depsc2 -r1200 fig8_16a2
% pause
% 
% subplot(2,1,1); plot(t,y2(:,3)); axis(pb_axis); 
% xlabel('time, s'); ylabel('pb/2V');
% subplot(2,1,2); plot(t,y2(:,4)); axis(rb_axis);
% xlabel('time, s'); ylabel('rb/2V');
% %print -depsc2 -r1200 fig8_16b2
% pause
% 
% % RESPONSE TO w_g
% subplot(2,1,1); plot(t,y3(:,1)); axis(beta_axis); 
% xlabel('time, s'); ylabel('beta');
% subplot(2,1,2); plot(t,y3(:,2)); axis(phi_axis);
% xlabel('time, s'); ylabel('phi');
% %print -depsc2 -r1200 fig8_17a3
% pause
% 
% subplot(2,1,1); plot(t,y3(:,3)); axis(pb_axis); 
% xlabel('time, s'); ylabel('pb/2V');
% subplot(2,1,2); plot(t,y3(:,4)); axis(rb_axis);
% xlabel('time, s'); ylabel('rb/2V');
% %print -depsc2 -r1200 fig8_17b3
% pause

% RESPONSE TO all together
figure(1);
subplot(5,1,1); plot(t,yt(:,1)); axis([0 60 -0.2 0.2]);
grid on;
xlabel('time [s]'); ylabel('beta [rad]');
subplot(5,1,2); plot(t,yt(:,2));axis([0 60 -0.2 0.2]);
grid on;
xlabel('time [s]'); ylabel('phi [rad]');
%print -depsc2 -r1200 fig8_17at
subplot(5,1,3); plot(t,yt(:,3));axis([0 60 -0.02 0.02]);
grid on;
xlabel('time [s]'); ylabel('pb/2V [rad]');
subplot(5,1,4); plot(t,yt(:,4));axis([0 60 -0.02 0.02]);
grid on;
xlabel('time [s]'); ylabel('rb/2V [rad]');
subplot(5,1,5); plot(t,ay); axis([0 60 -2.2 2.2]);
grid on;
xlabel('time [s]'); ylabel('ay [m/s^2]');
suptitle('Time responses of aircraft states : Full aircraft model')
figure(2);
subplot(3,1,1); plot(t,ytr(:,1)); axis([0 60 -0.2 0.2]);
grid on;
xlabel('time [s]'); ylabel('beta [rad]');
subplot(3,1,2); plot(t,ytr(:,2)); axis([0 60 -0.02 0.02]);
grid on;
xlabel('time [s]'); ylabel('rb/2V [rad]');
subplot(3,1,3); plot(t,ayr); axis([0 60 -1 1]);
grid on;
xlabel('time [s]'); ylabel('ay [m/s^2]');
suptitle('Time responses of aircraft states : Reduced aircraft model')

%%This is to check the calculated ay values
% figure(3)
% plot(t,ay,'o');
% hold on;
% plot(t,ytn(:,11));
% hold off;
%print -depsc2 -r1200 fig8_17bt
% pause