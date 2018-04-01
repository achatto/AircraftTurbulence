% Filename : timeresponses.m

% Abhishek Chatterjee (4743075)
% Assignment : AE4304P Stochastic Aerospace Systems Practical
% Delft University of Technology
% Simulation of aircraft asymmetric response to atmospheric turbulence.

% Time domain simulation of aircraft responses

clc, clf, clear, close all

% GET SYSTEM DYNAMICS
dynamics; %See dynamics.m
close all;

% TIME AXIS AND INPUT VECTOR DEFINITION
dt = 0.005; T  = 60; t = [0:dt:T]; N = length(t);
nn = zeros(1,N);
% TURBULENCE INPUTS
u_g = randn(1,N)/sqrt(dt);    % sqrt(dt) because of lsim characteristics
v_g = randn(1,N)/sqrt(dt);
w_g = randn(1,N)/sqrt(dt);

% INPUT VECTORS
u1 = [nn' nn' nn' nn'  nn'];
u2 = [nn' nn' nn'  nn' v_g'];
u3 = [nn' nn' nn'  w_g'  nn'];

% RESPONSE to u1
y1 = lsim(A2,B,C,D,u1,t);
% RESPONSE to u2
y2 = lsim(A2,B,C,D,u2,t);
% RESPONSE to u3
y3 = lsim(A2,B,C,D,u3,t);
% RESPONSE to all together (linear system!)
yt = y1+y2+y3;

% For reduced model

% RESPONSE to u1
y1r = lsim(Ar,Br,Cr,Dr,u1,t);
% RESPONSE to u2
y2r = lsim(Ar,Br,Cr,Dr,u2,t);
% RESPONSE to u3
y3r = lsim(Ar,Br,Cr,Dr,u3,t);
% RESPONSE to all together (linear system!)
ytr = y1r+y2r+y3r;

%Calculate lateral acceleration ay

ay = V*(A2(1,:)*yt'+ B(1,:)*u1' + B(1,:)*u2' + B(1,:)*u3' + (2*V/b)*yt(:,4)');    % ay= V*(beta_dot + r)

%Reduced Model

ayr = V*(Ar(1,:)*ytr' +Br(1,:)*u1' + Br(1,:)*u2' + Br(1,:)*u3'+ (2*V/b)*ytr(:,2)');


% Time response plots
figure(1);
subplot(5,1,1); plot(t,yt(:,1)); axis([0 60 -0.25 0.25]);
grid on;
xlabel('time [s]'); ylabel('beta [rad]');set(gca,'fontsize',15);
subplot(5,1,2); plot(t,yt(:,2));axis([0 60 -0.25 0.25]);
grid on;
xlabel('time [s]'); ylabel('phi [rad]');set(gca,'fontsize',15);
%print -depsc2 -r1200 fig8_17at
subplot(5,1,3); plot(t,yt(:,3));axis([0 60 -0.025 0.025]);
grid on;
xlabel('time [s]'); ylabel('pb/2V [rad]');set(gca,'fontsize',15);
subplot(5,1,4); plot(t,yt(:,4));axis([0 60 -0.025 0.025]);
grid on;
xlabel('time [s]'); ylabel('rb/2V [rad]');set(gca,'fontsize',15);
subplot(5,1,5); plot(t,ay); axis([0 60 -2.25 2.25]);
grid on;
xlabel('time [s]'); ylabel('ay [m/s^2]');set(gca,'fontsize',15);
suptitle('Time responses of aircraft states : Full aircraft model');set(gca,'fontsize',15);
figure(2);
subplot(3,1,1); plot(t,ytr(:,1)); axis([0 60 -0.25 0.25]);
grid on;
xlabel('time [s]'); ylabel('beta [rad]');set(gca,'fontsize',15);
subplot(3,1,2); plot(t,ytr(:,2)); axis([0 60 -0.025 0.025]);
grid on;
xlabel('time [s]'); ylabel('rb/2V [rad]');set(gca,'fontsize',15);
subplot(3,1,3); plot(t,ayr); axis([0 60 -1 1]);
grid on;
xlabel('time [s]'); ylabel('ay [m/s^2]');set(gca,'fontsize',15);
suptitle('Time responses of aircraft states : Reduced aircraft model');set(gca,'fontsize',15);

%%This is to check the calculated ay values
% figure(3)
% plot(t,ay,'o');
% hold on;
% plot(t,ytn(:,11));
% hold off;
%print -depsc2 -r1200 fig8_17bt
% pause