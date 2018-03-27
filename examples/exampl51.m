% Exampl51  shows the response of a CT first order linear system to a white
%           noise input signal and examine the statistical properties
%           of the response for a large number of realizations.
%
%           Check the outcome with the result that can be
%           analytically obtained through Tables 3.5 and 3.6.
%
% Chapter 5 of lecture notes ae4-304.
%
% (c) MM 2004

clc
close all
clear all

disp('   Example 5.1                                                     ');
disp('                                                                   ');  
disp('   Simulation of a continuous time first order system response     ');
disp('   to a white noise input signal and examine the mean and          ');
disp('   standard deviation of the response for a large number           ');
disp('   of realizations                                                 ');
disp('                                                                   ');
disp('   The outcome can be verified with the result obtained            ');
disp('   analytically using Tables 3.5 and 3.6                           ');
disp('                                                                   ');
disp('   This program can produce Figure 5.2 of the lecture              ');
disp('   notes: Aircraft Responses to Atmospheric Turbulence.            ');
disp('                                                                   ');

% SYSTEM DYNAMICS
%
% Assume we have a CT system, a first order low-pass filter K/(1+s tau)
%
K   = 2.00;         % system gain     [-]
tau = 0.50;         % system time lag [s]

num = K; den = [tau 1];    % system dynamics in rational polynomial form

sys = tf(num,den);         % Matlab LTI system

% TIME DEFINITION
fs = 100;           % sample rate [Hz], 100 Hz is very common.
dt = 1/fs;          % sample time [s]
t  = 0:dt:60;       % 60 seconds of data
NT = length(t);     % the number of time samples

% WHITE NOISE DEFINITION
Wn = 2 ;                    % intensity of the white noise

w  = sqrt(Wn)*randn(NT,1);  % random, normally distributed white 
                            % noise with intensity Wn
%                              
% COMPUTE SYSTEM TIME RESPONSE

% REPEAT THE SYSTEM RESPONSE N TIMES AND LOOK AT VARIANCES
% OF INPUT AND OUTPUT SIGNALS
% DO NOT FORGET TO SKIP THE TRANSIENT PART OF THE RESPONSE!
N = 500;   % the number of realizations

Var_in  = zeros(N,1);
Var_out = zeros(N,1);

% skip the first part of the response, because of the transient
% in this case this only depends on the value of tau, but let's be
% safe and skip the first 33% anyway
stationary_part = (floor(NT/3)+1):1:NT;

for kk=1:N
   %
    % create white noise input
    w = sqrt(Wn)*randn(NT,1);
    % compute system output
    y = lsim(sys,w/sqrt(dt),t,'zoh'); % Correct: DIVISION of Wn by dt!
    y_e = lsim(sys,w,t,'zoh');        % Wrong: NO DIVISION of Wn by dt!
    % compute variance of input and output
    Var_in(kk,1)  = cov(w(stationary_part,1));
    Var_out(kk,1) = cov(y(stationary_part,1));
    Var_out_e(kk,1) = cov(y_e(stationary_part,1));
   %
end

% NOW, NOTE THAT WE ARE INTERESTED IN THE VARIANCE OF
% THE WHITE NOISE INPUT SIGNAL AND THE VARIANCE OF THE SYSTEM
% OUTPUT SIGNAL. FOR EACH REALIZATION, THESE ARE DIFFERENT:
% THEY CAN BE CONSIDERED STOCHASTIC VARIABLES THEMSELVES!!!
%
%   SO, if  var(w; realization i) = vw_i    i = 1, ... N
%       and var(y; realization i) = vy_i    i = 1, ... N
%
% WE CAN THEN COMPUTE THE MEAN AND THE VARIANCES OF
% THESE SVs:
%
%           mean(vw_i, for all i) and std(vw_i, for all i)
%           mean(vy_i, for all i) and std(vy_i, for all i)
%
% AND THE MEANS SHOULD IN PRINCIPLE (for a large number of
% realizations N) BE IDENTICAL TO THE ANALYTICAL VALUES
% REPRESENTING THE ENSEMBLE.
%
% REMEMBER: AVERAGING OVER THE REALIZATIONS MEANS THAT WE
%           ARE ESTIMATING THE ENSEMBLE AVERAGE. WHEN THE
%           NUMBER OF REALIZATIONS INCREASES (N -> infinity)
%           THE ESTIMATIONS SHOULD CONVERGE TO THE ENSEMBLE
%           AVERAGE. IT CAN BE SHOWN THAT THE ESTIMATOR
%           FOR THE MEAN AND THE ESTIMATOR FOR THE VARIANCE
%           ARE UNBIASED AND ASYMPTOTICALLY RIGHT
%
%           (check in Matlab:  help mean
%                              help var  )
%
mean_var_in    = mean(Var_in);
mean_var_out   = mean(Var_out);
mean_var_out_e = mean(Var_out_e);
var_var_in     = var(Var_in); 
var_var_out    = var(Var_out);
var_var_out_e  = var(Var_out_e);

%
% THE ANALYTICAL VALUE OF THE SYSTEM OUTPUT SIGNAL VARIANCE
% CAN BE OBTAINED WITH TABLE 3.5
%
% For this system it is equal to:
%       K^2
%    W*-----       with: W   : the CT white noise intensity
%      2*tau             K   : the system gain
%                        tau : the system lag time constant
%
var_out_analytic = Wn*K*K/(2*tau);

% THE ANALYTICAL VALUE OF THE WHITE NOISE INPUT VARIANCE
% JUST EQUALS THE INTENSITY
var_in_analytic = Wn;

% PLOT THE RESULTS
figure(1)
% first the white noise input
subplot(3,1,1)
% show the values for each realization
plot([1:1:N],Var_in)
hold on
% show the mean +- var of the average over 
% all realizations in green
plot([1 N],mean_var_in*[1 1],'g');             
plot([1 N],mean_var_in*[1 1]+var_var_in,'g--'); 
% show the analytical value (the ensemble average) in red
plot([1 N],var_in_analytic*[1 1],'r');         
% plot the other var line here, otherwise it is 
% redundant in the legend
plot([1 N],mean_var_in*[1 1]-var_var_in,'g--');
hold off
title('variance of CT white noise input signal')
axis('tight')
ylabel('var(w)')
xlabel('realization')
legend('estimated','mean','variance','analytic')

% the system output without division of Wn by dt
subplot(3,1,2)
% show the values for each realization
plot([1:1:N],Var_out_e)
hold on
% show the mean +- var of the average over 
% all realizations in green
plot([1 N],mean_var_out_e*[1 1],'g');          
plot([1 N],mean_var_out_e*[1 1]+var_var_out_e,'g--'); 
% show the analytical value (the ensemble average) in red
plot([1 N],var_out_analytic*[1 1],'r');            
% plot the other var line here, otherwise it is 
% redundant in the legend
plot([1 N],mean_var_out_e*[1 1]-var_var_out_e,'g--');
hold off
title('variance of system output signal no dt')
ylabel('var(y)')
xlabel('realization')
legend('estimated','mean','var','analytic')

% the system output with Wn/dt
subplot(3,1,3)
% show the values for each realization
plot([1:1:N],Var_out)
hold on
% show the mean +- var of the average over 
% all realizations in green
plot([1 N],mean_var_out*[1 1],'g');              
plot([1 N],mean_var_out*[1 1]+var_var_out,'g--');
% show the analytical value (the ensemble average) in red
plot([1 N],var_out_analytic*[1 1],'r');          
% plot the other var line here, otherwise it is 
% redundant in the legend
plot([1 N],mean_var_out*[1 1]-var_var_out,'g--');
hold off
title('variance of system output signal')
axis('tight')
ylabel('var(y)')
xlabel('realization')
legend('estimated','mean','variance','analytic')
% EOF