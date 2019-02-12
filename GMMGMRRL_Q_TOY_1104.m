%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% This file contains a sample implementation of PoWER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% The algorithm is implemented with a number of simplifications:
% - the variance of the exploration is constant over trials
% - the exploration is constant during the trial
%   (as the motor primitives employ basis functions that are localized in time 
%   which are only active for a short period of time,
%   time-varying exploration does not have large effetcs)
% - always only one basis function is active
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
clc;
%% Definition of the number of components used in GMM.
nbStates = 9;
%% Load a dataset consisting of 3 demonstrations for motor commands.
load('motor.mat');
Data=motor;
nbVar = size(Data,1);
nbSamples=1257;
variance = 0.1.*ones(nbVar,nbStates);
%% Training of GMM by EM algorithm, initialized by K-means clustering.
[Priors, Mu, Sigma] = EM_init_kmeans(Data, nbStates);
[Priors, Mu, Sigma] = EM(Data, Priors, Mu, Sigma);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Number of iterations
n_iter = 50;

Return = zeros(1,n_iter+1);
Q = zeros(200,n_iter+1);
s_Return = zeros(n_iter+1,2);
param = zeros(4*nbStates,n_iter+1);

%Initialize parameters for the EM-RL
param(:,1) = reshape(Mu,4*nbStates,1); 
current_param = reshape(param(:,1),4,nbStates);

%% Use of GMR to retrieve a generalized version of the data and associated
%% constraints. A sequence of temporal values is used as input, and the 
%% expected distribution is retrieved. 
expData(1,:) = linspace(min(Data(1,:)), max(Data(1,:)), 200);
[expData(2:nbVar,:), expSigma] = GMR(Priors, current_param, Sigma, expData(1,:), [1], [2:nbVar]);
generalized_rm=expData(2,:);
generalized_lm=expData(3,:);
generalized_zb=expData(4,:);
%% Execute the motor commands and receive end-effector trajectory from EM tracker

%Calculate the reward function
for i=1:200
    Q(1:(end+1-i),1) = Q(1:(end+1-i),1) + exp(-abs(sqrt((x(i)-x0)^2+(y(i)-y0)^2)-30)+abs(z(i)-z0));
    Q(:,1) = Q(:,1)./200;
    plot3(x(i),y(i),z(i),'ro','MarkerSize',6);
    hold on;
    pause(0.025);
end
%Normalize the Q values

%% Do the Iterations of the Reinforcement Learning Algorithm
for iter=1:n_iter
    if (mod(iter,10)==0)
        disp(['Number of Iteration of the Reinforcement Learning=', num2str(iter)]);
    end
    Return(iter) = Q(1,iter);

	%This lookup table will be used for the importance sampling
    s_Return(1,:) = [Return(iter) iter];
    s_Return = sortrows(s_Return);
    
    %Update the policy parameters
    param_nom = zeros(nbVar*nbStates,1);
    param_dnom = zeros(nbVar*nbStates,1);

	%Calculate the expectations (the normalization is taken care of by the division)
    %As importance sampling we take the 10 best rollouts
    for i=1:min(iter,10)
        %get the rollout number for the 10 best rollouts
        j = s_Return(end+1-i,2);
		%calculate the exploration with respect to the current parameters
        %if you have time-varying exploration use 
        %temp_explore = (reshape(param(:,:,j),length(rt),n_rfs)-ones(length(rt),1)*current_param')';
        %instead
        temp_explore = (ones(200,nbStates)*(param(:,j)-reshape(current_param,nbVar*nbStates,1))')';
		%repeat the Q values
        temp_Q = (Q(:,j)*ones(1,nbVar*nbStates))';
        param_nom = param_nom + sum(temp_explore.*temp_Q,2);
        param_dnom = param_dnom + sum(temp_Q,2);
    end
    
    %Update the parameters
    param(:,iter+1) = current_param + param_nom./(param_dnom+1.e-10);
    %Set the new mean of the parameters
    current_param = param(:,:,iter+1);
    %In the last rollout we want to get the return without exploration
    if iter~=n_iter
        param(:,iter+1) = param(:,iter+1) + variance.^.5.*randn(nbVar,nbStates);
        current_param = reshape(param(:,1),4,nbStates);
    end
    %% Use of GMR to generalize motor commands for the TSM manipulator
    expData(1,:) = linspace(min(Data(1,:)), max(Data(1,:)), 200);
    [expData(2:nbVar,:), expSigma] = GMR(Priors, current_param, Sigma, expData(1,:), [1], [2:nbVar]);
    generalized_rm=expData(2,:);
    generalized_lm=expData(3,:);
    generalized_zb=expData(4,:);
    %% Execute the motor commands and receive end-effector trajectory from EM tracker

    
    %Calculate the reward function
    for i=1:200
        Q(1:(end+1-i),1) = Q(1:(end+1-i),1) + exp(-abs(sqrt((x(i)-x0)^2+(y(i)-y0)^2)-30)+abs(z(i)-z0));
        plot3(x(i),y(i),z(i),'ro','MarkerSize',6);
        hold on;
        pause(0.025);
    end
    %Normalize the Q values
    Q(:,iter+1) = Q(:,iter+1)./200;
    disp(['Number of Iteration=', num2str(iter), 'Return of the Rollout=', num2str(Q(1,iter+1))]);
end

%Calculate the return of the final rollout
Return(iter+1) = Q(1,iter+1);
%Plot the return over the rollouts
figure(1);
plot(Return);
ylabel('Return');
xlabel('Rollouts');
disp(['Return of the Final Rollout=', num2str(Return(end))]);