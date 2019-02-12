function demo1
%
%Learning and reproduction of a movement through a mixture of dynamical
%systems, where variability and correlation information along the movement
%and among the different examples is encapsulated as a full stiffness
%matrix in a set of mass-spring-damper systems. 
%For each primitive (or state), learning of the virtual attractor points 
%and associated stiffness matrices is done through least-squares regression.     
%

%% Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbData = 200; %Length of each trajectory
nbVar = 2; %Number of variables (Trajectory in a plane)
nbStates = 8; %Number of states (or primitives)

kPmax = 200; %Maximum stiffness gain
kPmin = 100; %Minimum stiffness gain
kP = 150; %Initial stiffness gain
kV = 20; %Damping gain
dt = .01; %Time step
alpha = 2.0; %Decay factor
%Gaussians equally distributed in time
Mu_t = linspace(0,nbData*dt,nbStates);
Sigma_t = (nbData*dt/nbStates)*8E-2;

%% Load dataset 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Data is composed of 3 concatenated demonstrations (3 trajectories in 2D with 200 time steps each).
%Data(posId,:), Data(velId,:) and Data(accId,:) correspond to position, velocity and acceleration variables.
posId=[1:nbVar]; velId=[nbVar+1:2*nbVar]; accId=[2*nbVar+1:3*nbVar]; 
load('data/datadmp02.mat'); 

%% Batch learning
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Compute weights - Eq.(3)
s = 1; %Initialization of decay term
for n=1:nbData
  s = s + (-alpha*s)*dt; %Update of decay term
  t = -log(s)/alpha; %Corresponding time (t=n*dt)
  for i=1:nbStates
    h(i) = gaussPDF(t,Mu_t(i),Sigma_t); %Probability to be in a given state
  end
  H(n,:) = h./sum(h); %Normalization
end
H = repmat(H,nbSamples,1); %Repete the process for each demonstration

%Batch least norm solution to find the centers of the states (or primitives) Mu_X (Y=Mu_x*H')
Y =  Data(accId,:).*(1/kP) + Data(posId,:) + Data(velId,:).*(kV/kP);
Mu_x = Y*pinv(H'); %Pseudoinverse solution Mu_x = [inv(H'*H)*H'*Y']'

%Compute residuals - Eq.(4)
RI = eye(nbVar,nbVar).*1E-3; %Regularization term for matrix inversion
%Fast computation
for i=1:nbStates
  Sigma_x(:,:,i) = cov(((Y-repmat(Mu_x(:,i),1,nbData*nbSamples))*diag(H(:,i)))');
  Wp(:,:,i) = inv(Sigma_x(:,:,i)+RI); %Use variation information to determine stiffness
end
% %Corresponding step-by-step computation (slow but easier to read)
% for i=1:nbStates
%   for j=1:nbData*nbSamples
%     Yp(:,j) = (Y(:,j)-Mu_x(:,i)) .* H(j,i);
%   end
%   for n=1:nbVar
%     Yc(n,:) = Yp(n,:) - mean(Yp(n,:));
%   end
%   covTmp=zeros(2,2);
%   for j=1:nbData*nbSamples
%     covTmp = covTmp + Yc(:,j)*Yc(:,j)';
%   end
%   Sigma_x(:,:,i) = covTmp./(nbData*nbSamples);
%   Wp(:,:,i) = inv(Sigma_x(:,:,i)+RI); %Use variation information to determine stiffness
% end;

%Rescale Wp to stay within the [kPmin,kPmax] range - Eq.(5)
for i=1:nbStates
  [V(:,:,i),Dtmp] = eig(Wp(:,:,i)); %Eigencomponents decomposition
  lambda(:,i) = diag(Dtmp); 
end
lambda_min = min(min(lambda));
lambda_max = max(max(lambda));
for i=1:nbStates
  %Rescale each eigenvalue such that they lie in the range [kPmin,kPmax]
  Dtmp = diag(((kPmax-kPmin) .* (lambda(:,i)-lambda_min)./(lambda_max-lambda_min)) + kPmin);
  Wp(:,:,i) = V(:,:,i) * Dtmp * inv(V(:,:,i)); %Reconstruction from the modified eigencomponents
end

%% Reproduction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %If desired, we can here generate variability for reproduction
% %(following the variability extracted from the training data)
% for i=1:nbStates
%   Mu_x(:,i) = Mu_x(:,i) + Sigma_x(:,:,i) * ((rand(2,1)-.5).*1E1);
% end
currPos = Data(1:nbVar,1); %Initial position
currVel = zeros(nbVar,1); %Initial velocity
s = 1; %Reinitialize the decay term
for n=1:nbData
  %Evaluate the current weights
  t = -log(s)/alpha; %Corresponding time (t=n*dt)
  for i=1:nbStates
    h(i) = gaussPDF(t,Mu_t(i),Sigma_t); %Probability to be in a given state
  end
  h = h./sum(h); %Normalization
  %Evaluate the current target and associated stiffness matrix
  currTar = zeros(nbVar,1);
  currWp = zeros(nbVar,nbVar);
  for i=1:nbStates
    currTar = currTar + Mu_x(:,i) .* h(i);
    currWp = currWp + Wp(:,:,i) .* h(i);
  end
  %Compute acceleration
  currAcc = currWp * (currTar-currPos) - kV * currVel; %Eq.(2)
  %currAcc = kP * (currTar-currPos) - kV * currVel;
  %Update veloctiy and position
  currVel = currVel + currAcc .* dt;
  currPos = currPos + currVel .* dt;
  %Update the decay term
  s = s + (-alpha*s)*dt; 
  %Keep a trace of data
  r.Data(:,n)=[currPos; currVel; currAcc];
  r.currTar(:,n)=currTar;
  r.currWp(:,:,n)=currWp;
  r.currDet(n)=det(currWp)^(1/nbVar);
end

%% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position',[10 100 1200 400]); 
%Create colormap to obtain a different color for each state 
clrmap = colormap('Jet');
xx = round(linspace(1,64,nbStates));
clrmap = clrmap(xx,:);

%Plot data and model
subplot(1,3,1); hold on; box on;
for i=1:nbStates
  plotGMM(Mu_x(1:2,i), Sigma_x(1:2,1:2,i), clrmap(i,:), 1, 3);
end
for n=1:nbSamples
 plot(Data(1,(n-1)*nbData+1:n*nbData),Data(2,(n-1)*nbData+1:n*nbData),'-','linewidth',1,'color',[.4 .4 .4]);
end
plot(r.Data(1,:),r.Data(2,:),'-','linewidth',2,'color',[0 0 0]);
plot(r.Data(1,1:3:end),r.Data(2,1:3:end),'.','markersize',16,'color',[1 0 0]);
%axis equal;
xlabel('x_1'); ylabel('x_2');

%Plot weights
subplot(1,3,2); hold on; box on;
for i=1:nbStates
  plot(H(1:nbData,i),'-','color',clrmap(i,:),'linewidth',2);
end
xlabel('t'); ylabel('h_i');

%Plot evolution of the determinant of the adaptive gain matrices
%(which reflects variability of the training set)
subplot(1,3,3); hold on; box on;
for i=1:nbStates
  plot(r.currDet,'k-','linewidth',2);
end
plot([1,nbData],[kPmin kPmin],':','color',[0 .8 0],'linewidth',2);
plot([1,nbData],[kPmax kPmax],':','color',[.8 0 0],'linewidth',2);
axis([1 nbData kPmin-20 kPmax+20]);
xlabel('t'); ylabel('|W^P|');

pause;
close all;

