% plot retrn over 10 runs 
clear all; close all;
R_all=[];
for i=1:10
ID=strcat('cr_', num2str(i), '.mat');  
tr{i}= load(ID);
R_all=[R_all;tr{i}.Return];
end
R_mean=mean(R_all);
R_min=min(R_all);
R_max=max(R_all);
index=1:1:length(R_max);
X=[index,index(end:-1:1),1];
Y=[R_max,R_min(end:-1:1),R_max(1)];
v=[X',Y'];
% v = [0 0; 1 0; 1 1; 0 1];
f = 1:1:length(X);
f1=figure(1) % ds1-6
f1.Units='normalized';
f1.Position=[0.1,0.1,0.5,0.4];

plot(index,R_mean,'b','LineWidth',2)
hold on
patch('Faces',f,'Vertices',v,'FaceColor','blue','FaceAlpha',.2,'EdgeColor','none')
axis([1 2000 0.5 1])
xlabel('Rollouts')
ylabel('Return')