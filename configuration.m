function [ theta,L ] = configuration( Pn,P0,N )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%this is to compute the confuguration of the TSM robot 
% given 
% Pn: the position vector of the end effector
% P0: position vector of the base node
% N: number of nodes, N=21
% L:N-by-3 store the 3D coordination of each node 
l=3.2; %length of each node
options = optimset('Display','none','TolFun',1e-8,'MaxFunEvals',1000);
%eq_x=@(ang)Pn(2)-(l*sin(N*ang/2)/sin(ang/2)*sin((N+1)*ang/2));
eq_x=@(ang)Pn(2)-l*sin(N*ang/2)/sin(ang/2)*sin((N+1)*ang/2)-P0(2);
theta=fsolve(eq_x,[1e-3],options);

L=zeros(N,3);
L(1,1)=l*cos(theta)+P0(1);
L(1,2)=l*sin(theta)+P0(2);
L(1,3)=Pn(3)+P0(3);
for i=2:N
    L(i,1)=L(i-1,1)+l*cos(theta*i);
    L(i,2)=L(i-1,2)+l*sin(theta*i);
    L(i,3)=Pn(3);
end
end

