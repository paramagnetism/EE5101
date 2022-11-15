%% Initial parameters
clc
clear
load('../ABCD.mat');
C = C(1:3,:);
disturbance = [-1; 1];
ysp1 = C/A*B*[2;5];
ysp2 = [1 10 2];
%% Get new ABCD for LQR
assert(rank([A B;C zeros(3,2)]) == min(size([A B;C zeros(3,2)])));
Abar=[A zeros(6,3);-C zeros(3,3)];
Bbar=[B;zeros(3,2)];
%% LQR for Kbar
Q = diag([1 1 1.5 1.5 1 1 1 1 1]); R = eye(2);
gamma=[Abar -Bbar/R*Bbar'
    -Q -Abar'];
[vector,value]=eig(gamma);
value=sum(value);
vec=vector(:,find(real(value)<0));
P=vec(10:18,:)/vec(1:9,:);
K=real(inv(R)*Bbar'*P);

%% observer pole placement
load('../pole/poles.mat')
%poles = [sols(1) sols(2) 3*real(sols(1)) 3*real(sols(1)) 4*real(sols(1)) 4*real(sols(1))];
%poles = [sols(1) sols(2) 4*real(sols(1)) 5*real(sols(1)) 6*real(sols(1)) 7*real(sols(1))];
poles = [sols(1) sols(2) 8*real(sols(1)) 8*real(sols(1)) 8*real(sols(1)) 8*real(sols(1))];
%% full order pole placement
q = [0;0;1]; B_q=C'*q;
Wc = [B_q A'*B_q A'^2*B_q A'^3*B_q A'^4*B_q A'^5*B_q];
assert(rank(Wc(:,1:6))==6);
phids = double((A'-poles(1)*eye(6))*(A'-poles(2)*eye(6))*(A'-poles(3)*eye(6))*(A'-poles(4)*eye(6))*(A'-poles(5)*eye(6))*(A'-poles(6)*eye(6)));
K_ = q*([0 0 0 0 0 1]/Wc*phids);
L = K_'; 