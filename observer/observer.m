%% Initial parameters
clc
clear
load('../ABCD.mat');
C = C(1:3,:);
%% pole placement
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
%% LQR for K
R = eye(2);
K = getK(diag([1 1 9 1 1 1]),R);

%% LQR
function K = getK(Q, R)
load('../ABCD.mat');
gamma=[A -B/R*B'
    -Q -A'];
[vector,value]=eig(gamma);
value=sum(value);
vec=vector(:,find(real(value)<0));
P=vec(7:12,:)/vec(1:6,:);
K=real(inv(R)*B'*P);
end
