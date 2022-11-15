%% Initial parameters
clc
clear
load('../ABCD.mat');
C = C(1:3,:);
disturbance = [-1; 1];
%% reachability
ctrlable = rank([A B
    C zeros(3,2)]) - rank(A);
%% Get ABCD for LQR
C_ = [C(1,:); C(3,:)];
assert(rank([A B;C_ zeros(2,2)]) == min(size([A B;C_ zeros(2,2)])));
Abar=[A zeros(6,2);-C_ zeros(2,2)];
Bbar=[B;zeros(2,2)];
Cbar = [C_ zeros(2,2)];
%% LQR for new
Q = diag([1 1 1.5 1.5 1 1 1 1]); R = eye(2);
gamma=[Abar -Bbar/R*Bbar'
    -Q -Abar'];
[vector,value]=eig(gamma);
value=sum(value);
vec=vector(:,find(real(value)<0));
P=vec(9:16,:)/vec(1:8,:);
K_=real(inv(R)*Bbar'*P);
%% y_sp can be arbitary
y_sp = [3 -5];