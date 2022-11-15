%% Initial parameters
clc
clear
load('../ABCD.mat');
%% plot eye
Q1 = eye(6);
Q2 = Q1*10;
Q3 = diag([1 1 4 5 1 4]);
Q4 = diag([1 1 9 1 1 1]);
R = eye(2);

K1 = getK(Q1,R); K2 = getK(Q2,R); K3 = getK(Q3,R); K4 = getK(Q4,R);
A_BK1 = A-B*K1; A_BK2 = A-B*K2; A_BK3 = A-B*K3; A_BK4 = A-B*K4;
sys1 = ss(A_BK1, B, C, D); sys2 = ss(A_BK2, B, C, D);
sys3 = ss(A_BK3, B, C, D); sys4 = ss(A_BK4, B, C, D);
initial(sys1,'r',sys2,'b',sys3, 'g', sys4, 'y',x0, t);
xlabel('t'),ylabel('Amplitude'),legend('eye(6)','eye(6)*10','114514', '119111')
%% plot feedback
[y1, t, x1]=initial(sys1, x0, t); u1 = -K1*x1';
[y2, t, x2]=initial(sys2, x0, t); u2 = -K2*x2';
[y3, t, x3]=initial(sys3, x0, t); u3 = -K3*x3';
[y4, t, x4]=initial(sys4, x0, t); u4 = -K4*x4';
figure(1)
plot(t,u1(1,:),t,u2(1,:),t,u3(1,:),t,u4(1,:));
legend('eye(6)','eye(6)*10','114514', '119111');
xlabel('t'),ylabel('Feedback1'),
figure(2)
plot(t,u1(2,:),t,u2(2,:),t,u3(2,:),t,u4(2,:));
legend('eye(6)','eye(6)*10','114514', '119111');
xlabel('t'),ylabel('Feedback2'),
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