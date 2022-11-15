%% Initial parameters
clc
clear
load('../ABCD.mat');
C = [C(1,:); C(3,:)];
%% get sigma
for i=1:2
    if C(1,:)*A^(i-1)*B==0
       continue
    end
    sigma1=i;
    break
end
for i=1:2
    if C(2,:)*A^(i-1)*B==0
       continue
    end
    sigma2=i;
    break
end
%% getting K and F
load('../pole/poles.mat')
Bstar=[C(1,:)*A^(sigma1-1)*B
       C(2,:)*A^(sigma2-1)*B];
F = Bstar^-1;
Cstar = [C(1,:)*(A^2+2*kesi*omega*A+omega^2*eye(6))
    C(2,:)*(A^2+2*kesi*omega*A+omega^2*eye(6))];
K = F*Cstar;
%% solve
A_BK=A-B*K;
sys=ss(A_BK,B*F,C,D);
%% plot initial
figure(1)
initial(sys, x0, t)
%% plot state transfer
figure(2)
step(sys)