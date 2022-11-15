%% Initial parameters
clc
clear
load('../ABCD.mat');
%% quality limitations
% Because Mp <= 0.1
kesi = sqrt(log(0.1)^2/(pi^2+log(0.1)^2));
kesi = ceil(kesi*10)/10; % bigger than max
% Because ts = 0.02
omega = 4/(5*kesi);
syms s;
eqn = s^2+2*kesi*omega*s+omega^2 == 0;
sols = vpasolve(eqn, s);
save('poles.mat','sols','kesi','omega')
%% plot Cart position, handle angle and bike angle
K1 = poleplace(3,3,3,3); K2= poleplace(2,3,4,5); K3=poleplace(5,5,5,5);
A_BK1 = A-B*K1; A_BK2 = A-B*K2; A_BK3 = A-B*K3;
sys1 = ss(A_BK1, B, C, D); sys2 = ss(A_BK2, B, C, D); sys3 = ss(A_BK3, B, C, D);
initial(sys1,'r',sys2,'b',sys3, 'g', x0, t);
xlabel('t');ylabel('Amplitude');legend('3333','2345','5555');
%% plot feedback
[y1, t, x1]=initial(sys1, x0, t); u1 = -K1*x1';
[y2, t, x2]=initial(sys2, x0, t); u2 = -K2*x2';
[y3, t, x3]=initial(sys3, x0, t); u3 = -K3*x3';
plot(t,u1(2,:),t,u2(2,:),t,u3(2,:));
legend('3333','2345','5555');
xlabel('t'),ylabel('Feedback'),
%% pole placement
function K = poleplace(pole3, pole4, pole5, pole6)
    load('../ABCD.mat');
    load('poles.mat');
    % check AB
    if rank([B A*B A^2*B A^3*B A^4*B A^5*B]) < 6
    error('AB not controlable');
    end
    % check ABq
    q = [0;1]; Bq = B*q;
    Wc = [Bq A*Bq A^2*Bq A^3*Bq A^4*Bq A^5*Bq];
    if rank(Wc)<6
        error('ABq not controlable');
    end
    % find K
    poles = [sols(1) sols(2) pole3*real(sols(1)) pole4*real(sols(1)) pole5*real(sols(1)) pole6*real(sols(1))];
    phids = double((A-poles(1)*eye(6))*(A-poles(2)*eye(6))*(A-poles(3)*eye(6))*(A-poles(4)*eye(6))*(A-poles(5)*eye(6))*(A-poles(6)*eye(6)));
    K = q*([0 0 0 0 0 1]*(Wc^-1)* phids);
end