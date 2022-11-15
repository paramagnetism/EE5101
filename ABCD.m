clc
clear
a = 0; b = 1; c = 6; d = 7;
Mf = 2.14 + c/20; Hf = 0.18;
Mr = 5.91 - b/10; Hr = 0.161;
Mc = 1.74; Hc = 0.098;
LFf = 0.05; LF = 0.133;
Lr = 0.128; LR = 0.308 +(a-d)/100;
Lc = 0.259; g = 9.8;
Jx = 0.5+(c-d)/100; mux = 3.333 - b/20 + a*c/60;
alpha = 15.5 -a/3 + b/2; beta = 27.5 - d/2;
gama = 11.5+(a-c)/(b+d+3); delta = 60 +(a-b)*c/10;
den = Mf*Hf*Hf+Mr*Hr*Hr+Mc*Hc*Hc+Jx;
a51 = -Mc*g/den; a52 = (Mf*Hf+Mr*Hr+Mc*Hc)*g/den;
a53 = (Mr*Lr*LF+Mc*Lc*LF+Mf*LFf*LR)*g/(LR+LF)*den;
a54 = -Mc*Hc*alpha/den; a55 = -mux/den; a56 = Mf*Hf*LFf*gama/den;
b51 = Mc*Hc*beta/den; b52 = -Mf*Hf*LFf*delta/den;
A = [0 0 0 1 0 0
    0 0 0 0 1 0
    0 0 0 0 0 1
    0 6.5 -10 -alpha 0 0
    a51 a52 a53 a54 a55 a56
    5 -3.6 0 0 0 -gama];
B = [0 0
    0 0 
    0 0 
    beta 11.2
    b51 b52
    40 delta];
C = eye(6);
D = [0];
t = 0:0.02:10;
x0 = [0.2 -0.1 0.15 -1 0.8 0];
para = [-0.5+(a-b)/20; 0.1+(b-c)/(a+d+10)];
ysp = -0.1*C(1:3,:)/A*B*para;
save('ABCD.mat','A','B','C','D','t','x0','ysp')