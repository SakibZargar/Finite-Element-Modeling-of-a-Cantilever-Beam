clc
clear all
close all

%% Beam parameters from the problem specification :------------------------

m= 100/32.2; % lbf-s^2/ft^2
EI1= (5e6)/144; EI2= 2*EI1; EI3= EI1;  EI4= 2*EI1;      %lbf-ft^2
l= 3/12;     % ft
k= 2400;     % lbf/ft

%% The 1st step is to define the shape functions  :------------------------

phi_1= @(zeta) 1-(3*zeta^2)+(2*zeta^3);
phi_2= @(zeta) (l*zeta)-(2*l*zeta^2)+(l*zeta^3);
phi_3= @(zeta) (3*zeta^2)-(2*zeta^3)  ;
phi_4= @(zeta) -(l*zeta^2)+(l*zeta^3) ;

%% Element mass matrices      (M1, M2, M3, M4)    :------------------------

M1= (m*l/420)*[156 22*l 54 -13*l;22*l 4*l^2 13*l -3*l^2;54 13*l 156 -22*l;-13*l -3*l^2 -22*l 4*l^2];
M2= M1;
M3= M1;
M4= M1;

%% Element stiffness matrices (K1, K2, K3, K4)    :------------------------

K1= (EI1/(l^3))*[12 6*l -12 6*l;6*l 4*l^2 -6*l 2*l^2;-12 -6*l 12 -6*l;6*l 2*l^2 -6*l 4*l^2];
K2= 2*K1;
K3= 1*K1;
K4= 2*K1;
 
%% Element load vectors       (Q1, Q2, Q3, Q4)    :------------------------

Q1= zeros(4,1);               
Q2= zeros(4,1);
Q3= zeros(4,1);
Q4_p0= @(p0) p0*[l/2;(l^2)/12;l/2;-(l^2)/12];     % due to the u.d.load
Q4_P = @(P) -1*P*[phi_1(2/3);phi_2(2/3);phi_3(2/3);phi_4(2/3)]; % due to the point load 
syms p0 P
Q4(p0,P)= Q4_p0(p0)+Q4_P(P);

%% Modifications to K3 & K4 due to the spring    :-------------------------

K3_spring= k*[phi_1(1/3);phi_2(1/3);phi_3(1/3);phi_4(1/3)]*[phi_1(1/3) phi_2(1/3) phi_3(1/3) phi_4(1/3)];
K3_mod   = K3+K3_spring;
          
K4_spring= k*[phi_1(1/2);phi_2(1/2);phi_3(1/2);phi_4(1/2)]*[phi_1(1/2) phi_2(1/2) phi_3(1/2) phi_4(1/2)];
K4_mod   = K4+K4_spring;

%% Transformation matrices to get the global [K], [M] and [Q] :------------
 
A1= [eye(4) zeros(4,6)];
A2= [zeros(4,2) eye(4) zeros(4,4)];
A3= [zeros(4,4) eye(4) zeros(4,2)];
A4= [zeros(4,6) eye(4)];

%% Assembling the element [M],[K],[Q] to find the global [M],[K],[Q] :-----

M = ((A1'*M1*A1)+(A2'*M2*A2)+(A3'*M3*A3)+(A4'*M4*A4));
K = (A1'*K1*A1)+(A2'*K2*A2)+(A3'*K3_mod*A3)+(A4'*K4_mod*A4);
Q = (A1'*Q1)+(A2'*Q2)+(A3'*Q3)+(A4'*Q4);

%% Since the left end is fixed, q1=q2=0, so the 1st two rows and columns 
% of the global [M],[K],[Q] matrices are elliminated            :----------

M_red = M(3:10,3:10);
K_red = K(3:10,3:10);

%% State-space representation of the equations of motion :-----------------

%% Part (1):displacement(@A) p0= 0, P= 500 lbf (0<=t<=0.01), P= 0 (t>0.01)

F = Q(0,1); F = eval(F(3:10));
AA= [zeros(8,8) eye(8);-inv(M_red)*K_red zeros(8,8)];
BB= [zeros(8,1);inv(M_red)*F];
CC= [zeros(1,6) 1 zeros(1,9)];
DD= zeros(1,1);

time= transpose(0:0.000001:0.3);
pforce= zeros(length(time),1);
pforce(1:10001)=500;
SYS_part_a= ss(AA,BB,CC,DD);
[u,t]= lsim(SYS_part_a,pforce,time);
plot(t,u,'Linewidth',2)
xlabel('time (s)'); 
ylabel('displacement of point A (ft)');
title('p0=0 and P= Impulse of 500 lbf for 0.01s');
grid on;

%% Part (1):Bending moment(@C) p0= 0, P= 500 lbf (0<=t<=0.01), P= 0 (t>0.01)

CC1= [1 zeros(1,15)]; S1= ss(AA,BB,CC1,DD); [u1,t]= lsim(S1,pforce,time);
CC2= [0 1 zeros(1,14)]; S2= ss(AA,BB,CC2,DD); [u2,t]= lsim(S2,pforce,time);
CC3= [0 0 1 zeros(1,13)]; S3= ss(AA,BB,CC3,DD); [u3,t]= lsim(S3,pforce,time);
CC4= [0 0 0 1 zeros(1,12)]; S4= ss(AA,BB,CC4,DD); [u4,t]= lsim(S4,pforce,time);
CC5= [0 0 0 0 1 zeros(1,11)]; S5= ss(AA,BB,CC5,DD); [u5,t]= lsim(S5,pforce,time);
CC6= [0 0 0 0 0 1 zeros(1,10)]; S6= ss(AA,BB,CC6,DD); [u6,t]= lsim(S6,pforce,time);
CC7= [0 0 0 0 0 0 1 zeros(1,9)]; S7= ss(AA,BB,CC7,DD); [u7,t]= lsim(S7,pforce,time);
CC8= [0 0 0 0 0 0 0 1 zeros(1,8)]; S8= ss(AA,BB,CC8,DD); [u8,t]= lsim(S8,pforce,time);

q= [zeros(300001,1) zeros(300001,1) u1 u2 u3 u4 u5 u6 u7 u8]';
BM_c= K(2,:)*q;
plot(t,BM_c,'Linewidth',2)
xlabel('time (s)'); 
ylabel('Bending moment at point C (lbf-ft)');
title('p0=0 and P= Impulse of 500 lbf for 0.01s');
grid on;

%% Part (2):displacement(@A) p0= 0, P= 500 sin(2*pi*10*t)

F = Q(0,1); F = eval(F(3:10));
BB= [zeros(8,1);inv(M_red)*F];        % AA, CC and DD will remain the same

time= transpose(0:0.000001:0.3);
pforce=500*sin(20*pi*time);
SYS_part_a= ss(AA,BB,CC,DD);
[u,t]= lsim(SYS_part_a, pforce, time);
plot(t,u,'Linewidth',2)
xlabel('time (s)'); 
ylabel('displacement of point A (ft)');
title('p0= 0 and P= 500 sin(2*pi*10*t)');
grid on;

%% Part (2):Bending moment(@C) p0= 0, P= 500 sin(2*pi*10*t)

CC1= [1 zeros(1,15)]; S1= ss(AA,BB,CC1,DD); [u1,t]= lsim(S1,pforce,time);
CC2= [0 1 zeros(1,14)]; S2= ss(AA,BB,CC2,DD); [u2,t]= lsim(S2,pforce,time);
CC3= [0 0 1 zeros(1,13)]; S3= ss(AA,BB,CC3,DD); [u3,t]= lsim(S3,pforce,time);
CC4= [0 0 0 1 zeros(1,12)]; S4= ss(AA,BB,CC4,DD); [u4,t]= lsim(S4,pforce,time);
CC5= [0 0 0 0 1 zeros(1,11)]; S5= ss(AA,BB,CC5,DD); [u5,t]= lsim(S5,pforce,time);
CC6= [0 0 0 0 0 1 zeros(1,10)]; S6= ss(AA,BB,CC6,DD); [u6,t]= lsim(S6,pforce,time);
CC7= [0 0 0 0 0 0 1 zeros(1,9)]; S7= ss(AA,BB,CC7,DD); [u7,t]= lsim(S7,pforce,time);
CC8= [0 0 0 0 0 0 0 1 zeros(1,8)]; S8= ss(AA,BB,CC8,DD); [u8,t]= lsim(S8,pforce,time);

q= [zeros(300001,1) zeros(300001,1) u1 u2 u3 u4 u5 u6 u7 u8]';
BM_c= K(2,:)*q;
plot(t,BM_c,'Linewidth',2)
xlabel('time (s)'); 
ylabel('Bending moment at point C (lbf-ft)');
title('p0= 0, P= 500 sin(2*pi*10*t)');
grid on;

%% Part (3): displacement(@A) p0= 100 lbf/in=1200 lbf/ft, P= 0

F = Q(1,0); F = eval(F(3:10));
BB= [zeros(8,1);inv(M_red)*F];        % AA, CC and DD will remain the same

time= transpose(0:0.000001:0.2);
pforce=1200*ones(length(time),1);
SYS_part_a= ss(AA,BB,CC,DD);
[u,t]= lsim(SYS_part_a, pforce, time);
plot(t,u,'Linewidth',2)
xlabel('time (s)'); 
ylabel('displacement of point A (ft)');
title('p0= 1200 lbf/ft and P= 0');
grid on;

%% Part (3): displacement(@B) p0= 100 lbf/in=1200 lbf/ft, P= 0

F = Q(1,0); F = eval(F(3:10));
BB= [zeros(8,1);inv(M_red)*F];        % AA and DD will remain the same
CC= [1 zeros(1,15)];
time= transpose(0:0.000001:0.2);
pforce=1200*ones(length(time),1);
SYS_part_a= ss(AA,BB,CC,DD);
[u,t]= lsim(SYS_part_a, pforce, time);
plot(t,u,'Linewidth',2)
xlabel('time (s)'); 
ylabel('displacement of point B (ft)');
title('p0= 1200 lbf/ft and P= 0');
grid on;

