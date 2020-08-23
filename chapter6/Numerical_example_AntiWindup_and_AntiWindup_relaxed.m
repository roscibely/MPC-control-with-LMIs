%%
% Anti windup based  Turner, Herrmann and Postlethwaite (2003) 
% Author: Rego, R. C. B.
% Year: 2019
%%
clear; clc;
%% Model parameters
B = [0.0935; 0.00478];
C = [0.333 -1];
D = [0];
N = 60;
alf=[1 1];
bet=[0.1 1];
alpha=alf(2)*rand(1,100) +alf(1);
beta=bet(2)*rand(1,100)+bet(1);
A1 = [0.872 -0.0623*alf(1); 0.0935 0.997];
A2 = [0.872 -0.0623*alf(2); 0.0935 0.997];
B1= bet(1)*B;
B2=bet(2)*B;
%Weighting matrix
Le = 1*eye(2);
R = 1;
%Constrain
umax = 1;
%Initial states
x = [-1.5; -0.2]; %initial 
xo = [-0.5; 1]; %observer
u=1;
%% 
n=size(B1,1); m=size(B1,2);
Qa=sdpvar(n,n,'symmetric');
La =  sdpvar(m,n, 'full');
Ua = sdpvar(m,m,  'symmetric');
xp = sdpvar(n,1);
mu=sdpvar(1);
 ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5);
%x=[0;0];
LMI=[[-Qa -La' zeros(n,m) (C*Qa+D*La)' (A1*Qa+B1*La)';
    -La -2*Ua eye(m) (D*Ua)' (B1*Ua)';
    zeros(m,n) eye(m) -mu*eye(m) zeros(m,n) zeros(m,m);
    (C*Qa+D*La) (D*Ua) zeros(m,m) -eye(m) zeros(m,n);
    (A1*Qa+B1*La) (B1*Ua) zeros(n,m) zeros(n,m) -Qa]<=0];
LMI=[LMI, [-Qa -La' zeros(n,m) (C*Qa+D*La)' (A2*Qa+B2*La)';
    -La -2*Ua eye(m) (D*Ua)' (B2*Ua)';
    zeros(m,n) eye(m) -mu*eye(m) zeros(m,n) zeros(m,m);
    (C*Qa+D*La) (D*Ua) zeros(m,m) -eye(m) zeros(m,n);
    (A2*Qa+B2*La) (B2*Ua) zeros(n,m) zeros(n,m) -Qa]<=0];
AW = optimizer(LMI,mu,ops,xp,{Qa,La,mu});
sol = AW{[0;0]};
Fa = sol{2}*inv(sol{1});
disp('Ganho AW'); fprintf('%f ', Fa);
%% With relaxation
Xa=sdpvar(n,n, 'full');
Qa=sdpvar(n,n,'symmetric');
La =  sdpvar(m,n, 'full');
Ua = sdpvar(m,m,  'symmetric');
xp = sdpvar(n,1);
mua=sdpvar(1);
LMI1=[[-(Xa+Xa'-Qa) -La' zeros(n,m) (C*Xa+D*La)' (A1*Xa+B1*La)';
    -La -2*Ua eye(m) (D*Ua)' (B1*Ua)';
    zeros(m,n) eye(m) -mua*eye(m) zeros(m,n) zeros(m,m);
    (C*Xa+D*La) (D*Ua) zeros(m,m) -eye(m) zeros(m,n);
    (A1*Xa+B1*La) (B1*Ua) zeros(n,m) zeros(n,m) -Qa]<=0];
LMI1=[LMI1, [-(Xa+Xa'-Qa) -La' zeros(n,m) (C*Xa+D*La)' (A2*Xa+B2*La)';
    -La -2*Ua eye(m) (D*Ua)' (B2*Ua)';
    zeros(m,n) eye(m) -mua*eye(m) zeros(m,n) zeros(m,m);
    (C*Xa+D*La) (D*Ua) zeros(m,m) -eye(m) zeros(m,n);
    (A2*Xa+B2*La) (B2*Ua) zeros(n,m) zeros(n,m) -Qa]<=0];
AW = optimizer(LMI1,mua,ops,xp,{Xa,La,mua});
sol = AW{[0;0]};
Fa1 = sol{2}*inv(sol{1});
disp('Ganho AW'); fprintf('%f ', Fa1);
x1=x;
u1=u;
for i = 1:N
A = [0.872 -0.0623*alpha(i); 0.0935 0.997];
x(:,i+1) = A*x(:,i)+beta(i)*B*u(i);
y(:,i)= C*x(:,i);
u(i+1) = Fa*x(:,i+1);
if(u(i+1)>umax) u(i+1)=umax; end
x1(:,i+1) = A*x1(:,i)+beta(i)*B*u1(i);
y1(:,i)= C*x1(:,i);
u1(i+1) = Fa1*x1(:,i+1);
if(u1(i+1)>umax) u1(i+1)=umax; end
end
t=0:1:N;
figure(1);
subplot(211); plot(t, x(1,:),'k','LineWidth',1.5); hold on;  plot(t, x1(1,:),'b-.','LineWidth',1.5);
subplot(212); plot(t, x(2,:),'k','LineWidth',1.5);  hold on;  plot(t, x1(2,:),'b-.','LineWidth',1.5);
figure(2); plot(t, u,'k','LineWidth',1.5); hold on;  plot(t, u1,'b-.','LineWidth',1.5);
