%%
% CONTROLADOR MPC OFFLINE COM ATUADOR ANTI-WINDUP
% Author: Rego, R. C. B.
% Year: 2019
%%

clear; clc; close all
%% Model parameters
B = [0.0935; 0.00478];
C = [0.333 -1];
D = [0];
% alpha = 1:0.05:5;
% beta=0.1:0.01:1;
N = 10;
alf=[1 5];
bet=[0.1 1];
alpha=alf(2)*rand(1,N) +alf(1);
beta=bet(2)*rand(1,N)+bet(1);
A1 = [0.872 -0.0623*alf(1); 0.0935 0.997];
A2 = [0.872 -0.0623*alf(2); 0.0935 0.997];
% B1= beta(1)*B;
% B2=beta(N)*B;
%Weighting matrix
Le = 1*eye(2);
R = 1;
%Constrain
umax = 1;
%Initial states
x = [-1.5; -0.2]; %initial 
xo = [-0.5; 1]; %observer
u=1;
%% x_set 
xset=x;
%alpha=5*rand(1,100);
for k=1:N
A = [0.872 -0.0623*alpha(k); 0.0935 0.997];
xset(:,k+1)=A*xset(:,k);
end
 
%% Off-line robust observer design
%p =sdpvar(1,1);
p = sqrt(0.6); 
Ge = sdpvar(2,2, 'full');
Pe = sdpvar(2,2, 'symmetric');
Ye = sdpvar(2,1);
Lmi= [Pe>=0, [p^2*(Ge+Ge'-Pe)-Le (Ge*A1-Ye*C)'; Ge*A1-Ye*C Pe]>=0]; 
Lmi = [Lmi, [p^2*(Ge+Ge'-Pe)-Le (Ge*A2-Ye*C)'; Ge*A2-Ye*C Pe]>=0];
ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5);
optimize(Lmi,-trace(Ge),ops);
Lp = inv(value(Ge))*value(Ye);
e = x-xo;% estimation error 
%%  OFF-line output feedback  MPC
%LMI variables
Q = sdpvar(2,2, 'symmetric');
gamma = sdpvar(1,1);
Y0 = sdpvar(1,2, 'full');
Y1 = sdpvar(1,2, 'full');
Y2 = sdpvar(1,2, 'full');

G = sdpvar(2,2, 'full');
X = sdpvar(1,1, 'full');
xp = sdpvar(2,1);
% LPV variables
alpha_=sdpvar(1);
beta_=sdpvar(1);
Y=Y0+alpha_*Y1+beta_*Y2;
A_ = [0.872 -0.0623*alpha_; 0.0935 0.997];
B_ = beta_*B;
% constraints and objective
objective = gamma;
%optimization object
ops = sdpsettings('solver','sedumi','sedumi.eps',1e-5);
gammav(1)=45;
% Optimization with LPV variables
LA_ = [G+G'-Q G*A_'+Y'*B_' G*sqrtm(Le) Y'*sqrtm(R);
      A_*G+B_*Y Q zeros(2,2) zeros(2,1);
      sqrtm(Le)*G zeros(2,2) gamma*eye(2) zeros(2,1);
      sqrtm(R)*Y zeros(1,2) zeros(1,2) gamma*eye(1)];
L3 = [[X Y; Y' G+G'-Q]>=0, X<=umax.^2];
Gn_=1*eye(2);
Qn_=1*eye(2);
for i = 1:N
LMIs = [LA_ >= 0, L3, gamma*(G+G'-Q)-0.001*(G'*G)>=0];
r=diag(xset(:,i));
LMIs=[LMIs, [G+G'-Q (G*r)'; (G*r) Q.^2 ]>=0]; %Rset
LMIs = [LMIs, alf(1)<=alpha_<=alf(2), uncertain(alpha_)]; 
LMIs = [LMIs, bet(1)<=beta_<=bet(2), uncertain(beta_)];
LMIs = [LMIs, G-Q>=0, G>=Gn_,Q>=Qn_, gamma<=gammav(i)];
controller = optimizer(LMIs,objective,ops,xp,{G,Y0,Y1,Y2,gamma,Q});
sol= controller{{xset(:,i)}};
Gn(:,:,i) = sol{1};
Yn0(:,:,i) = sol{2};
Yn1(:,:,i) = sol{3};
Yn2(:,:,i) = sol{4};
gammav(i+1) = sol{5};
Qn(:,:,i) = sol{6};
Gn_=Gn(:,:,i);
Qn_=Qn(:,:,i);
end
N_=10;
F0 = Yn0(:,:,N_)*inv(Gn(:,:,N_));
F1 = Yn1(:,:,N_)*inv(Gn(:,:,N_));
F2 = Yn2(:,:,N_)*inv(Gn(:,:,N_));

t_=0:10;
xhset=[xset(2,:); xset(1,:)];
figure(2);
for k=1:N
[xx(k,:), yy(k,:)] = elipse_matrix(Gn(:,:,k),40); 
hold on, plot(xx(k,:),yy(k,:),'k','linewidth' ,1.5), hold on,
end
%hold on, plot(xhset(2 ,:), xhset(1 ,:) ,'r','linewidth' ,2),
set(gca,'fontsize',30,'fontname','Times New Roman')
xlabel('x_1','fontsize',30,'fontname','Times New Roman','fontangle','normal'),
ylabel('x_2', 'fontsize',30,'fontname','Times New Roman','fontangle','normal'),
title('Elipsoides de Estabilidade');
% Ellipsoids and xset
figure(10);
for k =1:N
    plot3(yy(k ,:), t_(k)*ones(40), xx(k ,:) ,'k','linewidth' ,1.5), hold on
    n_str = int2str (k);
    text ( yy(k ,1)*1.1 , t_(k),xx(k ,1)*1.1 , n_str )
end
hold on,     plot3(xhset(1 ,:), t_(1:end) , xhset(2 ,:) ,'r','linewidth' ,2), grid on;
set(gca,'fontsize',30,'fontname','Times New Roman')
xlabel('x_2','fontsize',30,'fontname','Times New Roman','fontangle','normal'),
ylabel('Tempo','fontsize',30,'fontname','Times New Roman','fontangle','normal')
zlabel('x_1','fontsize',30,'fontname','Times New Roman','fontangle','normal')
title('Elipsoides de Estabilidade com r','fontsize',30,'fontname','Times New Roman','fontangle','normal')


x1 = [-1.5; -0.2]; %initial 
alpha=alf(2)*rand(1,100) +alf(1);
beta=bet(2)*rand(1,100)+bet(1);
for i=1:100
u(i) = (F0+alpha(i)*F1+beta(i)*F2)*x1(:,i);   
disp((F0+alpha(i)*F1+beta(i)*F2));
A = [0.872 -0.0623*alpha(i); 0.0935 0.997];
E(:,i) = eig(A+(F0+alpha(i)*F1+beta(i)*F2)*beta(i)*B);
x1(:,i+1) = A*x1(:,i)+beta(i)*B*u(i);
y(:,i)= C*x1(:,i);   
end
% figure(1)
% subplot(211); plot(x1(1,:),'k-o','LineWidth',1.5);
% subplot(212); plot(x1(2,:),'k-o','LineWidth',1.5);
figure(3), 
re=real(E);
im=imag(E);
plot(re',im','ko'); zgrid;

%% 
for i = 1:N
LMIs = [LA_ >= 0, L3, gamma*(G+G'-Q)-0.001*(G'*G)>=0];
LMIs = [LMIs, alf(1)<=alpha_<=alf(2), uncertain(alpha_)]; 
LMIs = [LMIs, bet(1)<=beta_<=bet(2), uncertain(beta_)];
LMIs = [LMIs, G-Q>=0, G>=Gn_,Q>=Qn_, gamma<=gammav(i)];
controller = optimizer(LMIs,objective,ops,xp,{G,Y0,Y1,Y2,gamma,Q});
sol= controller{{xset(:,i)}};
Gn(:,:,i) = sol{1};
Yn0(:,:,i) = sol{2};
Yn1(:,:,i) = sol{3};
Yn2(:,:,i) = sol{4};
gammav(i+1) = sol{5};
Qn(:,:,i) = sol{6};
Gn_=Gn(:,:,i);
Qn_=Qn(:,:,i);
end
N_=10;
F0 = Yn0(:,:,N_)*inv(Gn(:,:,N_));
F1 = Yn1(:,:,N_)*inv(Gn(:,:,N_));
F2 = Yn2(:,:,N_)*inv(Gn(:,:,N_));

x = [-1.5; -0.2]; %initial 
alpha=alf(2)*rand(1,100) +alf(1);
beta=bet(2)*rand(1,100)+bet(1);
for i=1:100
u(i) = (F0+alpha(i)*F1+beta(i)*F2)*x(:,i);   
A = [0.872 -0.0623*alpha(i); 0.0935 0.997];
E(:,i) = eig(A+(F0+alpha(i)*F1+beta(i)*F2)*beta(i)*B);
x(:,i+1) = A*x(:,i)+beta(i)*B*u(i);
y(:,i)= C*x(:,i);   
end
figure(1),
subplot(211); plot(x(1,:),'g-+','LineWidth',1.5); hold on, plot(x1(1,:),'k-o','LineWidth',1.5); legend('x_{sem r}', 'x_{com r}')
axis([1 100 -1.5 1.5])
set(gca,'fontsize',30,'fontname','Times New Roman')
xlabel('Tempo','fontsize',30,'fontname','Times New Roman','fontangle','normal'),
ylabel('x_1','fontsize',30,'fontname','Times New Roman','fontangle','normal'), grid on,
subplot(212); plot(x(2,:),'g-+','LineWidth',1.5), hold on, plot(x1(2,:),'k-o','LineWidth',1.5); legend('x_{sem r}', 'x_{com r}')
axis([1 100 -1.5 1.5])
set(gca,'fontsize',30,'fontname','Times New Roman')
xlabel('Tempo','fontsize',30,'fontname','Times New Roman','fontangle','normal'),
ylabel('x_2','fontsize',30,'fontname','Times New Roman','fontangle','normal'), grid on;
%% 

figure(3), hold on,
re=real(E);
im=imag(E);
plot(re',im','g+'); zgrid;
