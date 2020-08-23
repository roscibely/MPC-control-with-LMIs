%Simulation
clc, clear
N=60;
[xp,xop, rep, imp, gammap, up] = park2011();
[xr,xor, rer, imr, gammar, ur] = Rego2018();
[xk,xok, rek, imk, gammak, uk] = kim2017();
t=0:1:N;
%sates
figure(1);
subplot(211); stairs(t, xp(1,:),'r'); hold on;
stairs(t, xr(1,:),'b-.');  hold on;
stairs(t, xk(1,:),'c--');
leg= legend('${x}_{Park}$', '${x}_{proposed}$', '${x}_{Kim}$');
set(leg,'Interpreter','latex');
ylabel('State 1'); xlabel('Time step');
subplot(212); stairs(t, xp(2,:),'r'); hold on;
stairs(t, xr(2,:),'b-.');  hold on;
stairs(t, xk(2,:),'c--');
leg= legend('${x}_{Park}$', '${x}_{proposed}$', '${x}_{Kim}$');
set(leg,'Interpreter','latex');
ylabel('State 2'); xlabel('Time step');
%control signal
figure(2);
stairs(t, up, 'g--'); hold on; stairs(t, ur, 'k-'); hold on;  stairs(t, uk, 'r-.');
legend('u_{Park}', 'u_{proposed}', 'u_{Kim}');
ylabel('u(k)'); xlabel('Time step');
%upper bound
figure(3);
stairs(gammap, 'g--'); hold on; stairs(gammar, 'k-'); hold on;  stairs(gammak, 'r-.');
legend('u_{Park}', 'u_{proposed}', 'u_{Kim}');
ylabel('\zeta(k)'); xlabel('Time step');
%poles
figure(4);
plot(rep',imp','kx'); hold on; plot(rer',imr','ro'); hold on; plot(rek',imk','b*');
legend('{Park}', '{proposed}', '{Kim}'); zgrid;
figure(5);
subplot(211); stairs(t, xop(1,:),'r'); hold on;
stairs(t, xor(1,:),'b-.');  hold on;
stairs(t, xok(1,:),'c--');
leg= legend('$\hat{x}_{Park}$', '$\hat{x}_{proposed}$', '$\hat{x}_{Kim}$');
set(leg,'Interpreter','latex');
ylabel('State 1'); xlabel('Time step');
subplot(212); stairs(t, xop(2,:),'r'); hold on;
stairs(t, xor(2,:),'b-.');  hold on;
stairs(t, xok(2,:),'c--');
leg= legend('$\hat{x}_{Park}$', '$\hat{x}_{proposed}$', '$\hat{x}_{Kim}$');
set(leg,'Interpreter','latex');
ylabel('State 2'); xlabel('Time step');

