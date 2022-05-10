function [xx,yy]=elipse_matrix(C,num)
[U,L] = eig(C); %separa os autovalores e autovetores
M = [0 0]';
N = 1; %desvio unitário, sem distorção
radii = N*sqrt(diag(L)); %monta uma matriz diagonal extraindo-se a raíz
theta = linspace(0,2*pi,num); % circulo unitário
p(1,:) = (radii(1))*cos(theta); % resultado da raiz do autovalor 1 com cos
p(2,:) = (radii(2))*sin(theta);% resultado da raiz do autovalor 2 com sen
p = U*p; %multiplica este resultado com o respectivo autovetor
% plot(p(1,:),p(2,:),'g--','LineWidth',2)
% plot(0,0,'r+')
xx = p(1,:);
yy = p(2,:);