close all
clear all
clc

%%%% Domain discretization
dt=0.5;
Nx=100;%number of points in each side
maxX=50;%maximum
dx=maxX/Nx;
Ny=100;%number of points in each side
maxY=50;%maximum
dy=maxY/Ny;
Tmax=251;
T=0:dt:Tmax;
X=0:dx:maxX; 
Y=0:dy:maxY; 

%%%% g(x,y)
g=zeros(length(X),length(Y));
for x=10:20
    for y=20:30
        g(x,y)=1;
    end
end

%%%% Parameters
beta=0.2;
gamma=0.08;
D11=0.8; 
D12= 0; %% cross-Difusion
D21= 0; %% cross-Difusion
D22=0.4; 

%%%% Initial condition
U=zeros(2,length(X),length(Y),length(T));
Itot=zeros(1,length(T));
Rtot=zeros(1,length(T));
Stot=zeros(1,length(T));
Norm=zeros(2,length(T));

U(2, :, :, 1)=10^(-2)*g;
U(1, :, :, 1)=1-U(2, :, :, 1);

Itot1=sum(U(2, :, :, 1));
Itot(1,1)=sum(Itot1);
Stot1=sum(U(1, :, :, 1));
Stot(1,1)=sum(Stot1);
Rtot1=sum(1-U(1, :, :, 1)-U(2, :, :, 1));
Rtot(1,1)=sum(Rtot1);

%%%% Matrices
D=[D11 D12; D21 D22];
phi=dt/(2*dx^2)*D;
I=eye(2);

%%%% Iterações 

for n=1:length(T)-1

  %%%% Calcula U~
  U_tilde=zeros(2,length(X),length(Y));

%   %%% first order
%   U_tilde=U(:,:,:,n);

  %%% second order
  %%%% Calcula R
    R=SIR_Reaction(U(:,:,:,n), beta, gamma);

  U_u=zeros(length(X), length(Y));
  R_u=zeros(length(X), length(Y));
  U_v=zeros(length(X), length(Y));
  R_v=zeros(length(X), length(Y));
  U_u(:,:)=U(1, :, :, n);
  R_u(:,:)=R(1, :, :);
  U_v(:,:)=U(2, :, :, n);
  R_v(:,:)=R(2, :, :);
   

    U_tilde(1,:,:)=U_u+(dt/2)*(D11*laplacian2(dx,dy,U_u)+D12*laplacian2(dx,dy,U_v)+R_u);
    U_tilde(2,:,:)=U_v+(dt/2)*(D21*laplacian2(dx,dy,U_u)+D22*laplacian2(dx,dy,U_v)+R_v);
    
    %%%% Calcula F^n
    %%%% Calcula R
    R=SIR_Reaction(U_tilde, beta, gamma);
    Fn=zeros(2,length(X),length(Y));
    for i=1:length(X)
        for j=1:length(Y)
            if j==1
                Fn(:,i,j)=(I-2*phi)*U(:,i,j,n)+2*phi*U(:,i,j+1,n)+dt*R(:,i,j)/2;
            elseif j==length(Y)
                Fn(:,i,j)=2*phi*U(:,i,j-1,n)+(I-2*phi)*U(:,i,j,n)+dt*R(:,i,j)/2;
            else
                Fn(:,i,j)=phi*U(:,i,j-1,n)+(I-2*phi)*U(:,i,j,n)+phi*U(:,i,j+1,n)+dt*R(:,i,j)/2;
            end
        end
    end

    %%%% Calcula U^(n+1/2)
    U_n12=zeros(2, length(X),length(Y));
    for j=1:length(Y)
        U_n12(:,:,j)=max(0,Thomas_algorithm(-phi,Fn(:,:,j)));
    end

    %%%% Calcula F^(n+1/2)
    % %%%% Calcula R
    % R=SIR_Reaction(U_n12, beta, gamma);
    Fn12=zeros(2,length(X),length(Y));
    for j=1:length(Y)
        for i=1:length(X)
            if i==1
                Fn12(:,i,j)=(I-2*phi)*U_n12(:,i,j)+2*phi*U_n12(:,i+1,j)+dt*R(:,i,j)/2;
            elseif i==length(X)
                Fn12(:,i,j)=2*phi*U_n12(:,i-1,j)+(I-2*phi)*U_n12(:,i,j)+dt*R(:,i,j)/2;
            else
                Fn12(:,i,j)=phi*U_n12(:,i-1,j)+(I-2*phi)*U_n12(:,i,j)+phi*U_n12(:,i+1,j)+dt*R(:,i,j)/2;
            end
        end
    end

    %%%% Calcula U^(n+1)
    for i=1:length(X)
        U(:,i,:,n+1)=max(0,Thomas_algorithm(-phi,Fn12(:,i,:)));
    end

    Norm(1,n+1)=norm(U(1, i, j, n+1)-U_u);
    Norm(2,n+1)=norm(U(2, i, j, n+1)-U_v);

    Itoti=sum(U(2, :, :, n+1));
    Itot(1,n+1)=sum(Itoti);
    Stoti=sum(U(1, :, :, n+1));
    Stot(1,n+1)=sum(Stoti);
    Rtoti=sum(1-U(1, :, :, n+1)-U(2, :, :, n+1));
    Rtot(1,n+1)=sum(Rtoti);
end

U_u=zeros(length(X), length(Y));
U_v=zeros(length(X), length(Y));


U_u(:,:)=U(1, :, :, 2);
U_v(:,:)=U(2, :, :, 2);

figure
heatmap(U_u,'Colormap', jet);
grid off 
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clim([0, 1]);

figure
heatmap(U_v,'Colormap', jet);
grid off
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clim([0, 0.1]);

U_u(:,:)=U(1, :, :, 80);
U_v(:,:)=U(2, :, :, 80);


figure
heatmap(U_u,'Colormap', jet);
grid off 
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clim([0, 1]);

figure
heatmap(U_v,'Colormap', jet);
grid off
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clim([0, 0.1]);

U_u(:,:)=U(1, :, :, 160);
U_v(:,:)=U(2, :, :, 160);

figure
heatmap(U_u,'Colormap', jet);
grid off 
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clim([0, 1]);

figure
heatmap(U_v,'Colormap', jet);
grid off
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clim([0, 0.1]);

U_u(:,:)=U(1, :, :, 240);
U_v(:,:)=U(2, :, :, 240);

figure
heatmap(U_u,'Colormap', jet);
grid off 
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clim([0, 1]);

figure
heatmap(U_v,'Colormap', jet);
grid off
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clim([0, 0.1]);

U_u(:,:)=U(1, :, :, 320);
U_v(:,:)=U(2, :, :, 320);

figure
heatmap(U_u,'Colormap', jet);
grid off 
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clim([0, 1]);

figure
heatmap(U_v,'Colormap', jet);
grid off
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
clim([0, 0.1]);

figure 
hold on
%plot(T,Stot,'b');
plot(T,Itot,'r');
xlabel('t (days)');
ylabel('I');
title('Total of infectious over time');
