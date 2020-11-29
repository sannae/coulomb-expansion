clear all
close all
clc

colordef black

N=1e4;

R=1;

e=1;
QQ=1;
me=1;
M=100;
m1=2*M;
m2=10*M;

v=0.2665;

% elettroni
x=-1+2*rand(N,1);
y=-1+2*rand(N,1);
z=-1+2*rand(N,1);

vx=v*randn(N,1);
vy=v*randn(N,1);
vz=v*randn(N,1);

dentroe=x.^2+y.^2+z.^2<=R^2;
N0e=sum(dentroe);

qe=e*QQ/N0e;
QMe=-e/me*ones(N,1);
Qe=-qe*ones(N,1);

xe=x(dentroe);
ye=y(dentroe);
ze=z(dentroe);

vxe=vx(dentroe);
vye=vy(dentroe);
vze=vz(dentroe);

QMe=QMe(dentroe);
Qe=Qe(dentroe);

% ioni
x=-1+2*rand(N,1);
y=-1+2*rand(N,1);
z=-1+2*rand(N,1);

vx=v*zeros(N,1);
vy=v*zeros(N,1);
vz=v*zeros(N,1);

dentroi=x.^2+y.^2+z.^2<=R^2;
dentro1=x.^2+y.^2+z.^2<(R/2)^2;
dentro2=logical(dentroi-dentro1);
N0i=sum(dentroi);


x1=x(dentro1);
y1=y(dentro1);
z1=z(dentro1);    

qi=e*QQ/N0i;
QM1=e/m1*ones(N,1);
Qi=qi*ones(N,1);

vx1=vx(dentro1);
vy1=vy(dentro1);
vz1=vz(dentro1);

QM1=QM1(dentro1);
Qi=Qi(dentroi);

x2=x(dentro2);
y2=y(dentro2);
z2=z(dentro2);    

QM2=e/m2*ones(N,1);

vx2=vx(dentro2);
vy2=vy(dentro2);
vz2=vz(dentro2);

QM2=QM2(dentro2);

x=[xe;x1;x2];
y=[ye;y1;y2];
z=[ze;z1;z2];

vx=[vxe;vx1;vx2];
vy=[vye;vy1;vy2];
vz=[vze;vz1;vz2];

QM=[QMe;QM1;QM2];
Q=[Qe;Qi];

T=100;
dt=0.1;
t=[0:dt:T];

% fig=figure;
% set(fig,'color','k');
% mov = avifile('el_ion_density.avi','fps',25,'quality',100,'compression','Indeo5');

for it=1:length(t)
    t(it);
    r2=x.^2+y.^2+z.^2;
    [r2,ord]=sort(r2);
    x=x(ord);
    y=y(ord);
    z=z(ord);
    vx=vx(ord);
    vy=vy(ord);
    vz=vz(ord);
    QM=QM(ord);
    Q=Q(ord);
    r=sqrt(r2);
    E=e*cumsum(Q)./r2;
    vx=vx+dt*QM.*E.*x./r;
    vy=vy+dt*QM.*E.*y./r;
    vz=vz+dt*QM.*E.*z./r;
    x=x+vx*dt;
    y=y+vy*dt;
    z=z+vz*dt;
 
   figure(1)
    %plot3(x(Q<0),y(Q<0),z(Q<0),'.b','Markersize',3)
    %axis([-2 2 -2 2 -2 2])
    %hold on
    plot3(x(QM==e/m2),y(QM==e/m2),z(QM==e/m2),'.g','Markersize',3)
    axis([-5 5 -5 5 -5 5])
    hold on
    plot3(x(QM==e/m1),y(QM==e/m1),z(QM==e/m1),'.r','Markersize',3)
    axis([-5 5 -5 5 -5 5])
    hold off
    


   
end





