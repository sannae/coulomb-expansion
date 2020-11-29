clear all
close all
clc

colordef black

N=1e4;

R=1;

e=1;
QQ=1;
me=1;
M=1836;
m1=2*M;
m2=3*M;
m3=(m1+m2)/2;


v=0.5;

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
N0i=sum(dentroi);
l=length(dentroi);

if mod(l,2)==0
    M1=m1*ones(l/2,1);
    M2=m2*ones(l/2,1);
else
    M1=m1*ones(fix(l/2),1);
    M2=m2*ones(fix(l/2)+1,1);
end

mi=[M1;M2];
h=randperm(l);
mi=mi(h);

qi=e*QQ/N0i;
QMi=e./mi;
Qi=qi*ones(N,1);

xi=x(dentroi);
yi=y(dentroi);
zi=z(dentroi);

vxi=vx(dentroi);
vyi=vy(dentroi);
vzi=vz(dentroi);

QMi=QMi(dentroi);
Qi=Qi(dentroi);

x=[xe;xi];
y=[ye;yi];
z=[ze;zi];

vx=[vxe;vxi];
vy=[vye;vyi];
vz=[vze;vzi];

QM=[QMe;QMi];
Q=[Qe;Qi];

T=50;
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
    plot3(x(Q<0),y(Q<0),z(Q<0),'.b','Markersize',3)
    axis([-2 2 -2 2 -2 2])
    hold on
    plot3(x(QM==e/m2),y(QM==e/m2),z(QM==e/m2),'.g','Markersize',3)
    axis([-2 2 -2 2 -2 2])
    hold on
    plot3(x(QM==e/m1),y(QM==e/m1),z(QM==e/m1),'.r','Markersize',3)
    axis([-2 2 -2 2 -2 2])
    hold off
    


   
end





