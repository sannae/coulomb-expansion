clear all
close all
clc

colordef black

N=1e4;

R=1;

e=1;
QQ=1;
me=1;
mi=50;

v=0.01;

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

qi=e*QQ/N0i;
QMi=e/mi*ones(N,1);
Qi=qi*ones(N,1);

xi=x(dentroi);
yi=y(dentroi);
zi=z(dentroi);

vxi=vx(dentroi);
vyi=vy(dentroi);
vzi=vz(dentroi);

QMi=QMi(dentroi);
Qi=Qi(dentroi);

dt=0.1;
T=10;
t=[0:dt:T];

for it=1:length(t)
    t(it);
    
    
    %elettroni
    
    r2e=xe.^2+ye.^2+ze.^2;
    [r2e,orde]=sort(r2e);
    xe=xe(orde);
    ye=ye(orde);
    ze=ze(orde);
    vxe=vxe(orde);
    vye=vye(orde);
    vze=vze(orde);
    QMe=QMe(orde);
    Qe=Qe(orde);
    re=sqrt(r2e);
    Ee=e*cumsum(Qe)./r2e;
    vxe=vxe+dt*QMe.*Ee.*xe./re;
    vye=vye+dt*QMe.*Ee.*ye./re;
    vze=vze+dt*QMe.*Ee.*ze./re;
    xe=xe+vxe*dt;
    ye=ye+vye*dt;
    ze=ze+vze*dt;
    
    %ioni
    
    r2i=xi.^2+yi.^2+zi.^2;
    [r2i,ordi]=sort(r2i);
    xi=xi(ordi);
    yi=yi(ordi);
    zi=zi(ordi);
    vxi=vxi(ordi);
    vyi=vyi(ordi);
    vzi=vzi(ordi);
    QMi=QMi(ordi);
    Qi=Qi(ordi);
    ri=sqrt(r2i);
    Ei=e*cumsum(Qi)./r2i;
    vxi=vxi+dt*QMi.*Ei.*xi./ri;
    vyi=vyi+dt*QMi.*Ei.*yi./ri;
    vzi=vzi+dt*QMi.*Ei.*zi./ri;
    xi=xi+vxi*dt;
    yi=yi+vyi*dt;
    zi=zi+vzi*dt;
    
    figure(1)
plot3(xe,ye,ze,'.r','MarkerSize',3)
axis([-5 5 -5 5 -5 5])
hold on
plot3(xi,yi,zi,'.b','MarkerSize',3)
axis([-5 5 -5 5 -5 5])
hold on
[k,l,p]=sphere;
mesh ( k,l,p )
hold off

end





