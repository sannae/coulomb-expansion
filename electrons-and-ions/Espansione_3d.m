clear all
close all
clc

colordef white

N=1e4;

R=1;

e=1;
QQ=1;
me=1;
mi=50;

v=1;

k=[-5:5];
l=[-5:5];
p=[-5:5];

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


x=[xe;xi];
y=[ye;yi];
z=[ze;zi];

vx=[vxe;vxi];
vy=[vye;vyi];
vz=[vze;vzi];

QM=[QMe;QMi];
Q=[Qe;Qi];

T=10;
dt=0.1;
t=[0:dt:T];

% fig=figure;
% set(fig,'color','k');
% mov = avifile('el_ion_density.avi','fps',25,'quality',100,'compression','Indeo5');

for it=1:length(t)
    %t(it);
r2=x.^2+y.^2+z.^2;
[r2,ord]=sort(r2);
Q=Q(ord);
j=1;
s=1;   
    %elettroni
    
    %xe=xe(ord);
    %ye=ye(ord);
    %ze=ze(ord);
    %vxe=vxe(ord);
    %vye=vye(ord);
    %vze=vze(ord);
   
r2e=(xe.^2+ye.^2+ze.^2);
r2e=sort(r2e);
re=sqrt(r2e);
re=sort(re);
    
    for i=1:length(Q)
        if Q(i)<0
            w(j)=i;
            j=j+1;
        else
            j=j;
        end
    end
    
    Qt=cumsum(Q);
    Qtn=Qt(w);
    Ee=e*Qtn./r2e;
    vxe=vxe+dt*QMe.*Ee.*xe./re;
    vye=vye+dt*QMe.*Ee.*ye./re;
    vze=vze+dt*QMe.*Ee.*ze./re;
    xe=xe+vxe*dt;
    ye=ye+vye*dt;
    ze=ze+vze*dt;
    
    %ioni
     
r2i=(xi.^2+yi.^2+zi.^2);
r2i=sort(r2i);
ri=sqrt(r2i);
ri=sort(ri);
    
     
    for k=1:length(Q)
        if Q(k)>0
            h(s)=k;
            s=s+1;
        else
            s=s;
        end
    end
    
    %Qt=cumsum(Q);
    Qsn=Qt(h);
    Ei=e*Qsn./r2i
    vxi=vxi+dt*QMi.*Ei.*xi./ri;
    vyi=vyi+dt*QMi.*Ei.*yi./ri;
    vzi=vzi+dt*QMi.*Ei.*zi./ri;
    xi=xi+vxi*dt;
    yi=yi+vyi*dt;
    zi=zi+vzi*dt;
    
    
    %rc=((0.5+(0:49))/50)*5;
    %d_rc=rc(2)-rc(1);
    %Vol=4*pi/3*((rc+d_rc/2).^3-(rc-d_rc/2).^3);
    
figure(1)
plot3(xe,ye,ze,'.r','MarkerSize',3)
axis([-5 5 -5 5 -5 5])
hold on
plot3(xi,yi,zi,'.b','MarkerSize',3)
axis([-5 5 -5 5 -5 5])
%hold on
%[k,l,p]=sphere;
%mesh ( k,l,p )
hold off
%axis('equal')

    %Nel=Ne;
    %Nio=Ni;
   
end


