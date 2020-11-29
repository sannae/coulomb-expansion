close all
clear all
clc

colordef white

R=1;
Q=1;
M=1;
Mi=1836*M;
rho=500;

NN=(6/pi*4/3*pi*R^3*rho);

v=0.5;

% elettroni
x=-1+2*rand(NN,1);
y=-1+2*rand(NN,1);
z=-1+2*rand(NN,1);

dentroe=x.^2+y.^2+z.^2<=R^2;
N0e=sum(dentroe);

q=-Q/N0e*ones(N0e,1);
qm=-Q/M*ones(N0e,1);
ma=M/N0e*ones(N0e,1);

x=x(dentroe);
y=y(dentroe);
z=z(dentroe);

vx=v*randn(N0e,1);
vy=v*randn(N0e,1);
vz=v*randn(N0e,1);

% ioni
xi=-1+2*rand(NN,1);
yi=-1+2*rand(NN,1);
zi=-1+2*rand(NN,1);

dentroi=xi.^2+yi.^2+zi.^2<=R^2;
N0i=sum(dentroi);

qi=Q/N0i*ones(N0i,1);
qmi=Q/Mi*ones(N0i,1);
mai=Mi/N0i*ones(N0i,1);

xi=xi(dentroi);
yi=yi(dentroi);
zi=zi(dentroi);

vxi=v*zeros(N0i,1);
vyi=v*zeros(N0i,1);
vzi=v*zeros(N0i,1);

x=[x;xi];
y=[y;yi];
z=[z;zi];

vx=[vx;vxi];
vy=[vy;vyi];
vz=[vz;vzi];

Q=[q;qi];
QM=[qm;qmi];
MA=[ma;mai];

Nt=N0e+N0i;

t_end=100*R/v;
passi=5000;
t=linspace(0,t_end,passi+1);
dt=t(2)-t(1);

for it=1:length(t)
    it;
    
    r2=x.^2+y.^2+z.^2;
    [r2,ord]=sort(r2);
    
    x=x(ord);
    y=y(ord);
    z=z(ord);
    
    vx=vx(ord);
    vy=vy(ord);
    vz=vz(ord);
    
    Q=Q(ord);
    QM=QM(ord);
    MA=MA(ord);
    
    r=sqrt(r2);
    E=cumsum(Q)./r2;
    
    phi(Nt)=0;
    for j=Nt-1:-1:1
        phi(j)=phi(j+1)+0.5*(E(j)+E(j+1))*(r(j+1)-r(j));
    end
    
    rmax(it)=max(r(Q>0)); % distanza max a ogni iterazione di una particella positiva
    
    K=1/2*MA.*(vx.^2+vy.^2+vz.^2);
    Ktot(it)=sum(K);
    U=Q.*phi';
    Utot(it)=sum(U)/2;
    eps=K+U;
    epstot(it)=Ktot(it)+Utot(it);
    
    vx=vx+dt*QM.*E.*x./r;
    vy=vy+dt*QM.*E.*y./r;
    vz=vz+dt*QM.*E.*z./r;
    
    x=x+vx*dt;
    y=y+vy*dt;
    z=z+vz*dt;
end

figure
plot(t,rmax,'Linewidth',1.5)
xlabel('t\omega_p','fontsize',14)
ylabel('r/R','fontsize',14)

 figure
 plot(t,epstot)
 xlabel('t\omega_p','fontsize',14)
 ylabel('r/R','fontsize',14)
 hold on
 plot(t,Ktot,'r')
 hold on
 plot(t,Utot,'g')
 legend('epstot','Ktot','Utot')
 hold off

figure
Ki=K(Q>0);
Kim=mean(Ki);
Kis=std(Ki);
nc=20;
vett=zeros(1,nc);
DKi=(4*Kis+4*Kis)/nc;
KK=[Kim-4*Kis:DKi:Kim+4*Kis];
KK=KK+DKi/2;
KK=KK(1:nc);
Ki=Ki(Ki>=Kim-4*Kis);
Ki=Ki(Ki<=Kim+4*Kis);
for ik=1:length(Ki)
    cell=ceil((Ki(ik)-(Kim-4*Kis))/DKi);
    vett(cell)=vett(cell)+1;
end
plot(KK,vett/N0i*100,'Linewidth',1.5)
xlabel('K/(e^2N_0/R)','fontsize',14)
ylabel('N/N_0[%]','fontsize',14)
