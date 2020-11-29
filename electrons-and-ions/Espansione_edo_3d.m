clear all
close all
clc

colordef black

N=1e3;

R=1;

e=1;
QQ=1;
me=1;
mi=50;

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
    t(it);
    
    %elettroni
    
    r2e=xe.^2+ye.^2+ze.^2;
    [r2e,ord]=sort(r2e);
    xe=xe(ord);
    ye=ye(ord);
    ze=ze(ord);
    vxe=vxe(ord);
    vye=vye(ord);
    vze=vze(ord);
    QMe=QMe(ord);
    Qe=Qe(ord);
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
    [r2i,ord]=sort(r2i);
    xi=xi(ord);
    yi=yi(ord);
    zi=zi(ord);
    vxi=vxi(ord);
    vyi=vyi(ord);
    vzi=vzi(ord);
    QMi=QMi(ord);
    Qi=Qi(ord);
    ri=sqrt(r2i);
    Ei=e*cumsum(Qi)./r2i;
    vxi=vxi+dt*QMi.*Ei.*xi./ri;
    vyi=vyi+dt*QMi.*Ei.*yi./ri;
    vzi=vzi+dt*QMi.*Ei.*zi./ri;
    xi=xi+vxi*dt;
    yi=yi+vyi*dt;
    zi=zi+vzi*dt;
    
    
    rc=((0.5+(0:49))/50)*5;
    d_rc=rc(2)-rc(1);
    Vol=4*pi/3*((rc+d_rc/2).^3-(rc-d_rc/2).^3);
    Ne=hist(r(Q<0),rc)/N0e;
    Ni=hist(r(Q>0),rc)/N0i;
    if it>1
       Nel=(Ne+Nel)/2;
       Nio=(Ni+Nio)/2;
    else
       Nel=Ne;
       Nio=Ni;
      
    end
    rho_e=Nel./Vol;
    rho_i=Nio./Vol;
    if it==1
        %figure(1)
        %plot(rc,rho_e,'r','Linewidth',3)
        %hold on
        %plot(rc,rho_i,'m','Linewidth',3)
        %legend('electrons','ions')
        %title ('istante iniziale')
    end
        
    if it==(T/dt)
        %figure(2)
        %plot(rc,rho_e,'r','Linewidth',3)
        %hold on
        %plot(rc,rho_i,'m','Linewidth',3)
        %legend('electrons','ions')
        %title('istante finale')
        
    end
    %figure(3)
    %plot(rc,rho_e,'y','Linewidth',3)
    %hold on
   %plot(rc,rho_i,'m','Linewidth',3)
    %hold off
    %axis([0 3 0 .4])
    %xlabel('Radius','fontsize',16)
    %ylabel('Density','fontsize',16)
    %legend('electrons','ions')
%     F=getframe(gcf);
%     mov=addframe(mov,F);
figure
plot3(xe,ye,ze,'.r','MarkerSize',1)
hold on
plot3(xi,yi,zi,'.b','MarkerSize',1)
    drawnow
    Nel=Ne;
    Nio=Ni;
   
end





