clear all
close all
clc

colordef black

N=1e5;

R=1;

e=1;
QQ=1;
me=1;
m=50;
mi=abs(m*randn(N,1));

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
QMi=(e./mi);
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
    %
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
        figure(1)
        plot(rc,rho_e,'r','Linewidth',3)
        hold on
        plot(rc,rho_i,'m','Linewidth',3)
        legend('electrons','ions')
        title ('istante iniziale')
    end
        
    if it==(T/dt)
        figure(2)
        plot(rc,rho_e,'r','Linewidth',3)
        hold on
        plot(rc,rho_i,'m','Linewidth',3)
        legend('electrons','ions')
        title('istante finale')
        
    end
    figure(3)
    plot(rc,rho_e,'y','Linewidth',3)
    hold on
   plot(rc,rho_i,'m','Linewidth',3)
    hold off
    axis([0 3 0 .4])
    xlabel('Radius','fontsize',16)
    ylabel('Density','fontsize',16)
    legend('electrons','ions')
%     F=getframe(gcf);
%     mov=addframe(mov,F);
    drawnow
    Nel=Ne;
    Nio=Ni;
   
end





