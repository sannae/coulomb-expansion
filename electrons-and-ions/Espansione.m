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

v=0.9;

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

dt=0.1;
t=[0:dt:100];

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

%PARTE SULLE DENSITA'
%    rc=((0.5+(0:49))/50)*5;
%    d_rc=rc(2)-rc(1);
%    Vol=4*pi/3*((rc+d_rc/2).^3-(rc-d_rc/2).^3);
%    Ne=hist(r(Q<0),rc)/N0e;
%    Ni=hist(r(Q>0),rc)/N0i;
%    if it>1
%       Nel=(Ne+Nel)/2;
%       Nio=(Ni+Nio)/2;
 %   else
 %      Nel=Ne;
 %      Nio=Ni;
%   end
%    rho_e=Nel./Vol;
%    rho_i=Nio./Vol;
%    plot(rc,rho_e,'y','Linewidth',3)
%    hold on
%    plot(rc,rho_i,'m','Linewidth',3)
%    hold off
%    axis([0 3 0 .4])
%    xlabel('Radius','fontsize',16)
%    ylabel('Density','fontsize',16)
%     F=getframe(gcf);
%     mov=addframe(mov,F);
%    Nel=Ne;
%    Nio=Ni;

%GRAFICO 3D
plot3(x(Q<0),y(Q<0),z(Q<0),'.b')
axis([-5 5 -5 5 -5 5])
grid on
hold on
plot3(x(Q>0),y(Q>0),z(Q>0),'.r')
axis([-5 5 -5 5 -5 5])
grid on
%SFERA DI RIFERIMENTO
%hold on
%[k,l,p]=sphere;
%mesh ( k,l,p )
hold off
drawnow

end





