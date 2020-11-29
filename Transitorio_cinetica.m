clear all
close all
clc

colordef white

N=1e5;

R=1;
v=sqrt(7.2e-3);

x=-1+2*rand(N,1);
y=-1+2*rand(N,1);
z=-1+2*rand(N,1);

vx=v*randn(N,1);
vy=v*randn(N,1);
vz=v*randn(N,1);

dentro=x.^2+y.^2+z.^2<=R^2;
N0=sum(dentro);

x=x(dentro);
y=y(dentro);
z=z(dentro);

vx=vx(dentro);
vy=vy(dentro);
vz=vz(dentro);

dt=0.1;
t=[0:dt:10];

% fig=figure;
% set(fig,'color','k');
% mov = avifile('en_distr.avi','fps',25,'quality',100,'compression','Indeo5');

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
    r=sqrt(r2);
    N=((0:1:N0-1)+0.5)';
    E=(min(1,r).^3-N/N0)./r2;
    vx=vx-dt*E.*x./r;
    vy=vy-dt*E.*y./r;
    vz=vz-dt*E.*z./r;
    x=x+vx*dt;
    y=y+vy*dt;
    z=z+vz*dt;
    phi(N0)=0;
    for j=N0-1:-1:1
        phi(j)=phi(j+1)+0.5*(E(j)+E(j+1))*(r(j+1)-r(j));
    end
        
        Energia=(vx.^2+vy.^2+vz.^2)/2-phi';
        plot(r,Energia,'.r',r,-phi,'Linewidth',1.5,'MarkerSize',1)
        xlabel('x/R','fontsize',16)
        grid on
        axis([0 5 -0.5 .1])
        drawnow
%         F=getframe(gcf);
%         mov=addframe(mov,F);

         interno=x.^2+y.^2+z.^2<=R^2;
         int(it)=sum(interno);
         
         K(it)=sum((vx.^2+vy.^2+vz.^2)*.5);
end

figure (2)
plot(r,phi,'Linewidth',3)
grid on
xlabel('t','fontsize',16)
ylabel('Potential','fontsize',16)

figure(3)
plot(t,K,'r')


