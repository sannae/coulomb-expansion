clear all
close all
clc

colordef white

N=1e4;

R=1;

rand('twister',5489)
x=-1+2*rand(N,1);
y=-1+2*rand(N,1);
z=-1+2*rand(N,1);

dentro=x.^2+y.^2+z.^2<=R^2;
N0=sum(dentro);

Q=1;
M=1;

q=Q/N0;
m=M/N0;

rho_i=Q/(4/3*pi*R^3); %densità di carica

x=x(dentro);
y=y(dentro);
z=z(dentro);

Etot=11.187e-5/N0; %???
v=sqrt(2*Etot/(3*m));

randn('state',0)
vx=v*randn(N0,1);
randn('state',2)
vy=v*randn(N0,1);
randn('state',4)
vz=v*randn(N0,1);

t_end=2*R/v;
passi=50;
t=linspace(0,t_end,passi+1);
t=2*t;
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
    
    K(it)=sum(0.5*m*(vx.^2+vy.^2+vz.^2));
    
    r=sqrt(r2);
    
    N=((0:1:N0-1)+0.5)';
    E=q*(min(1,r).^3-N/N0)./r2;
    
    vx=vx-dt*q/m*E.*x./r;
    vy=vy-dt*q/m*E.*y./r;
    vz=vz-dt*q/m*E.*z./r;
    
    x=x+vx*dt;
    y=y+vy*dt;
    z=z+vz*dt;
    
    phi(N0)=0;
    for j=N0-1:-1:1
        phi(j)=phi(j+1)+0.5*(E(j)+E(j+1))*(r(j+1)-r(j));
    end
    
    int=4*pi*r.^2.*phi';
    int=(int(r<=R));
    rr=r(r<=R);
    U_i=rho_i/2*trapz(rr,int); %???
    U_e=-1/2*q*sum(phi);        %???
    U(it)=U_i+U_e;
    
    eps(it)=U(it)+K(it);
end

Etot=sum(eps(1:end))/length(eps(1:end));

figure
plot(t,U,'Linewidth',3)
hold on
plot(t,K,'r','Linewidth',3)
plot(t,eps,'g','Linewidth',3)
text(t(end-15),eps(end-15),num2str(Etot))
grid on
xlabel('\omega t','fontsize',16)
ylabel('\epsilon/(e^2N_0/R^2)','fontsize',16)
legend('U','K','TOT')

