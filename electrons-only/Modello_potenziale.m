clear all;
close all;
N0 = 10000;
N = round(N0*8/(4/3*pi));
R = 1;
v = 1.2665;
e = 1;
m = 1;
Ne = 1;
x = -R + 2*R*rand(N,1);
y = -R + 2*R*rand(N,1);
z = -R + 2*R*rand(N,1);
vx = v.*randn(N,1);
vy = v.*randn(N,1);
vz = v.*randn(N,1);
dentro = x.^2 + y.^2 + z.^2 <= R^2;
x = x(dentro);
y = y(dentro);
z = z(dentro);
vx = vx(dentro);
vy = vy(dentro);
vz = vz(dentro);
passi = 1000;
tB = linspace(0,2*R/v,passi+1);
tB = tB(2:end);
t = 5*tB;
N0 = sum(dentro);
q = e*Ne/N0;
dt = t(end)/passi;
w = e/m*dt;
for ii = 1 : passi
r2 = x.^2 + y.^2 + z.^2;
[r2,ord] = sort(r2);
x = x(ord); 
y = y(ord); 
z = z(ord);
vx = vx(ord); 
vy = vy(ord); 
vz = vz(ord);
r = sqrt(r2);
E = (e*Ne*min(R,r).^3/R^3-q*(0:(N0-1)+.5)')./r2;
vx = vx-w*E.*x./r;
vy = vy-w*E.*y./r;
vz = vz-w*E.*z./r;
x = x + vx*dt;
y = y + vy*dt;
z = z + vz*dt;
dentro = x.^2 + y.^2 + z.^2 <= R^2;
rimasti(ii) = sum(dentro); 
%Blocco Potenziale%
if mod(ii,10)==1
    phi(N0)=0;
    for jj = N0:-1:2
        phi(jj-1) = phi(jj) + .5*(E(jj) + E(jj-1))*(r(jj)-r(jj-1)); %Incremento del potenziale
        end
Energia = (vx.^2 + vy.^2 + vz.^2)/2-phi';
ET = Energia + phi(jj-1);
sum(ET < 0);
plot(r,-phi,r,Energia,'.r','MarkerSize',1)
xlabel('t(s)')
ylabel('Phi(J)')
axis ([0 30 -2 10])
drawnow
end
end
