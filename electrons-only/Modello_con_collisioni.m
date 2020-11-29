clear all;
close all;
N0 =5000;
N = round(N0*8/(4/3*pi));
R = 1;
v = 0.2665;
e = 1;
m = 1;
Ne = 1;
epsb = 0.003;
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
E = (e*Ne*min(R,r).^3/R^3-q*(0:(N0-1))')./r2;
vx = vx-w*E.*x./r;
vy = vy-w*E.*y./r;
vz = vz-w*E.*z./r;
%Intermezzo sulle Collisioni
ib = (rand(N0,1) < epsb); 
    %Genera numeri a caso: ib contiene il True o False dell'espressione
vmod = sqrt(vx(ib).^2 + vy(ib).^2 + vz(ib).^2);
    %Di ogni vx,vy,vz prende solo i numeri che soddisfano l'espressione e
    %ne compone la vmod (velocità post-collisione)
nb = length(vmod);
    %Dimensione di vmod
cb = -1 + 2*rand(nb,1); 
phib = 2*pi*rand(nb,1); 
sb = sqrt(1 - cb.^2);
vx(ib) = vmod.*sb.*cos(phib);
vy(ib) = vmod.*sb.*sin(phib);
vz(ib) = vmod.*cb;
%Fine intermezzo
x = x + vx*dt;
y = y + vy*dt;
z = z + vz*dt;
dentro = x.^2 + y.^2 + z.^2 <= R^2;
rimasti(ii) = sum(dentro); 
end

figure(3)
plot(t,rimasti/N0*100,'b')
ylim([50 100]);