% Q4 - A shooting solution - CATAM 23.4
X = 0.7;
Y = 1-X;
M = 1.9891*10^(30);
R = 6.9598*10^(8);
G = 6.6726*10^(-11);
L = 3.8515*10^(26);
sfb = 5.6704*10^(-8);
mf = 6.1752*10^(11);
rc = 8314.5;
nmu = 2*X + 0.75*Y;
mu = 1/nmu;
gm = 5/3;
k = 0.02*(1+X);         % Constants defined

% Two equations to solve now by Euler method as below
% Guess of central T and P values to calculate constant from ideal gas,
% employing the brotropic relationship between P and rho

Tc = 1.539*10^(7);
Pc = 1.619*10^(15);

K = (((Tc)^(gm))*((rc/mu)^(gm-1)))/(Pc^(gm-1));
Rhoc = ((Pc*mu)/(rc*K))^(1/gm);

h=1/(10000); % Step wise increse in mass
        poutvector = zeros(1,1.5/h+999);
        routvector = zeros(1,1.5/h+999);
        moutvector = zeros(1,1.5/h+999);    % Defining vectors here

HE = @(m,r,p) (G*m)/((4*pi*(r^4)));                          %  Hydrostatic equilibirum equation 
MC = @(m,r,p) 1/(4*pi*r^2*(((mu*p)/(rc*K))^(1/gm)));         %  Mass continuity equation

for n = 1:1:(1.5/h+999)
    if n==1
        moutvector(1) = 0;
        poutvector(1) = Pc;
        routvector(1) = 0;   
    elseif n==2
     routvector(n) = ((3*0.001*h*M)/(4*pi*Rhoc))^(1/3); 
     poutvector(n) = Pc - (2/3)*(pi)*G*((Rhoc)^2)*(routvector(n)^2);
     moutvector(n) = 0.001*h*M;
   elseif n>2 && n<=1000
        poutvector(n) = poutvector(n-1)-0.001*h*M*HE(moutvector(n-1),routvector(n-1), poutvector(n-1));
        routvector(n) = routvector(n-1)+0.001*h*M*MC(moutvector(n-1),routvector(n-1), poutvector(n-1));
        moutvector(n) = moutvector(n-1)+0.001*h*M;
   else 
        poutvector(n) = poutvector(n-1)-h*M*HE(moutvector(n-1),routvector(n-1), poutvector(n-1));
        routvector(n) = routvector(n-1)+h*M*MC(moutvector(n-1),routvector(n-1), poutvector(n-1));
        moutvector(n) = moutvector(n-1)+h*M;
   end 
end 
        pinvector = zeros(1,1.5/h+999);
        rinvector = zeros(1,1.5/h+999);
        minvector = zeros(1,1.5/h+999);
        
for n=(1.5/h+999):-1:1
    if n == 1.5/h+999
        rinvector(1.5/h+999) = 1.5*R;
        minvector(1.5/h+999) = 3*M;
        pinvector(1.5/h+999) = ((2/3)*(G*3*M))/(((1.5*R)^2)*k);
    elseif n<(1.5/h)+999 && n>((1.5/h)-1)
        pinvector(n) = pinvector(n+1)+0.001*h*M*HE(minvector(n+1),rinvector(n+1), pinvector(n+1));
        rinvector(n) = rinvector(n+1)-0.001*h*M*MC(minvector(n+1),rinvector(n+1), pinvector(n+1));
        minvector(n) = minvector(n+1)-0.001*h*M;
    else 
        pinvector(n) = pinvector(n+1)+h*M*HE(minvector(n+1),rinvector(n+1), pinvector(n+1));
        rinvector(n) = rinvector(n+1)-h*M*MC(minvector(n+1),rinvector(n+1), pinvector(n+1));
        minvector(n) = minvector(n+1)-h*M;
    end 
end 

plot(moutvector./M,poutvector./Pc);
hold on
plot(minvector./M,pinvector./Pc)

% plot(moutvector./M,routvector/R);
% hold on
% plot(minvector./M,rinvector/R)

title('Mass-pressure plot','interpreter','latex')
xlabel('Mass ($M/M_{\odot}$)','interpreter','latex')
ylabel('Presssure ($P/P_c$)', 'interpreter','latex')
legend('Shooting outward from centre', 'Shooting inward from surface')
