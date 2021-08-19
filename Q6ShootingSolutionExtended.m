% Q5 - A shooting solution: Extended - CATAM 23.4

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
k = 0.02*(1+X);         
a = 7.5646*10^(-16);
c = 2.9979*10^(8);           % Constants defined (to accuracy 
                             % given in project brief)

Tc = 2.895*10^(7);
Pc = 7.210*10^(15);
RO = 1.536*R;
LM = 157.1*L;
rhoc = (Pc*mu)/(rc*Tc);          % Initial conditions defined

    % Defining functions for use in Euler method below
HE = @(m,r) (G*m)/((4*pi*(r^4)));                                                                % Hydrostatic equilibirum equation 
MC = @(r) 1/(4*pi*(r^2));                                                                        % Mass continuity 
EG = @(T)   (((0.25*(X^2)*exp(-33.8*T^(-1/3)))+8.8*(10^(18))*X*exp(-152.28*T^(-1/3))))*T^(-2/3); % Energy generation
ET = @(r,t,l)  (3*k*l)/((((64*(pi^2)*a*c*(r^(4))*(t^3)))));                                      % Energy transport 
rho = @(p,T) ((p*mu)/(rc*T));

h=1/(10000);        % Euler stepsize

        poutvector = zeros(1,1.5/h+1001);
        routvector = zeros(1,1.5/h+1001);
        moutvector = zeros(1,1.5/h+1001);
        loutvector = zeros(1,1.5/h+1001);
        toutvector = zeros(1,1.5/h+1001);        % Vectors defined for use
                                                 % below
                                                 
    % Shooting outward from centre to surface
for n = 1:1:(1.5/h+1001)
    if n==1
        moutvector(1) = 0;
        poutvector(1) = Pc;
        routvector(1) = 0;
        loutvector(1) = 0;
        toutvector(1) = Tc;
    elseif n==2
     moutvector(n) = 0.001*h*M;
     routvector(n) = ((3*0.001*h*M)/(4*pi*rhoc))^(1/3); 
     poutvector(n) = Pc - (2/3)*(pi)*(G)*((rhoc)^2)*(routvector(n)^2);
     loutvector(n) = (4/3)*pi*(routvector(n)^3)*rhoc^2*EG(Tc/(10^6));
     toutvector(n) = ((( (Tc)^(4) - ((k*rhoc^(3)*EG(Tc/(10^6))*(routvector(n)^2)))/((2*a*c)))))^(1/4);
   elseif n>2 && n<=1002
        poutvector(n) = poutvector(n-1)-0.001*h*M*HE(moutvector(n-1),routvector(n-1));
        toutvector(n) = toutvector(n-1)-0.001*h*M*ET(routvector(n-1),toutvector(n-1), loutvector(n-1));
        routvector(n) = routvector(n-1)+0.001*h*M*MC(routvector(n-1))*(1/(rho(poutvector(n-1),toutvector(n-1))));
        moutvector(n) = moutvector(n-1)+0.001*h*M;
        loutvector(n) = loutvector(n-1)+0.001*h*M*EG(toutvector(n-1)/(10^6))*(rho(poutvector(n-1),toutvector(n-1)));
    else
        poutvector(n) = poutvector(n-1)-h*M*HE(moutvector(n-1),routvector(n-1));
        toutvector(n) = toutvector(n-1)-h*M*ET(routvector(n-1), toutvector(n-1), loutvector(n-1));
        routvector(n) = routvector(n-1)+h*M*MC(routvector(n-1))*(1/(rho(poutvector(n-1),toutvector(n-1))));
        moutvector(n) = moutvector(n-1)+h*M;
        loutvector(n) = loutvector(n-1)+h*M*EG(toutvector(n-1)/(10^6))*(rho(poutvector(n-1),toutvector(n-1)));
   end 
end 
        pinvector = zeros(1,1.5/h+999);
        rinvector = zeros(1,1.5/h+999);
        minvector = zeros(1,1.5/h+999);
        linvector = zeros(1,1.5/h+999);
        tinvector = zeros(1,1.5/h+999); % Vectors defined for use
                                        % below
                                                
        % Shooting inward from surface to meeting point
for n=(1.5/h+999):-1:1
    if n == 1.5/h+999
        rinvector(1.5/h+999) = RO;
        minvector(1.5/h+999) = 3*M;
        pinvector(1.5/h+999) = ((2/3)*(G*3*M))/(((rinvector(1.5/h+999))^2)*k);
        linvector(1.5/h+999) = LM;
        tinvector(1.5/h+999) = ((((linvector(1.5/h+999))/((4*pi*(rinvector(1.5/h+999)^2)*0.25*a*c)))))^(1/4);
    elseif n<(1.5/h)+999 && n>((1.5/h)-1)
        pinvector(n) = pinvector(n+1)+0.001*h*M*HE(minvector(n+1),rinvector(n+1));
        tinvector(n) = tinvector(n+1)+0.001*h*M*ET(rinvector(n+1), tinvector(n+1), linvector(n+1));
        rinvector(n) = rinvector(n+1)-0.001*h*M*MC(minvector(n+1))*(1/(rho(pinvector(n+1),tinvector(n+1))));
        minvector(n) = minvector(n+1)-0.001*h*M;
        linvector(n) = linvector(n+1)-0.001*h*M*EG(tinvector(n+1)/(10^6))*(rho(pinvector(n+1),tinvector(n+1)));
    else 
        pinvector(n) = pinvector(n+1)+h*M*HE(minvector(n+1),rinvector(n+1));
        tinvector(n) = tinvector(n+1)+h*M*ET(rinvector(n+1), tinvector(n+1), linvector(n+1));
        rinvector(n) = rinvector(n+1)-h*M*MC(rinvector(n+1))*(1/(rho(pinvector(n+1),tinvector(n+1))));
        minvector(n) = minvector(n+1)-h*M;
        linvector(n) = linvector(n+1)-h*M*EG(tinvector(n+1)/(10^6))*(rho(pinvector(n+1),tinvector(n+1)));
    end 
end 

% plot(moutvector/M, loutvector/L)
% hold on
% plot(minvector/M,linvector/L)
% 
% plot(moutvector/M, routvector/R)
% hold on
% plot(minvector/M,rinvector/R)
% 
% plot(moutvector/M, toutvector/Tc)
% hold on
% plot(minvector/M,tinvector/Tc)

plot(moutvector/M, poutvector/Pc)  
hold on
plot(minvector/M,pinvector/Pc)

title('Mass-pressure plot','interpreter','latex')
xlabel('Mass ($M/M_{\odot}$)','interpreter','latex')
ylabel('Pressure ($P/P_c$)','interpreter','latex')
legend('Shooting outward from centre', 'Shooting inward from surface')
% 
% title('Mass-radius plot','interpreter','latex')
% xlabel('Mass ($M/M_{\odot}$)','interpreter','latex')
% ylabel('Radius ($R/R_{\odot}$)','interpreter','latex')
% legend('Shooting outward from centre', 'Shooting inward from surface')
