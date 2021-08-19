%  Q7 Correction Matrix method - CATAM 23.4
 
X = 0.7;
Y = 1-X;
M = 1.9891*10^(30);
R = 6.9598*10^(8);
G = 6.6726*10^(-11);
L = 3.8515*10^(26);
P = 10^(15);
T = 10^(7);
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
                             
    % Defining functions for use in Euler method below              
HE = @(m,r) (G*m)/((4*pi*(r^4)));                                                                % Hydrostatic equilibirum equation 
MC = @(r) 1/(4*pi*(r^2));                                                                        % Mass continuity 
EG = @(T)   (((0.25*(X^2)*exp(-33.8*T^(-1/3)))+8.8*(10^(18))*X*exp(-152.28*T^(-1/3))))*T^(-2/3); % Energy generation
ET = @(r,t,l)  (3*k*l)/((((64*(pi^2)*a*c*(r^(4))*(t^3)))));                                      % Energy transport 
rho = @(p,T) ((p*mu)/(rc*T)); 

h=1/(10000); % Euler step-size 
d=1/(10000); % Delta change for parameters to calculate entries of 
             % the Jacobian numerically
Tc =2.89*10^(7);
Pc =7.15*10^(15);
RO=1.5011*R;
LM=L*155.4580; % Initial conditions (varied here from question 6 
               % to demonstate that this method is general)
dR=0;
dL=0;
dP=0;
dT=0; 
x1=1;
x2=1;
x3=1;
x4=1; % Values defined here for differences between variables 
      % at meeting point and corrections 
      % such that we can enter the while loop below

xvector=zeros(4,1);
dvector=zeros(4,1);
J=zeros(4,4); % Empty arrays defined to be filled during iteration below

        poutvector = zeros(1,1.5/h+1001);
        routvector = zeros(1,1.5/h+1001);
        moutvector = zeros(1,1.5/h+1001);
        loutvector = zeros(1,1.5/h+1001);
        toutvector = zeros(1,1.5/h+1001); 
        pinvector = zeros(1,1.5/h+999);
        rinvector = zeros(1,1.5/h+999);
        minvector = zeros(1,1.5/h+999);
        linvector = zeros(1,1.5/h+999);
        tinvector = zeros(1,1.5/h+999);    
        
        poutvectordR = zeros(1,1.5/h+1001);
        routvectordR = zeros(1,1.5/h+1001);
        moutvectordR = zeros(1,1.5/h+1001);
        loutvectordR = zeros(1,1.5/h+1001);
        toutvectordR = zeros(1,1.5/h+1001); 
        pinvectordR = zeros(1,1.5/h+999);
        rinvectordR = zeros(1,1.5/h+999);
        minvectordR = zeros(1,1.5/h+999);
        linvectordR = zeros(1,1.5/h+999);
        tinvectordR = zeros(1,1.5/h+999);   
        
        poutvectordL = zeros(1,1.5/h+1001);
        routvectordL = zeros(1,1.5/h+1001);
        moutvectordL = zeros(1,1.5/h+1001);
        loutvectordL = zeros(1,1.5/h+1001);
        toutvectordL = zeros(1,1.5/h+1001); 
        pinvectordL = zeros(1,1.5/h+999);
        rinvectordL = zeros(1,1.5/h+999);
        minvectordL = zeros(1,1.5/h+999);
        linvectordL = zeros(1,1.5/h+999);
        tinvectordL = zeros(1,1.5/h+999);   
        
        poutvectordP = zeros(1,1.5/h+1001);
        routvectordP = zeros(1,1.5/h+1001);
        moutvectordP = zeros(1,1.5/h+1001);
        loutvectordP = zeros(1,1.5/h+1001);
        toutvectordP = zeros(1,1.5/h+1001); 
        pinvectordP = zeros(1,1.5/h+999);
        rinvectordP = zeros(1,1.5/h+999);
        minvectordP = zeros(1,1.5/h+999);
        linvectordP = zeros(1,1.5/h+999);
        tinvectordP = zeros(1,1.5/h+999);   
        
        poutvectordT = zeros(1,1.5/h+1001);
        routvectordT = zeros(1,1.5/h+1001);
        moutvectordT = zeros(1,1.5/h+1001);
        loutvectordT = zeros(1,1.5/h+1001);
        toutvectordT = zeros(1,1.5/h+1001); 
        pinvectordT = zeros(1,1.5/h+999);
        rinvectordT = zeros(1,1.5/h+999);
        minvectordT = zeros(1,1.5/h+999);
        linvectordT = zeros(1,1.5/h+999);
        tinvectordT = zeros(1,1.5/h+999);    % All required arrays defined here
                                             % require a new set for each
                                             % row colum of the Jacobian as
                                             % we vary each parameter to
                                             % calculate each entry
                                             % numerically 
        
while abs(x4)>10^(-5)   % Condition for iterations until desired accuracy met
while abs(x3)>10^(-5)       
while abs(x2)>10^(-5)       
while abs(x1)>10^(-5) 

       RO=RO+dR; 
       Tc=Tc+dT;
       LM=LM+dL;
       Pc=Pc+dP;        % Correctio calculated on previous run of the 
                        % loop are applied 

rhoc = (Pc*mu)/(rc*Tc);     % New denisty defined 

% Beginnings of main loop  

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
end  % Normal loop
    x1=(routvector(1.5/h+1001)-rinvector(1))/RO;
    x2=(loutvector(1.5/h+1001)-linvector(1))/LM;
    x3=(poutvector(1.5/h+1001)-pinvector(1))/Pc;
    x4=(toutvector(1.5/h+1001)-tinvector(1))/Tc;    % Differences between the values at the meeting point defined here - we aim to minimise these values

xvector(1)=x1;
xvector(2)=x2;
xvector(3)=x3;
xvector(4)=x4;      % Put into vector for use in matrix equation below

% Beginning of loop where R* is varied  for column 1 of Jacobian
for n = 1:1:(1.5/h+1001)
    if n==1
        moutvectordR(1) = 0;
        poutvectordR(1) = Pc;
        routvectordR(1) = 0;
        loutvectordR(1) = 0;
        toutvectordR(1) = Tc;

    elseif n==2
     moutvectordR(n) = 0.001*h*M;
     routvectordR(n) = ((3*0.001*h*M)/(4*pi*rhoc))^(1/3); 
     poutvectordR(n) = Pc - (2/3)*(pi)*(G)*((rhoc)^2)*(routvectordR(n)^2);
     loutvectordR(n) = (4/3)*pi*(routvectordR(n)^3)*rhoc^2*EG(Tc/(10^6));
     toutvectordR(n) = ((( (Tc)^(4) - ((k*rhoc^(3)*EG(Tc/(10^6))*(routvectordR(n)^2)))/((2*a*c)))))^(1/4);
   elseif n>2 && n<=1002
        poutvectordR(n) = poutvectordR(n-1)-0.001*h*M*HE(moutvectordR(n-1),routvectordR(n-1));
        toutvectordR(n) = toutvectordR(n-1)-0.001*h*M*ET(routvectordR(n-1),toutvectordR(n-1), loutvectordR(n-1));
        routvectordR(n) = routvectordR(n-1)+0.001*h*M*MC(routvectordR(n-1))*(1/(rho(poutvectordR(n-1),toutvectordR(n-1))));
        moutvectordR(n) = moutvectordR(n-1)+0.001*h*M;
        loutvectordR(n) = loutvectordR(n-1)+0.001*h*M*EG(toutvectordR(n-1)/(10^6))*(rho(poutvectordR(n-1),toutvectordR(n-1)));
    else
        poutvectordR(n) = poutvectordR(n-1)-h*M*HE(moutvectordR(n-1),routvectordR(n-1));
        toutvectordR(n) = toutvectordR(n-1)-h*M*ET(routvectordR(n-1), toutvectordR(n-1), loutvectordR(n-1));
        routvectordR(n) = routvectordR(n-1)+h*M*MC(routvectordR(n-1))*(1/(rho(poutvectordR(n-1),toutvectordR(n-1))));
        moutvectordR(n) = moutvectordR(n-1)+h*M;
        loutvectordR(n) = loutvectordR(n-1)+h*M*EG(toutvectordR(n-1)/(10^6))*(rho(poutvectordR(n-1),toutvectordR(n-1)));
   end 
end 

for n=(1.5/h+999):-1:1
    if n == 1.5/h+999
        rinvectordR(1.5/h+999) = RO+d*RO;
        minvectordR(1.5/h+999) = 3*M;
        pinvectordR(1.5/h+999) = ((2/3)*(G*3*M))/(((rinvectordR(1.5/h+999))^2)*k);
        linvectordR(1.5/h+999) = LM;
        tinvectordR(1.5/h+999) = ((((linvectordR(1.5/h+999))/((4*pi*(rinvectordR(1.5/h+999)^2)*0.25*a*c)))))^(1/4);
    elseif n<(1.5/h)+999 && n>((1.5/h)-1)
        pinvectordR(n) = pinvectordR(n+1)+0.001*h*M*HE(minvectordR(n+1),rinvectordR(n+1));
        tinvectordR(n) = tinvectordR(n+1)+0.001*h*M*ET(rinvectordR(n+1), tinvectordR(n+1), linvectordR(n+1));
        rinvectordR(n) = rinvectordR(n+1)-0.001*h*M*MC(minvectordR(n+1))*(1/(rho(pinvectordR(n+1),tinvectordR(n+1))));
        minvectordR(n) = minvectordR(n+1)-0.001*h*M;
        linvectordR(n) = linvectordR(n+1)-0.001*h*M*EG(tinvectordR(n+1)/(10^6))*(rho(pinvectordR(n+1),tinvectordR(n+1)));
    else 
        pinvectordR(n) = pinvectordR(n+1)+h*M*HE(minvectordR(n+1),rinvectordR(n+1));
        tinvectordR(n) = tinvectordR(n+1)+h*M*ET(rinvectordR(n+1), tinvectordR(n+1), linvectordR(n+1));
        rinvectordR(n) = rinvectordR(n+1)-h*M*MC(rinvectordR(n+1))*(1/(rho(pinvectordR(n+1),tinvectordR(n+1))));
        minvectordR(n) = minvectordR(n+1)-h*M;
        linvectordR(n) = linvectordR(n+1)-h*M*EG(tinvectordR(n+1)/(10^6))*(rho(pinvectordR(n+1),tinvectordR(n+1)));
    end 
end % dR Loop
    x1dR=(routvectordR(1.5/h+1001)-rinvectordR(1))/RO;
    x2dR=(loutvectordR(1.5/h+1001)-linvectordR(1))/LM;
    x3dR=(poutvectordR(1.5/h+1001)-pinvectordR(1))/Pc;
    x4dR=(toutvectordR(1.5/h+1001)-tinvectordR(1))/Tc;

% Beginning of loop where L* is varied for column 2 of Jacobian
for n = 1:1:(1.5/h+1001)
    if n==1
        moutvectordL(1) = 0;
        poutvectordL(1) = Pc;
        routvectordL(1) = 0;
        loutvectordL(1) = 0;
        toutvectordL(1) = Tc;
    elseif n==2
     moutvectordL(n) = 0.001*h*M;
     routvectordL(n) = ((3*0.001*h*M)/(4*pi*rhoc))^(1/3); 
     poutvectordL(n) = Pc - (2/3)*(pi)*(G)*((rhoc)^2)*(routvectordL(n)^2);
     loutvectordL(n) = (4/3)*pi*(routvectordL(n)^3)*rhoc^2*EG(Tc/(10^6));
     toutvectordL(n) = ((( (Tc)^(4) - ((k*rhoc^(3)*EG(Tc/(10^6))*(routvectordL(n)^2)))/((2*a*c)))))^(1/4);
   elseif n>2 && n<=1002
        poutvectordL(n) = poutvectordL(n-1)-0.001*h*M*HE(moutvectordL(n-1),routvectordL(n-1));
        toutvectordL(n) = toutvectordL(n-1)-0.001*h*M*ET(routvectordL(n-1),toutvectordL(n-1), loutvectordL(n-1));
        routvectordL(n) = routvectordL(n-1)+0.001*h*M*MC(routvectordL(n-1))*(1/(rho(poutvectordL(n-1),toutvectordL(n-1))));
        moutvectordL(n) = moutvectordL(n-1)+0.001*h*M;
        loutvectordL(n) = loutvectordL(n-1)+0.001*h*M*EG(toutvectordL(n-1)/(10^6))*(rho(poutvectordL(n-1),toutvectordL(n-1)));
    else
        poutvectordL(n) = poutvectordL(n-1)-h*M*HE(moutvectordL(n-1),routvectordL(n-1));
        toutvectordL(n) = toutvectordL(n-1)-h*M*ET(routvectordL(n-1), toutvectordL(n-1), loutvectordL(n-1));
        routvectordL(n) = routvectordL(n-1)+h*M*MC(routvectordL(n-1))*(1/(rho(poutvectordL(n-1),toutvectordL(n-1))));
        moutvectordL(n) = moutvectordL(n-1)+h*M;
        loutvectordL(n) = loutvectordL(n-1)+h*M*EG(toutvectordL(n-1)/(10^6))*(rho(poutvectordL(n-1),toutvectordL(n-1)));
   end 
end 

for n=(1.5/h+999):-1:1
    if n == 1.5/h+999
        rinvectordL(1.5/h+999) = RO;
        minvectordL(1.5/h+999) = 3*M;
        pinvectordL(1.5/h+999) = ((2/3)*(G*3*M))/(((rinvectordL(1.5/h+999))^2)*k);
        linvectordL(1.5/h+999) = LM+d*LM;
        tinvectordL(1.5/h+999) = ((((linvectordL(1.5/h+999))/((4*pi*(rinvectordL(1.5/h+999)^2)*0.25*a*c)))))^(1/4);
    elseif n<(1.5/h)+999 && n>((1.5/h)-1)
        pinvectordL(n) = pinvectordL(n+1)+0.001*h*M*HE(minvectordL(n+1),rinvectordL(n+1));
        tinvectordL(n) = tinvectordL(n+1)+0.001*h*M*ET(rinvectordL(n+1), tinvectordL(n+1), linvectordL(n+1));
        rinvectordL(n) = rinvectordL(n+1)-0.001*h*M*MC(minvectordL(n+1))*(1/(rho(pinvectordL(n+1),tinvectordL(n+1))));
        minvectordL(n) = minvectordL(n+1)-0.001*h*M;
        linvectordL(n) = linvectordL(n+1)-0.001*h*M*EG(tinvectordL(n+1)/(10^6))*(rho(pinvectordL(n+1),tinvectordL(n+1)));
    else 
        pinvectordL(n) = pinvectordL(n+1)+h*M*HE(minvectordL(n+1),rinvectordL(n+1));
        tinvectordL(n) = tinvectordL(n+1)+h*M*ET(rinvectordL(n+1), tinvectordL(n+1), linvectordL(n+1));
        rinvectordL(n) = rinvectordL(n+1)-h*M*MC(rinvectordL(n+1))*(1/(rho(pinvectordR(n+1),tinvectordL(n+1))));
        minvectordL(n) = minvectordL(n+1)-h*M;
        linvectordL(n) = linvectordL(n+1)-h*M*EG(tinvectordL(n+1)/(10^6))*(rho(pinvectordL(n+1),tinvectordL(n+1)));
    end 
end  % dL Loop
    x1dL=(routvectordL(1.5/h+1001)-rinvectordL(1))/RO;
    x2dL=(loutvectordL(1.5/h+1001)-linvectordL(1))/LM;
    x3dL=(poutvectordL(1.5/h+1001)-pinvectordL(1))/Pc;
    x4dL=(toutvectordL(1.5/h+1001)-tinvectordL(1))/Tc;


% Beginning of loop where Pc is varied for column 3 of Jacobian
for n = 1:1:(1.5/h+1001)
Pcd =Pc*(1+d);
rhocd = (Pcd*mu)/(rc*Tc);    

    if n==1
        moutvectordP(1) = 0;
        poutvectordP(1) = Pcd;
        routvectordP(1) = 0;
        loutvectordP(1) = 0;
        toutvectordP(1) = Tc;
    elseif n==2
     moutvectordP(n) = 0.001*h*M;
     routvectordP(n) = ((3*0.001*h*M)/(4*pi*rhocd))^(1/3); 
     poutvectordP(n) = Pcd - (2/3)*(pi)*(G)*((rhocd)^2)*(routvectordP(n)^2);
     loutvectordP(n) = (4/3)*pi*(routvectordP(n)^3)*rhocd^2*EG(Tc/(10^6));
     toutvectordP(n) = ((( (Tc)^(4) - ((k*rhocd^(3)*EG(Tc/(10^6))*(routvectordP(n)^2)))/((2*a*c)))))^(1/4);
   elseif n>2 && n<=1002
        poutvectordP(n) = poutvectordP(n-1)-0.001*h*M*HE(moutvectordP(n-1),routvectordP(n-1));
        toutvectordP(n) = toutvectordP(n-1)-0.001*h*M*ET(routvectordP(n-1),toutvectordP(n-1), loutvectordP(n-1));
        routvectordP(n) = routvectordP(n-1)+0.001*h*M*MC(routvectordP(n-1))*(1/(rho(poutvectordP(n-1),toutvectordP(n-1))));
        moutvectordP(n) = moutvectordP(n-1)+0.001*h*M;
        loutvectordP(n) = loutvectordP(n-1)+0.001*h*M*EG(toutvectordP(n-1)/(10^6))*(rho(poutvectordP(n-1),toutvectordP(n-1)));
    else
        poutvectordP(n) = poutvectordP(n-1)-h*M*HE(moutvectordP(n-1),routvectordP(n-1));
        toutvectordP(n) = toutvectordP(n-1)-h*M*ET(routvectordP(n-1), toutvectordP(n-1), loutvectordP(n-1));
        routvectordP(n) = routvectordP(n-1)+h*M*MC(routvectordP(n-1))*(1/(rho(poutvectordP(n-1),toutvectordP(n-1))));
        moutvectordP(n) = moutvectordP(n-1)+h*M;
        loutvectordP(n) = loutvectordP(n-1)+h*M*EG(toutvectordP(n-1)/(10^6))*(rho(poutvectordL(n-1),toutvectordP(n-1)));
   end 
end 

for n=(1.5/h+999):-1:1
    if n == 1.5/h+999
        rinvectordP(1.5/h+999) = RO;
        minvectordP(1.5/h+999) = 3*M;
        pinvectordP(1.5/h+999) = ((2/3)*(G*3*M))/(((rinvectordP(1.5/h+999))^2)*k);
        linvectordP(1.5/h+999) = LM;
        tinvectordP(1.5/h+999) = ((((linvectordP(1.5/h+999))/((4*pi*(rinvectordP(1.5/h+999)^2)*0.25*a*c)))))^(1/4);
    elseif n<(1.5/h)+999 && n>((1.5/h)-1)
        pinvectordP(n) = pinvectordP(n+1)+0.001*h*M*HE(minvectordP(n+1),rinvectordP(n+1));
        tinvectordP(n) = tinvectordP(n+1)+0.001*h*M*ET(rinvectordP(n+1), tinvectordP(n+1), linvectordP(n+1));
        rinvectordP(n) = rinvectordP(n+1)-0.001*h*M*MC(minvectordP(n+1))*(1/(rho(pinvectordP(n+1),tinvectordP(n+1))));
        minvectordP(n) = minvectordP(n+1)-0.001*h*M;
        linvectordP(n) = linvectordP(n+1)-0.001*h*M*EG(tinvectordP(n+1)/(10^6))*(rho(pinvectordP(n+1),tinvectordP(n+1)));
    else 
        pinvectordP(n) = pinvectordP(n+1)+h*M*HE(minvectordP(n+1),rinvectordP(n+1));
        tinvectordP(n) = tinvectordP(n+1)+h*M*ET(rinvectordP(n+1), tinvectordP(n+1), linvectordP(n+1));
        rinvectordP(n) = rinvectordP(n+1)-h*M*MC(rinvectordP(n+1))*(1/(rho(pinvectordP(n+1),tinvectordP(n+1))));
        minvectordP(n) = minvectordP(n+1)-h*M;
        linvectordP(n) = linvectordP(n+1)-h*M*EG(tinvectordP(n+1)/(10^6))*(rho(pinvectordP(n+1),tinvectordP(n+1)));
    end 
end % dP Loop
    x1dP=(routvectordP(1.5/h+1001)-rinvectordP(1))/RO;
    x2dP=(loutvectordP(1.5/h+1001)-linvectordP(1))/LM;
    x3dP=(poutvectordP(1.5/h+1001)-pinvectordP(1))/Pc;
    x4dP=(toutvectordP(1.5/h+1001)-tinvectordP(1))/Tc;


% Beginning of loop where Tc is varied for column 4 of Jacobian
for n = 1:1:(1.5/h+1001)
Tcd =Tc*(1+d);
rhocD = (Pc*mu)/(rc*Tcd); 
    if n==1
        moutvectordT(1) = 0;
        poutvectordT(1) = Pc;
        routvectordT(1) = 0;
        loutvectordT(1) = 0;
        toutvectordT(1) = Tcd;
       
    elseif n==2
     moutvectordT(n) = 0.001*h*M;
     routvectordT(n) = ((3*0.001*h*M)/(4*pi*rhocD))^(1/3); 
     poutvectordT(n) = Pc - (2/3)*(pi)*(G)*((rhocD)^2)*(routvectordT(n)^2);
     loutvectordT(n) = (4/3)*pi*(routvectordT(n)^3)*rhocD^2*EG(Tcd/(10^6));
     toutvectordT(n) = ((( (Tcd)^(4) - ((k*rhocD^(3)*EG(Tcd/(10^6))*(routvectordT(n)^2)))/((2*a*c)))))^(1/4);
   elseif n>2 && n<=1002
        poutvectordT(n) = poutvectordP(n-1)-0.001*h*M*HE(moutvectordT(n-1),routvectordT(n-1));
        toutvectordT(n) = toutvectordT(n-1)-0.001*h*M*ET(routvectordT(n-1),toutvectordT(n-1), loutvectordT(n-1));
        routvectordT(n) = routvectordT(n-1)+0.001*h*M*MC(routvectordT(n-1))*(1/(rho(poutvectordT(n-1),toutvectordT(n-1))));
        moutvectordT(n) = moutvectordT(n-1)+0.001*h*M;
        loutvectordT(n) = loutvectordT(n-1)+0.001*h*M*EG(toutvectordT(n-1)/(10^6))*(rho(poutvectordT(n-1),toutvectordT(n-1)));
    else
        poutvectordT(n) = poutvectordT(n-1)-h*M*HE(moutvectordT(n-1),routvectordT(n-1));
        toutvectordT(n) = toutvectordT(n-1)-h*M*ET(routvectordT(n-1), toutvectordT(n-1), loutvectordT(n-1));
        routvectordT(n) = routvectordT(n-1)+h*M*MC(routvectordT(n-1))*(1/(rho(poutvectordT(n-1),toutvectordT(n-1))));
        moutvectordT(n) = moutvectordT(n-1)+h*M;
        loutvectordT(n) = loutvectordT(n-1)+h*M*EG(toutvectordT(n-1)/(10^6))*(rho(poutvectordT(n-1),toutvectordT(n-1)));
   end 
end 

for n=(1.5/h+999):-1:1
    if n == 1.5/h+999
        rinvectordT(1.5/h+999) =RO;
        minvectordT(1.5/h+999) = 3*M;
        pinvectordT(1.5/h+999) = ((2/3)*(G*3*M))/(((rinvectordT(1.5/h+999))^2)*k);
        linvectordT(1.5/h+999) = LM;
        tinvectordT(1.5/h+999) = ((((linvectordT(1.5/h+999))/((4*pi*(rinvectordT(1.5/h+999)^2)*0.25*a*c)))))^(1/4);
    elseif n<(1.5/h)+999 && n>((1.5/h)-1)
        pinvectordT(n) = pinvectordT(n+1)+0.001*h*M*HE(minvectordT(n+1),rinvectordT(n+1));
        tinvectordT(n) = tinvectordT(n+1)+0.001*h*M*ET(rinvectordT(n+1), tinvectordT(n+1), linvectordT(n+1));
        rinvectordT(n) = rinvectordT(n+1)-0.001*h*M*MC(minvectordT(n+1))*(1/(rho(pinvectordT(n+1),tinvectordT(n+1))));
        minvectordT(n) = minvectordT(n+1)-0.001*h*M;
        linvectordT(n) = linvectordT(n+1)-0.001*h*M*EG(tinvectordT(n+1)/(10^6))*(rho(pinvectordT(n+1),tinvectordT(n+1)));
    else 
        pinvectordT(n) = pinvectordT(n+1)+h*M*HE(minvectordT(n+1),rinvectordT(n+1));
        tinvectordT(n) = tinvectordT(n+1)+h*M*ET(rinvectordT(n+1), tinvectordT(n+1), linvectordT(n+1));
        rinvectordT(n) = rinvectordT(n+1)-h*M*MC(rinvectordT(n+1))*(1/(rho(pinvectordT(n+1),tinvectordT(n+1))));
        minvectordT(n) = minvectordT(n+1)-h*M;
        linvectordT(n) = linvectordT(n+1)-h*M*EG(tinvectordT(n+1)/(10^6))*(rho(pinvectordT(n+1),tinvectordT(n+1)));
    end 
end  % dT Loop
    x1dT=(routvectordT(1.5/h+1001)-rinvectordT(1))/RO;
    x2dT=(loutvectordT(1.5/h+1001)-linvectordT(1))/LM;
    x3dT=(poutvectordT(1.5/h+1001)-pinvectordT(1))/Pc;
    x4dT=(toutvectordT(1.5/h+1001)-tinvectordT(1))/Tc;

J(1,1) = (-(x1dR-x1)/(d*(RO/R)));
J(1,2) = (-(x1dL-x1)/(d*(LM/L)));
J(1,3) = (-(x1dP-x1)/(d*(Pc/P)));
J(1,4) = (-(x1dT-x1)/(d*(Tc/T)));

J(2,1) = (-(x2dR-x2)/(d*(RO/R)));
J(2,2) = (-(x2dL-x2)/(d*(LM/L)));
J(2,3) = (-(x2dP-x2)/(d*(Pc/P)));
J(2,4) = (-(x2dT-x2)/(d*(Tc/T)));

J(3,1) = (-(x3dR-x3)/(d*(RO/R)));
J(3,2) = (-(x3dL-x3)/(d*(LM/L)));
J(3,3) = (-(x3dP-x3)/(d*(Pc/P)));
J(3,4) = (-(x3dT-x3)/(d*(Tc/T)));

J(4,1) = (-(x4dR-x4)/(d*(RO/R)));
J(4,2) = (-(x4dL-x4)/(d*(LM/L)));
J(4,3) = (-(x4dP-x4)/(d*(Pc/P)));
J(4,4) = (-(x4dT-x4)/(d*(Tc/T)));       % Jacobian elements defined here 
                                        % and rescaled by typical order of 
                                        % magnitude for a stable inversion
    

dvector=J\xvector;      % Matrix equation solved here using 
                        % \ MATLAB function for 4x4 matrix inversion

dR=dvector(1)*R/10;
dL=dvector(2)*L/10;
dP=dvector(3)*P/10;
dT=dvector(4)*T/10;     % Corrections to the initial conditions hence 
                        % defined here to be implemented on next iteration
end
end
end 
end 

    % Differnt plots needed commented in/out below (useful to plot in
    % conjuction with loop to show the convergence of results)

% plot(moutvector/M, routvector/R)
% hold on
% plot(minvector/M,rinvector/R)

% plot(moutvector/M, toutvector/Tc)
% hold on
% plot(minvector/M,tinvector/Tc)

% plot(moutvector/M, loutvector/L)
% hold on
% plot(minvector/M,linvector/L)

% plot(moutvector/M, poutvector/Pc)
% hold on
% plot(minvector/M,pinvector/Pc)


% title('Mass-radius plot','interpreter','latex')
% xlabel('Mass ($M/M_{\odot}$)','interpreter','latex')
% ylabel('Radius ($R/R_{\odot}$)','interpreter','latex')
% legend('Shooting outward from centre', 'Shooting inward from surface')     






