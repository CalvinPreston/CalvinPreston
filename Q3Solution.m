    % Defining constants for use throughout
X = 0.7;
Y = 1-X;
M = 1.9891*10^(30);
R = 6.9598*10^(8);
G = 6.6726*10^(-11);
L = 3.8515*10^(26);
rhoa = 1408.567712;
sfb = 5.6704*10^(-8);
mf = 6.1752*10^(11);
r = 8314.5;
nmu = 2*X + 0.75*Y;
mu = 1/nmu;
gm = 5/3;
k = 0.02*(1+X);         % Constants defined

    % Defining function
e = @(t) (((((0.25*(X^2)*exp(-33.8*t^(-1/3)))+8.8*10^(18)*X*exp(-152.28*t^(-1/3))))*t^(-2/3))*(10^(6)*t)^3)*rhoa*(0.25*6.22*10^(11)) - L;

% Use of simple zero finder to give route 
 
fzero(e,10)
