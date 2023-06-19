function [critShields] = zanke03(D,rhof,rhos,temp,g,phi,pores,depth,damp)
%Zanke 2003 model for the Shields curve: [Shcr,D] = zanke03(optional input);
%Zanke, U.C.E. 2003. On the influence of turbulence on the initiation of 
% sediment motion, Int. J. of Sediment Research 18(1), 1-15.
%Function for computation of critical Shields parameter for input parameters:
% D diameter D (m), rhos sediment and rhof flow density (2650 or 3400 for silica
% or basalt, 1000 kg/m^3 for water), flow kinematic viscosity
% (1.6e-6 for nearly freezing water), gravity (9.81 for Earth, 3.74 for Mars m/s^2)
% phi angle of repose of sediment (30° sand 45° angular gravel), pores fraction 0.3
% depth is flow depth (m) which affects turbulence, damp is factor for turbulence
% damping (<1 is damped)

%TRY DIFFERENT ROUGHNESS PREDICTORS!

%input
if nargin < 2
   rhof = 1000;
   rhos = 2650;
   temp = 15; %°C
   g = 9.81; %3.74
   phi = 35;
   pores = 0.35;
   depth = 100;
   damp = 1;
end;
if nargin == 0
   psiD = [-8:0.25:9]';
   Dmm = 2.^psiD;
   D = Dmm./1000;
end;

%constants
visco = 4e-5/(20+temp); %kinematic viscosity of water (1.2e-6 for cold water)
R = (rhos-rhof)/rhof;		  %relative submerged density of sediment
Ddim = D.*(R*g/visco^2).^(1/3); %Bonnefille dimensionless grain size
phir = 2*pi*phi/360; %deg2rad

%approximation function for first step by Kleinhans (for large water depths and D<5mm)
Shini=(0.145.*Ddim.^(-1/2) + 0.045*10.^(-1100.*Ddim.^(-9/4)));
%inilarge = find(D>0.005);
%if length(inilarge)>=1
%   Shini(inilarge) = Shini(inilarge(1)-1);
%end;
%critShields=Shini;

%Empirical Shields curves
Rep = D.^(1.5).*sqrt(R*g)/visco;
Rstar = Rep.^(-0.6);
Brownlie = 0.22.*Rstar+0.06.*10.^(-7.7.*Rstar);
Soulsby = 0.3./(1+1.2.*Ddim)+0.055.*(1-exp(-0.02.*Ddim)); 
%loglog(D,Shini,D,Brownlie,D,Soulsby);

%Zanke's equations
Kohesion = 1+3e-8./((rhos-rhof).*D.^2); % eq. 37
n = damp*1.8; %turbulence damping factor, normal stand dev for turb=1.8
ks = 1.*D;

minimise = repmat(1,size(D));
tel = 0;
while max(minimise)>1e-8%loop starts
tel=tel+1;
Tauini = Shini.*((rhos-rhof).*g.*D);
ustarini = sqrt(Tauini./rhof);
ksplusini = ustarini.*ks./visco;
Pt = 1-exp(-0.08.*ksplusini); %probability that flow is turbulent
Bsand = (1-Pt).*(2.5.*log(ksplusini)+5.25)+8.5.*Pt; %eq. 31
Bnatrough = 2.5.*log(1./(0.33+0.11./ksplusini)); %eq. 31a
%eq. 28, water depth effect and Driest damping:
uprmsstarini = 0.31.*ksplusini.*exp(-0.1.*ksplusini) +...
   1.8.*exp(-0.88.*D./depth).*(1-exp(-0.1.*ksplusini));
uystarini = (((1-Pt)./ksplusini.^2)+Pt./Bsand.^2).^-0.5; %eq. 29
%uystarini = (((1-Pt)./ksplusini.^2)+Pt./Bnatrough.^2).^-0.5; %eq. 29
ubstarini = 0.8+0.9.*uystarini; %eq. 32
upuini = uprmsstarini./ubstarini;
%eq. 15c the final (initial) answer:
Shsecond = ((1-pores)*tan(phir/1.5).*Kohesion)./... %denominator
   ( ((1+n.*upuini).^2).*... %part one of numerator
   (1+(1/2.5).*((n.*uprmsstarini).^2).*tan(phir/1.5).*Kohesion) );
minimise = abs(Shsecond - Shini);
Shini = Shsecond;
end; %loop ends
critShields=Shini;