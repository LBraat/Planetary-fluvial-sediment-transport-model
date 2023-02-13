%% Flow on Mars 2 - Initiation of Motion
% Author: Lisanne Braat (lisannebraat@gmail.com)
% Last update: 2023-02-13
% Created in Matlab version: R2022b

%% Initialize
clear variables
close all
clc

output = 'FlowOnMars_exportfig';
addpath(genpath('Checkout'));

%% Input parameters
for a = 1
    Drough = 0.01; %[m] Grain size, only used for roughness
    h = 2:0.5:15; %[m] Water depth
    hselect = 7; %[1 27];
    Q = 500:500:15000; %[m3/s] Discharge
    Qselect = 4; %[1 30];
    W = 200; %[m] Channel width
    %grav = 1:0.5:12; %[m/s2] Gravitational acceleration
    grav = [3.7 9.8]; %[m/s2] Gravitational acceleration
    S = 0.001; %[m/m] Slope
    TC = 4; %[degrees C] Temperature
    rho = 1000; %[kg/m3] Water density
    
    % Sediment parameters
    rhos = 2900; %[kg/m3] Sediment density
    D50 = logspace(-6,0,50); %[m] Grain size vector, this is not the model range
    Ha = 1e-20; %[J | Nm] Hamaker constant, order of e-20
    CDryB = 1750; %[kg/m3] Dry bed density
    C1 = 20; %18-24, smooth to rough sphere, affects small diameters, Ferguson and Church 2004
    C2 = 1; %0.4-1.2, smooth to rough sphere, affects big diameters, Ferguson and Church 2004
    phi = 30; %angle of repose in degrees
    
    % Derived input parameters
    TK = TC+273.15; %[degrees K] Temperature
    mu = 1.856e-14*exp(4209/TK+0.04527*TK-3.376e-5*TK^2); %[Pa | Ns/m2] Dynamic viscosity
    v = mu/rho; %[m2/s] Viscosity
    
    R = (rhos-rho)/rho; %[-] Relative density
    n = 1-CDryB/rhos; %[-] Porosity
    kk = ((1-(1-pi/(3*2^0.5)))/(1-n))^(1/3)-1; %[-] ?
end

%% allocate
for a = 1
    Lh = length(h);
    LQ = length(Q);
    Lg = length(grav);
    LD = length(D50);
    
    Rw = NaN(1,Lh);
    C = NaN(Lg,Lh);
    u = NaN(Lg,Lh);
    Q_h = NaN(Lg,Lh);
    Fr = NaN(Lg,Lh);
    Re = NaN(Lg,Lh);
    tau = NaN(Lg,Lh);
    ust = NaN(Lg,Lh);
    %lamlyr = NaN(Lg,Lh);
    h_Q = NaN(Lg,LQ);
    Rw_Q = NaN(Lg,LQ);
    C_Q = NaN(Lg,LQ);
    u_Q = NaN(Lg,LQ);
    Fr_Q = NaN(Lg,LQ);
    Re_Q = NaN(Lg,LQ);
    tau_Q = NaN(Lg,LQ);
    ust_Q = NaN(Lg,LQ);
    %lamlyr_Q = NaN(Lg,LQ);
    
    %Fc = NaN(1,LD);
    %Fw = NaN(Lg,LD);
    ws = NaN(Lg,LD);
    Rep = NaN(Lg,LD);
    Dst = NaN(Lg,LD);
    Rest = NaN(Lg,LD,Lh);
    shield = NaN(Lg,LD,Lh);
    lambda = NaN(Lg,LD,Lh);
    AdvL = NaN(Lg,LD,Lh);
    Rest_Q = NaN(Lg,LD,LQ);
    shield_Q = NaN(Lg,LD,LQ);
    lambda_Q = NaN(Lg,LD,LQ);
    AdvL_Q = NaN(Lg,LD,LQ);
    
    Rest_ini = NaN(Lg,LD);
    shield_cr_Z2003 = NaN(Lg,LD);
    shield_cr_Z2003b = NaN(Lg,LD);
    shield_cr_S1997 = NaN(Lg,LD);
    shield_cr_S1997b = NaN(Lg,LD);
    shield_cr_P2001 = NaN(Lg,LD);
    shield_cr_P2001b = NaN(Lg,LD);
    shield_cr_B1981 = NaN(Lg,LD);
    shield_cr_B1981b = NaN(Lg,LD);
    shield_cr_vR2007 = NaN(Lg,LD);
    shield_cr_C2006 = NaN(Lg,LD);
    shield_cr_M1977 = NaN(Lg,LD);
    shield_cr_S2014 = NaN(Lg,LD);
    shield_cr_B2008 = NaN(Lg,LD);
    shield_cr_K1986 = NaN(Lg,LD);
    shield_cr_K1986b = NaN(Lg,LD);
    shield_cr_K1986c = NaN(Lg,LD);
    shield_cr_C1982 = NaN(Lg,LD);
    shield_cr_K1986d = NaN(Lg,LD);
    
    lambda_cr_Z2003 = NaN(Lg,LD);
    lambda_cr_Z2003b = NaN(Lg,LD);
    lambda_cr_S1997 = NaN(Lg,LD);
    lambda_cr_S1997b = NaN(Lg,LD);
    lambda_cr_P2001 = NaN(Lg,LD);
    lambda_cr_P2001b = NaN(Lg,LD);
    lambda_cr_B1981 = NaN(Lg,LD);
    lambda_cr_B1981b = NaN(Lg,LD);
    lambda_cr_vR2007 = NaN(Lg,LD);
    lambda_cr_C2006 = NaN(Lg,LD);
    lambda_cr_M1977 = NaN(Lg,LD);
    lambda_cr_S2014 = NaN(Lg,LD);
    lambda_cr_B2008 = NaN(Lg,LD);
    lambda_cr_K1986 = NaN(Lg,LD);
    lambda_cr_K1986b = NaN(Lg,LD);
    lambda_cr_K1986c = NaN(Lg,LD);
    lambda_cr_C1982 = NaN(Lg,LD);
    lambda_cr_K1986d = NaN(Lg,LD);
    
    tau_cr_Z2003 = NaN(Lg,LD);
    tau_cr_Z2003b = NaN(Lg,LD);
    tau_cr_S1997 = NaN(Lg,LD);
    tau_cr_S1997b = NaN(Lg,LD);
    tau_cr_P2001 = NaN(Lg,LD);
    tau_cr_P2001b = NaN(Lg,LD);
    tau_cr_B1981 = NaN(Lg,LD);
    tau_cr_B1981b = NaN(Lg,LD);
    tau_cr_vR2007 = NaN(Lg,LD);
    tau_cr_C2006 = NaN(Lg,LD);
    tau_cr_M1977 = NaN(Lg,LD);
    tau_cr_S2014 = NaN(Lg,LD);
    tau_cr_B2008 = NaN(Lg,LD);
    tau_cr_K1986 = NaN(Lg,LD);
    tau_cr_K1986b = NaN(Lg,LD);
    tau_cr_K1986c = NaN(Lg,LD);
    tau_cr_C1982 = NaN(Lg,LD);
    tau_cr_K1986d = NaN(Lg,LD);
    
    ust_cr_Z2003 = NaN(Lg,LD);
    ust_cr_Z2003b = NaN(Lg,LD);
    ust_cr_S1997 = NaN(Lg,LD);
    ust_cr_S1997b = NaN(Lg,LD);
    ust_cr_P2001 = NaN(Lg,LD);
    ust_cr_P2001b = NaN(Lg,LD);
    ust_cr_B1981 = NaN(Lg,LD);
    ust_cr_B1981b = NaN(Lg,LD);
    ust_cr_vR2007 = NaN(Lg,LD);
    ust_cr_C2006 = NaN(Lg,LD);
    ust_cr_M1977 = NaN(Lg,LD);
    ust_cr_S2014 = NaN(Lg,LD);
    ust_cr_B2008 = NaN(Lg,LD);
    ust_cr_K1986 = NaN(Lg,LD);
    ust_cr_K1986b = NaN(Lg,LD);
    ust_cr_K1986c = NaN(Lg,LD);
    ust_cr_C1982 = NaN(Lg,LD);
    ust_cr_K1986d = NaN(Lg,LD);
end

%% Hydro parameters Q_mars = Q_earth
h_guess = 3;
for i = 1:LQ
    for g = 1:Lg
        h_Q(g,i) = h_guess+1;
        while abs(h_guess-h_Q(g,i))>0.0001
            h_guess = h_Q(g,i);
            Rw_Q(g,i) = (h_guess*W)/(h_guess+h_guess+W);                %[m] Hydraulic radius
            C_Q(g,i) = 5.75*grav(g)^0.5*log10(12*h_guess/(3*Drough));   %[m0.5/s] Chezy roughness
            u_Q(g,i) = C_Q(g,i)*(Rw_Q(g,i)*S)^0.5;                 %[m/s] Velocity
            h_Q(g,i) = Q(i)/(W*u_Q(g,i));                               %[m] Water depth
        end
        Fr_Q(g,i) = u_Q(g,i)/(grav(g)*h_Q(g,i))^0.5;                    %[-] Froude number
        Re_Q(g,i) = u_Q(g,i)*h_Q(g,i)/v;                                %[-] Reynolds number
        tau_Q(g,i) = rho*grav(g)*Rw_Q(g,i)*S;                      %[N/m2] Bed shear stress
        ust_Q(g,i) = (tau_Q(g,i)/rho)^0.5;                              %[m/s] Shear velocity
        %lamlyr_Q(g) = 11.63*v/ust_Q(g);                                 %[m] Laminar sublayer thickness
    end
end

%% Hydro parameters h_mars = h_earth
for i = 1:Lh
    Rw(i) = (h(i)*W)/(h(i)+h(i)+W);                             %[m] Hydraulic radius
    for g = 1:Lg
        C(g,i) = 5.74*grav(g)^0.5*log10(12.2*h(i)/(3*Drough));  %[m0.5/s] Chezy roughness
        u(g,i) = C(g,i)*(Rw(i)*S)^0.5;                     %[m/s] Velocity
        Q_h(g,i) = h(i)*W*u(g,i);                               %[m3/s] Discharge
        Fr(g,i) = u(g,i)/(grav(g)*h(i))^0.5;                    %[-] Froude number
        Re(g,i) = u(g,i)*h(i)/v;                                %[-] Reynolds number
        tau(g,i) = rho*grav(g)*Rw(i)*S;                    %[N/m2] Bed shear stress
        ust(g,i) = (tau(g,i)/rho)^0.5;                          %[m/s] Shear velocity
        %lamlyr(g) = 11.63*v/ust(g);                             %[m] Laminar sublayer thickness
    end
end

%% Sediment parameters
for d = 1:LD
    %Fc(d) = (3*Ha/(12*kk^2*D50(d)))*(1-cos(52.5*pi/180)); %Cohesive force. Lapotre 2019
    for g = 1:Lg
        %Fw(g,d) = (pi/6)*(rhos-rho)*grav(g)*D50(d)^3; %? force, Lapotre 2019
        ws(g,d) = R*grav(g)*D50(d)^2/(C1*v+(0.75*C2*R*grav(g)*D50(d)^3)^0.5); %[m/s] Settling velocity (Ferguson and Church 2004)
        Rep(g,d) = D50(d)^1.5*(R*grav(g))^0.5/v; %[-] Particle Reynolds number (Klein05,Leeuw20,Miedema), >70 turbulent <3.5 laminar
        Dst(g,d) = D50(d)*(R*grav(g)/v^2)^(1/3); %[-] Dimensionless particle parameter/Bonnefille number
        
        % Q_mars = Q_earth
        for i = 1:LQ
            Rest_Q(g,d,i) = D50(d)*ust_Q(g,i)/v;                            %[-] Particle Reynolds number (vR84,Leeuw20), Reynolds shear velocity number (Klein05), Boundary Reynolds number (Miedema10)
            shield_Q(g,d,i) = tau_Q(g,i)/((rhos-rho)*grav(g)*D50(d));       %[-] Shields parameter/Particle mobility parameter
            lambda_Q(g,d,i) = ust_Q(g,i)/ws(g,d);                           %[-] Movability number (Liu 1958)
            AdvL_Q(g,d,i) = u_Q(g,i)*h_Q(g,i)/ws(g,d);                      %[m] Advection length
        end
        % h_mars = h_earth
        for i = 1:Lh
            Rest(g,d,i) = D50(d)*ust(g,i)/v;                                %[-] Particle Reynolds number (vR84,Leeuw20,nino03), Reynolds shear velocity number (Klein05), Boundary Reynolds number (Miedema10), fluid reynolds number (Bagnold1966)
            shield(g,d,i) = tau(g,i)/((rhos-rho)*grav(g)*D50(d));           %[-] Shields parameter/Particle mobility parameter
            lambda(g,d,i) = ust(g,i)/ws(g,d);                               %[-] Movability number (Liu 1958)
            AdvL(g,d,i) = u(g,i)*h(i)/ws(g,d);                              %[m] Advection length
        end
    end
end

%% Critical Shields curves
i = 1;
acc = 1e-5;
shield_cr_L2006 = 0.15*S^0.25; %Lamb et al 2006 [gravel]
for d = 1:LD
    for g = 1:Lg
        %Zanke 2003 fit as in Kleinhans 2005
        if D50(d)>0.00005 && D50(d)<0.005
            shield_cr_Z2003(g,d) = 0.145*Rep(g,d)^-0.33+0.045*10^(-1100*Rep(g,d)^-1.5);
        else
            shield_cr_Z2003(g,d) = NaN;
        end
        
        %Zanke 2003 iterated
        shieldini = shield_Q(g,d,i);
        minimise = 1;
        while max(minimise)>acc
            tau_Q_ini = shieldini *((rhos-rho)*grav(g)*D50(d));
            ust_Q_ini = (tau_Q_ini/rho)^0.5;
            Rest_ini(g,d) = D50(d)*ust_Q_ini/v;
            uprmsb_ust = 0.31*Rest_ini(g,d)*exp(-0.1*Rest_ini(g,d))+1.8*exp(-0.88*D50(d)/(shieldini*R*D50(d)/sin(S)))*(1-exp(-0.1*Rest_ini(g,d)));
            Pt = 1 - exp(-0.08*Rest_ini(g,d));
            B = (1-Pt)*(2.5.*log(Rest_ini(g,d))+5.25)+8.5*Pt;
            uy_ust = ((1-Pt)/Rest_ini(g,d)^2+Pt/(2.5*log(1)+B^2))^-0.5;
            ub_ust = 0.8+0.9*uy_ust;
            uprmsb_ub = uprmsb_ust / ub_ust;
            K = 1 + 3e-8/((rhos-rho)*D50(d)^2);
            shield2 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
            minimise = abs(shield2 - shieldini);
            shieldini = shield2;
        end
        shield_cr_Z2003b(g,d) = shieldini;
        
        %Soulsby 1997 as in Kleinhans 2017 + Soulsby and Whitehouse 1997 as in Lapotre 2019 & Miedema 2010
        shield_cr_S1997(g,d) = 0.3/(1+1.2*Dst(g,d))+0.055*(1-exp(-0.02*Dst(g,d)));
        
        %Soulsby 1997 as in Kleinhans 2017 minus typo brackets
        shield_cr_S1997b(g,d) = 0.5*(0.3/(1+1.2*Dst(g,d))+0.055*(1-exp(-0.02*Dst(g,d))));
        
        %Paphitis 2001 - 1
        shieldini = shield_Q(g,d,i);
        minimise = 1;
        while max(minimise)>acc
            tau_Q_ini = shieldini *((rhos-rho)*grav(g)*D50(d));
            ust_Q_ini = (tau_Q_ini/rho)^0.5;
            Rest_ini(g,d) = D50(d)*ust_Q_ini/v;
            shield2 = 0.188/(1+Rest_ini(g,d))+0.0475*(1-0.699*exp(-0.015*Rest_ini(g,d)));
            minimise = abs(shield2 - shieldini);
            shieldini = shield2;
        end
        shield_cr_P2001(g,d) = shieldini;
        shield_cr_P2001(Rest_ini>10^5) = NaN;
        shield_cr_P2001(Rest_ini<0.01) = NaN;
        
        %Paphitis 2001 - 2
        shieldini = shield_Q(g,d,i);
        minimise = 1;
        while max(minimise)>acc
            tau_Q_ini = shieldini *((rhos-rho)*grav(g)*D50(d));
            ust_Q_ini = (tau_Q_ini/rho)^0.5;
            Rest_ini(g,d) = D50(d)*ust_Q_ini/v;
            shield2 = 0.273/(1+1.2*Dst(g,d))+0.046*(1-0.576*exp(-0.02*Dst(g,d)));
            minimise = abs(shield2 - shieldini);
            shieldini = shield2;
        end
        shield_cr_P2001b(g,d) = shieldini;
        shield_cr_P2001b(Rest_ini>10^4) = NaN;
        shield_cr_P2001b(Rest_ini<0.1) = NaN;
        
        %Brownlie 1981
        shield_cr_B1981(g,d) = 0.22*Rep(g,d)^-0.6+0.06*10^(-7.7*Rep(g,d)^-0.6);
        shield_cr_B1981(Rep<3) = NaN;
        shield_cr_B1981(Rep>9000) = NaN;
        
        %Brownlie 1981 as in Miedema 2010 & Righetti 2007
        shield_cr_B1981b(g,d) = 0.22*Rep(g,d)^-0.9+0.06*exp(-17.77*Rep(g,d)^-0.9);
        
        %van Rijn 2007
        if Dst(g,d)<4
            shield_cr_vR2007(g,d) = 0.115*Dst(g,d)^-0.5;
        elseif Dst(g,d)<10 && Dst(g,d)>=4
            shield_cr_vR2007(g,d) = 0.14*Dst(g,d)^-0.64;
        else
            shield_cr_vR2007(g,d) = NaN;
        end
        
        %Cao et al 2006
        if Rep(g,d)>=6.61 && Rep(g,d)<=282.84
            shield_cr_C2006(g,d) = (1+(0.0223*Rep(g,d))^2.8358)^0.3542/(3.0946*Rep(g,d)^0.6769);
        elseif Rep(g,d)<6.61
            shield_cr_C2006(g,d) = 0.1414*Rep(g,d)^(-0.2306);
        elseif Rep(g,d)>282.84
            shield_cr_C2006(g,d) = 0.045;
        else
            shield_cr_C2006(g,d) = NaN;
        end
        
        %Mantz 1977 as in Komar 1986 & Paphitis 2001
        shieldini = shield_Q(g,d,i);
        minimise = 1;
        while max(minimise)>acc
            tau_Q_ini = shieldini *((rhos-rho)*grav(g)*D50(d));
            ust_Q_ini = (tau_Q_ini/rho)^0.5;
            Rest_ini(g,d) = D50(d)*ust_Q_ini/v;
            shield2 = 0.1*Rest_ini(g,d)^-0.3;
            minimise = abs(shield2 - shieldini);
            shieldini = shield2;
        end
        shield_cr_M1977(g,d) = shieldini;
        shield_cr_M1977(Rest_ini>1) = NaN;
        shield_cr_M1977(Rest_ini<0.03) = NaN;
    end
end

for d = 1:LD %based on mobility number
    for g = 1:Lg
        %Simoes 2014
        lambda_cr_S2014(g,d) = 0.215 +(6.79/Dst(g,d)^1.7) - (0.075*exp(-2.62*10^-3*Dst(g,d)));
        shield_cr_S2014(g,d) = (lambda_cr_S2014(g,d)*ws(g,d))^2 / (R*grav(g)*D50(d));
        
        %Beheshti amd Ataie-Ashtiani 2008
        if Dst(g,d)>10 && Dst(g,d)<500
            lambda_cr_B2008(g,d) = 0.4738*Dst(g,d)^-0.226;
        elseif Dst(g,d)<=10 && Dst(g,d)>0.4
            lambda_cr_B2008(g,d) = 9.6674*Dst(g,d)^-1.57;
        else
            lambda_cr_B2008(g,d) = NaN;
        end
        shield_cr_B2008(g,d) = (lambda_cr_B2008(g,d)*ws(g,d))^2 / (R*grav(g)*D50(d));
        
        %Komar 1986 - 1
        lambdaP = lambda_Q(g,d,i);
        lambdaini = lambdaP+0.1;
        while abs(lambdaP-lambdaini)<acc
            lambdaini = lambdaP;
            ust_Q_ini = lambdaini*ws(g,d);
            Rest_ini(g,d) = D50(d)*ust_Q_ini/v;
            lambdaP = 1.8*Rest_ini(g,d)^-1.3;
        end
        lambda_cr_K1986(g,d) = lambdaP;
        lambda_cr_K1986(Rest_ini>1) = NaN;
        shield_cr_K1986(g,d) = (lambda_cr_K1986(g,d)*ws(g,d))^2 / (R*grav(g)*D50(d));
        
        %Komar 1986 - 2
        lambdaP = lambda_Q(g,d,i);
        lambdaini = lambdaP+0.1;
        while abs(lambdaP-lambdaini)<acc
            lambdaini = lambdaP;
            ust_Q_ini = lambdaini*ws(g,d);
            Rest_ini(g,d) = D50(d)*ust_Q_ini/v;
            lambdaP = 1.14*Rest_ini(g,d)^-1.37;
        end
        lambda_cr_K1986b(g,d) = lambdaP;
        lambda_cr_K1986b(Rest_ini>1) = NaN;
        shield_cr_K1986b(g,d) = (lambda_cr_K1986b(g,d)*ws(g,d))^2 / (R*grav(g)*D50(d));
        
        %%Komar 1986 - 3
        lambda_cr_K1986c(g,d) = 5.54*Rep(g,d)^-1.09;
        lambda_cr_K1986c(Rest_ini>1) = NaN;
        shield_cr_K1986c(g,d) = (lambda_cr_K1986c(g,d)*ws(g,d))^2 / (R*grav(g)*D50(d));
    end
end

for d = 1:LD %based on shear stress
    for g = 1:Lg 
        %Collins 1982 (quartz, Komar 1986)
        tau_cr_C1982(g,d) = 1.24*ws(g,d)^0.33;
        tau_cr_C1982(ws>0.1) = NaN;
        shield_cr_C1982(g,d) = tau_cr_C1982(g,d) / ((rhos-rho)*grav(g)*D50(d));       
    end
end

for d = 1:LD %based on shear velocity
    for g = 1:Lg       
        %Komar 1986 - 4 (generalised form of Collins 1982)
        ust_cr_K1986d(g,d) = 0.482*(R*grav(g)*v)^0.282*ws(g,d)^0.154;
        ust_cr_K1986d(ws>0.1) = NaN;
        tau_cr_K1986d(g,d) = ust_cr_K1986d(g,d)^2*rho;
        shield_cr_K1986d(g,d) = tau_cr_K1986d(g,d) / ((rhos-rho)*grav(g)*D50(d));      
    end
end

%% Critical shear stress
for d = 1:LD %omgerekend
    for g = 1:Lg
        tau_cr_Z2003(g,d) = shield_cr_Z2003(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_Z2003b(g,d) = shield_cr_Z2003b(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_S1997(g,d) = shield_cr_S1997(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_S1997b(g,d) = shield_cr_S1997b(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_P2001(g,d) = shield_cr_P2001(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_P2001b(g,d) = shield_cr_P2001b(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_B1981(g,d) = shield_cr_B1981(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_B1981b(g,d) = shield_cr_B1981b(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_vR2007(g,d) = shield_cr_vR2007(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_C2006(g,d) = shield_cr_C2006(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_M1977(g,d) = shield_cr_M1977(g,d) * (rhos-rho)*grav(g)*D50(d);
        
        tau_cr_S2014(g,d) = shield_cr_S2014(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_B2008(g,d) = shield_cr_B2008(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_K1986(g,d) = shield_cr_K1986(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_K1986b(g,d) = shield_cr_K1986b(g,d) * (rhos-rho)*grav(g)*D50(d);
        tau_cr_K1986c(g,d) = shield_cr_K1986c(g,d) * (rhos-rho)*grav(g)*D50(d);
    end
end
%In previous sections: tau_cr_C1982 tau_cr_K1986d

%% Critical shear velocity
for d = 1:LD %omgerekend
    for g = 1:Lg
        ust_cr_Z2003(g,d) = (tau_cr_Z2003(g,d) / rho)^0.5;
        ust_cr_Z2003b(g,d) = (tau_cr_Z2003b(g,d) / rho)^0.5;
        ust_cr_S1997(g,d) = (tau_cr_S1997(g,d) / rho)^0.5;
        ust_cr_S1997b(g,d) = (tau_cr_S1997b(g,d) / rho)^0.5;
        ust_cr_P2001(g,d) = (tau_cr_P2001(g,d) / rho)^0.5;
        ust_cr_P2001b(g,d) = (tau_cr_P2001b(g,d) / rho)^0.5;
        ust_cr_B1981(g,d) = (tau_cr_B1981(g,d) / rho)^0.5;
        ust_cr_B1981b(g,d) = (tau_cr_B1981b(g,d) / rho)^0.5;
        ust_cr_vR2007(g,d) = (tau_cr_vR2007(g,d) / rho)^0.5;
        ust_cr_C2006(g,d) = (tau_cr_C2006(g,d) / rho)^0.5;
        ust_cr_M1977(g,d) = (tau_cr_M1977(g,d) / rho)^0.5;
        
        ust_cr_S2014(g,d) = (tau_cr_S2014(g,d) / rho)^0.5;
        ust_cr_B2008(g,d) = (tau_cr_B2008(g,d) / rho)^0.5;
        ust_cr_K1986(g,d) = (tau_cr_K1986(g,d) / rho)^0.5;
        ust_cr_K1986b(g,d) = (tau_cr_K1986b(g,d) / rho)^0.5;
        ust_cr_K1986c(g,d) = (tau_cr_K1986c(g,d) / rho)^0.5;
        
        ust_cr_C1982(g,d) = (tau_cr_C1982(g,d) / rho)^0.5;
    end
end
%In previous sections: ust_cr_K1986d

%% Mobility parameter
for d = 1:LD %omgerekend
    for g = 1:Lg
        lambda_cr_Z2003(g,d) = ust_cr_Z2003(g,d) / ws(g,d);
        lambda_cr_Z2003b(g,d) = ust_cr_Z2003b(g,d) / ws(g,d);
        lambda_cr_S1997(g,d) = ust_cr_S1997(g,d) / ws(g,d);
        lambda_cr_S1997b(g,d) = ust_cr_S1997b(g,d) / ws(g,d);
        lambda_cr_P2001(g,d) = ust_cr_P2001(g,d) / ws(g,d);
        lambda_cr_P2001b(g,d) = ust_cr_P2001b(g,d) / ws(g,d);
        lambda_cr_B1981(g,d) = ust_cr_B1981(g,d) / ws(g,d);
        lambda_cr_B1981b(g,d) = ust_cr_B1981b(g,d) / ws(g,d);
        lambda_cr_vR2007(g,d) = ust_cr_vR2007(g,d) / ws(g,d);
        lambda_cr_C2006(g,d) = ust_cr_C2006(g,d) / ws(g,d);
        lambda_cr_M1977(g,d) = ust_cr_M1977(g,d) / ws(g,d);
        
        lambda_cr_C1982(g,d) = ust_cr_C1982(g,d) / ws(g,d);
        lambda_cr_K1986d(g,d) = ust_cr_K1986d(g,d) / ws(g,d);
    end
end
%In previous sections: lambda_S2014 lambda_B2008 lambda_K1986
%lambda_K1986b lambda_K1986c

%% Plot settings
for a = 1
    pl.width = 18;
    pl.height = 16;
    pl.line = 1;
    pl.line_ax = 0.75;
    pl.fsz = 7;
    pl.fsz2 = 8.5;
    CMAPOBJ = clrmap('read','Gravity.clrmap');
    clr1 = clrmap(CMAPOBJ,Lg);
    CMAPOBJ = clrmap('read','Discharge.clrmap');
    clr2 = clrmap(CMAPOBJ,Lh);
    clr3 = clrmap(CMAPOBJ,LQ);
    clrM = [1 0.804 0.804];
    clrE = [0.8 0.8 0.9804];
    pl.vertical = 2;
    pl.horizantal = 2;
    pl.lmarge = 0.065;
    pl.rmarge = 0.025;
    pl.bmarge = 0.06;
    pl.tmarge = 0.01;
    pl.intv = 0.065;
    pl.inth = 0.065;
    pl.wt = (1 - pl.lmarge - pl.rmarge - ((pl.horizantal-1)*pl.inth)) / pl.horizantal;
    pl.ht = (1 - pl.bmarge - pl.tmarge - ((pl.vertical-1)*pl.intv)) / pl.vertical;
    y = 0.94;
    xd = 10^(log10(D50(1))+( log10(D50(end))-log10(D50(1)) )*0.04);
end

%% plot Q_mars = Q_earth
for a = 1
    close all; clc;
    f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[10 1 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');
    
    s1(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    %plot(D50,ones(size(D50)).*shield_cr_L2006,':','color','b','linewidth',pl.line) %Lamb et al 2006 [gravel], makes no sense for this situation
    for g = 1:Lg
        plot(D50,(ws(g,:).^2*rho)./((rhos-rho)*grav(g)*D50),'--','color',clr1(g,:),'linewidth',pl.line)
        
        %plot(D50,shield_cr_Z2003(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Zanke 2003 as in Kleinhans 2005, overbodig door overlap Z2003b en kleinere range
        plot(D50,shield_cr_Z2003b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Zanke 2003 interated
        plot(D50,shield_cr_S1997(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Soulsby 1997 as in Kleinhans 2017 + Soulsby and Whitehouse 1997 as in Lapotre 2019
        %plot(D50,shield_cr_S1997b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Soulsby 1997 as in Kleinhans 2017 minus typo brackets, confusing compared to others
        plot(D50,shield_cr_P2001(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Paphitis 2001
        plot(D50,shield_cr_P2001b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Paphitis 2001
        plot(D50,shield_cr_B1981(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Brownlie 1981
        %plot(D50,shield_cr_B1981b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Brownlie 1981 as in Miedema 2010 & Righetti 2007, this seems wrong because tau_cr increases for smaller grain sizes.
        plot(D50,shield_cr_vR2007(g,:),':','color',clr1(g,:),'linewidth',pl.line); %van Rijn 2007
        plot(D50,shield_cr_C2006(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Cao et al 2006
        plot(D50,shield_cr_M1977(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Mantz 1977 as in Komar 1986 & Paphitis 2001
        plot(D50,shield_cr_S2014(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Simoes 2014
        %plot(D50,shield_cr_B2008(g,:),'-','color',clr1(g,:),'linewidth',pl.line); %Beheshti amd Ataie-Ashtiani 2008, part shows increasing tau_cr for smaller grains
        %plot(D50,shield_cr_K1986(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,shield_cr_K1986b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,shield_cr_K1986c(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,shield_cr_C1982(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Collins 1982 (quartz, Komar 1986), bad with gravity, Komar 1986 - 4 is better version
        plot(D50,shield_cr_K1986d(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986 (generalised form of Collins 1982)
        count = 1;
        for i = Qselect
            plot(D50,shield_Q(g,:,i),'color',[clr1(g,:) 1-(1/length(Qselect))*(count-1)],'linewidth',pl.line)
            count = count + 1;
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-3 10^3]);
    l(1) = ylabel('Shields parameter (\theta) [-]');
    l(2) = xlabel('Grain size (D_{50}) [m]');
    t(1) = text(xd,10^(-3+(3--3)*y),' a');
   
    s1(2) = axes('Position',[pl.lmarge+pl.wt+pl.inth 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    %plot(D50,ones(size(D50)).*shield_cr_L2006,':','color','b','linewidth',pl.line) %Lamb et al 2006 [gravel], makes no sense for this situation
    for g = 1:Lg
        plot(D50,ws(g,:).^2*rho,'--','color',clr1(g,:),'linewidth',pl.line)
        
        %plot(D50,tau_cr_Z2003(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Zanke 2003 as in Kleinhans 2005, overbodig door overlap Z2003b en kleinere range
        plot(D50,tau_cr_Z2003b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Zanke 2003 interated
        plot(D50,tau_cr_S1997(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Soulsby 1997 as in Kleinhans 2017 + Soulsby and Whitehouse 1997 as in Lapotre 2019
        %plot(D50,tau_cr_S1997b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Soulsby 1997 as in Kleinhans 2017 minus typo brackets, confusing compared to others
        plot(D50,tau_cr_P2001(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Paphitis 2001
        plot(D50,tau_cr_P2001b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Paphitis 2001
        plot(D50,tau_cr_B1981(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Brownlie 1981
        %plot(D50,tau_cr_B1981b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Brownlie 1981 as in Miedema 2010 & Righetti 2007, this seems wrong because tau_cr increases for smaller grain sizes.
        plot(D50,tau_cr_vR2007(g,:),':','color',clr1(g,:),'linewidth',pl.line); %van Rijn 2007
        plot(D50,tau_cr_C2006(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Cao et al 2006
        plot(D50,tau_cr_M1977(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Mantz 1977 as in Komar 1986 & Paphitis 2001        
        plot(D50,tau_cr_S2014(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Simoes 2014
        %plot(D50,tau_cr_B2008(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Beheshti amd Ataie-Ashtiani 2008, part shows increasing tau_cr for smaller grains
        %plot(D50,tau_cr_K1986(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,tau_cr_K1986b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,tau_cr_K1986c(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,tau_cr_C1982(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Collins 1982 (quartz, Komar 1986), bad with gravity, Komar 1986 - 4 is better version
        plot(D50,tau_cr_K1986d(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986 (generalised form of Collins 1982)
        count = 1;
        for i = Qselect
            plot([D50(1) D50(end)],[tau_Q(g,i) tau_Q(g,i)],'-','color',[clr1(g,:) 1-(1/length(Qselect))*(count-1)],'linewidth',pl.line);
            count = count+1;
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-3 10^3]);
    l(3) = ylabel('Bed shear stress (\tau) [N/m^2]');
    l(4) = xlabel('Grain size (D_{50}) [m]');
    t(2) = text(xd,10^(-3+(3--3)*y),' b');
    
    s1(3) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    %plot(D50,ones(size(D50)).*ust_cr_L2006,':','color','b','linewidth',pl.line) %Lamb et al 2006 [gravel], makes no sense for this situation
    for g = 1:Lg
        plot(D50,ws(g,:),'--','color',clr1(g,:),'linewidth',pl.line);
        
        %plot(D50,ust_cr_Z2003(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_Z2003b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_S1997(g,:),':','color',clr1(g,:),'linewidth',pl.line);
       % plot(D50,ust_cr_S1997b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_P2001(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_P2001b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_B1981(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_B1981b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_vR2007(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_C2006(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_M1977(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_S2014(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_B2008(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_K1986(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_K1986b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_K1986c(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_C1982(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_K1986d(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        count = 1;
        for i = Qselect
            plot([D50(1) D50(end)],[ust_Q(g,i) ust_Q(g,i)],'-','color',[clr1(g,:) 1-(1/length(Qselect))*(count-1)],'linewidth',pl.line)
            count = count+1;
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-3 10^0]);
    l(5) = ylabel('Shear velocity (u_{*}) [m/s]');
    l(6) = xlabel('Grain size (D_{50}) [m]');
    t(3) = text(xd,10^(-3+(0--3)*y),' c');
    
    s1(4) = axes('Position',[pl.lmarge+pl.wt+pl.inth pl.bmarge pl.wt pl.ht]);
    hold on
    %plot(D50,ones(size(D50)).*lambda_cr_L2006,':','color','b','linewidth',pl.line) %Lamb et al 2006 [gravel], makes no sense for this situation
    for g = 1:Lg
        count = 1;
        for i = Qselect
            plot(D50,lambda_Q(g,:,i),'-','color',[clr1(g,:) 1-(1/length(Qselect))*(count-1)],'linewidth',pl.line)
            count = count+1;
        end
    end
    plot([D50(1) D50(end)],[1 1],'--k','linewidth',pl.line)
    plot([D50(1) D50(end)],[100000 100000],':k','linewidth',pl.line)
    for g = 1:Lg
        %plot(D50,lambda_cr_Z2003(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_Z2003b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_S1997(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_S1997b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_P2001(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_P2001b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_B1981(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_B1981b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_vR2007(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_C2006(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_M1977(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_S2014(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_B2008(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_K1986(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_K1986b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_K1986c(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_C1982(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_K1986d(g,:),':','color',clr1(g,:),'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-1 10^3]);
    l(7) = ylabel('Movability number (k) [-]');
    l(8) = xlabel('Grain size (D_{50}) [m]');
    t(4) = text(xd,10^(-1+(3--1)*y),' d');
    
    if length(Qselect) == 1
    legend('Mars Q=2000 m^3/s','Earth Q=2000 m^3/s','Suspension threshold','Thresholds of motion', ...
        'location','northeast')
    else
    legend('Mars Q=500 m^3/s','Mars Q=15000 m^3/s','Earth Q=500 m^3/s','Earth Q=15000 m^3/s','Suspension threshold','Thresholds of motion', ...
        'location','northeast')
    end
    
    set(s1,'box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);
    
    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);
    
    set(gcf,'renderer','painters');
    print(f1,'-dpng',[output '/FlowOnMars2_QD50'],'-r400');
    print(f1,'-dpdf',[output '/FlowOnMars2_QD50'],'-r400');
end

%% plot h_mars = h_earth
for a = 1
    close all; clc;
    f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[10 1 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');
    
    s1(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    %plot(D50,ones(size(D50)).*shield_cr_L2006,':','color','b','linewidth',pl.line) %Lamb et al 2006 [gravel], makes no sense for this situation
    for g = 1:Lg
        plot(D50,(ws(g,:).^2*rho)./((rhos-rho)*grav(g)*D50),'--','color',clr1(g,:),'linewidth',pl.line)
        
        %plot(D50,shield_cr_Z2003(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Zanke 2003 as in Kleinhans 2005, overbodig door overlap Z2003b en kleinere range
        plot(D50,shield_cr_Z2003b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Zanke 2003 interated
        plot(D50,shield_cr_S1997(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Soulsby 1997 as in Kleinhans 2017 + Soulsby and Whitehouse 1997 as in Lapotre 2019
        %plot(D50,shield_cr_S1997b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Soulsby 1997 as in Kleinhans 2017 minus typo brackets, confusing compared to others
        plot(D50,shield_cr_P2001(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Paphitis 2001
        plot(D50,shield_cr_P2001b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Paphitis 2001
        plot(D50,shield_cr_B1981(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Brownlie 1981
        %plot(D50,shield_cr_B1981b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Brownlie 1981 as in Miedema 2010 & Righetti 2007, this seems wrong because tau_cr increases for smaller grain sizes.
        plot(D50,shield_cr_vR2007(g,:),':','color',clr1(g,:),'linewidth',pl.line); %van Rijn 2007
        plot(D50,shield_cr_C2006(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Cao et al 2006
        plot(D50,shield_cr_M1977(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Mantz 1977 as in Komar 1986 & Paphitis 2001
        plot(D50,shield_cr_S2014(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Simoes 2014
        %plot(D50,shield_cr_B2008(g,:),'-','color',clr1(g,:),'linewidth',pl.line); %Beheshti amd Ataie-Ashtiani 2008, part shows increasing tau_cr for smaller grains
        %plot(D50,shield_cr_K1986(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,shield_cr_K1986b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,shield_cr_K1986c(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,shield_cr_C1982(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Collins 1982 (quartz, Komar 1986), bad with gravity, Komar 1986 - 4 is better version
        plot(D50,shield_cr_K1986d(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986 (generalised form of Collins 1982)
        count = 1;
        for i = hselect
            plot(D50,shield(g,:,i),'color',[clr1(g,:) 1-(1/length(hselect))*(count-1)],'linewidth',pl.line)
            count = count + 1;
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-3 10^3]);
    l(1) = ylabel('Shields parameter (\theta) [-]');
    l(2) = xlabel('Grain size (D_{50}) [m]');
    t(1) = text(xd,10^(-3+(3--3)*y),' a');
   
    s1(2) = axes('Position',[pl.lmarge+pl.wt+pl.inth 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    %plot(D50,ones(size(D50)).*shield_cr_L2006,':','color','b','linewidth',pl.line) %Lamb et al 2006 [gravel], makes no sense for this situation
    for g = 1:Lg
        plot(D50,ws(g,:).^2*rho,'--','color',clr1(g,:),'linewidth',pl.line)
        
        %plot(D50,tau_cr_Z2003(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Zanke 2003 as in Kleinhans 2005, overbodig door overlap Z2003b en kleinere range
        plot(D50,tau_cr_Z2003b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Zanke 2003 interated
        plot(D50,tau_cr_S1997(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Soulsby 1997 as in Kleinhans 2017 + Soulsby and Whitehouse 1997 as in Lapotre 2019
        %plot(D50,tau_cr_S1997b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Soulsby 1997 as in Kleinhans 2017 minus typo brackets, confusing compared to others
        plot(D50,tau_cr_P2001(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Paphitis 2001
        plot(D50,tau_cr_P2001b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Paphitis 2001
        plot(D50,tau_cr_B1981(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Brownlie 1981
        %plot(D50,tau_cr_B1981b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Brownlie 1981 as in Miedema 2010 & Righetti 2007, this seems wrong because tau_cr increases for smaller grain sizes.
        plot(D50,tau_cr_vR2007(g,:),':','color',clr1(g,:),'linewidth',pl.line); %van Rijn 2007
        plot(D50,tau_cr_C2006(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Cao et al 2006
        plot(D50,tau_cr_M1977(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Mantz 1977 as in Komar 1986 & Paphitis 2001        
        plot(D50,tau_cr_S2014(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Simoes 2014
        %plot(D50,tau_cr_B2008(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Beheshti amd Ataie-Ashtiani 2008, part shows increasing tau_cr for smaller grains
        %plot(D50,tau_cr_K1986(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,tau_cr_K1986b(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,tau_cr_K1986c(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986, Komar 1986 - 4 is better version
        %plot(D50,tau_cr_C1982(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Collins 1982 (quartz, Komar 1986), bad with gravity, Komar 1986 - 4 is better version
        plot(D50,tau_cr_K1986d(g,:),':','color',clr1(g,:),'linewidth',pl.line); %Komar 1986 (generalised form of Collins 1982)
        count = 1;
        for i = hselect
            plot([D50(1) D50(end)],[tau(g,i) tau(g,i)],'-','color',[clr1(g,:) 1-(1/length(hselect))*(count-1)],'linewidth',pl.line);
            count = count+1;
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-3 10^3]);
    l(3) = ylabel('Bed shear stress (\tau) [N/m^2]');
    l(4) = xlabel('Grain size (D_{50}) [m]');
    t(2) = text(xd,10^(-3+(3--3)*y),' b');
    
    s1(3) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    %plot(D50,ones(size(D50)).*ust_cr_L2006,':','color','b','linewidth',pl.line) %Lamb et al 2006 [gravel], makes no sense for this situation
    for g = 1:Lg
        plot(D50,ws(g,:),'--','color',clr1(g,:),'linewidth',pl.line);
        
        %plot(D50,ust_cr_Z2003(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_Z2003b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_S1997(g,:),':','color',clr1(g,:),'linewidth',pl.line);
       % plot(D50,ust_cr_S1997b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_P2001(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_P2001b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_B1981(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_B1981b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_vR2007(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_C2006(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_M1977(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_S2014(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_B2008(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_K1986(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_K1986b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_K1986c(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,ust_cr_C1982(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,ust_cr_K1986d(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        count = 1;
        for i = hselect
            plot([D50(1) D50(end)],[ust(g,i) ust(g,i)],'-','color',[clr1(g,:) 1-(1/length(hselect))*(count-1)],'linewidth',pl.line)
            count = count+1;
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-3 10^0]);
    l(5) = ylabel('Shear velocity (u_{*}) [m/s]');
    l(6) = xlabel('Grain size (D_{50}) [m]');
    t(3) = text(xd,10^(-3+(0--3)*y),' c');
    
    s1(4) = axes('Position',[pl.lmarge+pl.wt+pl.inth pl.bmarge pl.wt pl.ht]);
    hold on
    %plot(D50,ones(size(D50)).*lambda_cr_L2006,':','color','b','linewidth',pl.line) %Lamb et al 2006 [gravel], makes no sense for this situation
    for g = 1:Lg
        count = 1;
        for i = hselect
            plot(D50,lambda(g,:,i),'-','color',[clr1(g,:) 1-(1/length(hselect))*(count-1)],'linewidth',pl.line)
            count = count+1;
        end
    end
    plot([D50(1) D50(end)],[1 1],'--k','linewidth',pl.line)
    plot([D50(1) D50(end)],[100000 100000],':k','linewidth',pl.line)
    for g = 1:Lg
        %plot(D50,lambda_cr_Z2003(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_Z2003b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_S1997(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_S1997b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_P2001(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_P2001b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_B1981(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_B1981b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_vR2007(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_C2006(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_M1977(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_S2014(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_B2008(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_K1986(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_K1986b(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_K1986c(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        %plot(D50,lambda_cr_C1982(g,:),':','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,lambda_cr_K1986d(g,:),':','color',clr1(g,:),'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-1 10^3]);
    l(7) = ylabel('Movability number (k) [-]');
    l(8) = xlabel('Grain size (D_{50}) [m]');
    t(4) = text(xd,10^(-1+(3--1)*y),' d');
    
    if length(hselect) == 1
    legend('Mars h=5 m','Earth h=5 m','Suspension threshold','Thresholds of motion', ...
        'location','northeast')
    else
    legend('Mars h=3 m^3/s','Mars h=15 m','Earth h=3 m^3/s','Earth h=15 m','Suspension threshold','Thresholds of motion', ...
        'location','northeast')
    end
    
    set(s1,'box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);
    
    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);
    
    set(gcf,'renderer','painters');
    print(f1,'-dpng',[output '/FlowOnMars2_hD50'],'-r400');
    print(f1,'-dpdf',[output '/FlowOnMars2_hD50'],'-r400');
end
