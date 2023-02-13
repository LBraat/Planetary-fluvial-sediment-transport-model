%% Flow on Mars 4 - Transport Q (1 Q equal for Mars and Earth)
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
    Q = 500:500:15000; %[m3/s] Discharge
    Qselect = [4];
    W = 200; %[m] Channel width
    %grav = 1:0.5:12; %[m/s2] Gravitational acceleration
    grav = [3.7 9.8]; %[m/s2] Gravitational acceleration
    S = 0.001; %[m/m] Slope
    TC = 4; %[degrees C] Temperature
    rho = 1000; %[kg/m3] Water density
    
    % Sediment parameters
    rhos = 2900; %[kg/m3] Sediment density
    D50 = logspace(-6,0,100); %[m] Grain size vector, this is not the model range
    Ha = 1e-20; %[J | Nm] Hamaker constant, order of e-20
    CDryB = 1750; %[kg/m3] Dry bed density
    C1 = 20; %18-24, smooth to rough sphere, affects small diameters, Ferguson and Church 2004
    C2 = 1; %0.4-1.2, smooth to rough sphere, affects big diameters, Ferguson and Church 2004
    phi = 30; %angle of repose in degrees
    Kappa = 0.41; %[-] Kappa
    
    % Derived input parameters
    TK = TC+273.15; %[degrees K] Temperature
    mu = 1.856e-14*exp(4209/TK+0.04527*TK-3.376e-5*TK^2); %[Pa | Ns/m2] Dynamic viscosity
    v = mu/rho; %[m2/s] Viscosity    
    R = (rhos-rho)/rho; %[-] Relative density
    n = 1-CDryB/rhos; %[-] Porosity
    kk = ((1-(1-pi/(3*2^0.5)))/(1-n))^(1/3)-1; %[-] ?
    ks = 3*Drough; %[m] Nikuradse
end

%% allocate
for a = 1
    LQ = length(Q);
    Lg = length(grav);
    LD = length(D50);
    
    h = NaN(Lg,LQ);
    Rw = NaN(Lg,LQ);
    C = NaN(Lg,LQ);
    u = NaN(Lg,LQ);
    Fr = NaN(Lg,LQ);
    Re = NaN(Lg,LQ);
    tau = NaN(Lg,LQ);
    ust = NaN(Lg,LQ);
    %lamlyr = NaN(Lg,Lh);
    
    %Fc = NaN(1,LQ);
    %Fw = NaN(Lg,LQ);
    ws = NaN(Lg,LD);
    Rep = NaN(Lg,LD);
    Dst = NaN(Lg,LD);
    Rest = NaN(Lg,LD);
    shield = NaN(Lg,LD);
    lambda = NaN(Lg,LD);
    AdvL = NaN(Lg,LD);
    
    shield_cr = NaN(Lg,LD);
    tau_cr = NaN(Lg,LD);
    ust_cr = NaN(Lg,LD);
    lambda_cr = NaN(Lg,LD);
    S0 = NaN(Lg,LD);
    
    Qb_MPM = NaN(Lg,LD);
    qb_MPM = NaN(Lg,LD);
    Qb_MP = NaN(Lg,LD);
    qb_MP = NaN(Lg,LD);
    Qb_E50 = NaN(Lg,LD);
    qb_E50 = NaN(Lg,LD);
    Qb_FB = NaN(Lg,LD);
    qb_FB = NaN(Lg,LD);
    Qb_R = NaN(Lg,LD);
    qb_R = NaN(Lg,LD);
    Qb_HJ = NaN(Lg,LD);
    qb_HJ = NaN(Lg,LD);
    Qb_WP1 = NaN(Lg,LD);
    Qb_WP2 = NaN(Lg,LD);
    qb_WP1 = NaN(Lg,LD);
    qb_WP2 = NaN(Lg,LD);
    Qb_Ca = NaN(Lg,LD);
    qb_Ca = NaN(Lg,LD);
    Qb_C = NaN(Lg,LD);
    qb_C = NaN(Lg,LD);
    qb_N = NaN(Lg,LD);
    Qb_N = NaN(Lg,LD);
    Qb_S = NaN(Lg,LD);
    qb_S = NaN(Lg,LD);
    Qb_AM = NaN(Lg,LD);
    qb_AM = NaN(Lg,LD);
    p_EF = NaN(Lg,LD);
    Qb_EF = NaN(Lg,LD);
    qb_EF = NaN(Lg,LD);
    Qb_P = NaN(Lg,LD);
    qb_P = NaN(Lg,LD);
    Qb_vR = NaN(Lg,LD);
    qb_vR = NaN(Lg,LD);
    Qb_vR2 = NaN(Lg,LD);
    qb_vR2 = NaN(Lg,LD);
    qb_EH = NaN(Lg,LD);
    Qb_EH = NaN(Lg,LD);
    qb_W = NaN(Lg,LD);
    Qb_W = NaN(Lg,LD);
    Qb_E42 = NaN(Lg,LD);
    qb_E42 = NaN(Lg,LD);
    G2 = NaN(Lg,LD);
    eb = NaN(Lg,LD);
    qb_B = NaN(Lg,LD);
    Qb_B = NaN(Lg,LD);
    
    a_E50 = NaN(1,LD);
    a_EF = NaN(1,LD);
    a_IK = NaN(1,Lg);
    a_CR = NaN(1,Lg);
    a_AF = NaN(1,Lg);
    a_GP = NaN(1,Lg);
    a_WP = NaN(1,Lg);
    a_dL = NaN(1,Lg);
    a_vR = NaN(1,Lg);
    a_SM = NaN(Lg,LD);
    a_ML = NaN(Lg,LD);
    
    Rouse = NaN(Lg,LD);
    Rouse_dL = NaN(Lg,LD);
    
    Ca_E50 = NaN(Lg,LD);
    labda_EF = NaN(Lg,LD);
    Ca_EF = NaN(Lg,LD);
    Ca_SM = NaN(Lg,LD);
    Ca_ML = NaN(Lg,LD);
    Ca_vR = NaN(Lg,LD);
    Z_AF = NaN(Lg,LD);
    Ca_AF = NaN(Lg,LD);
    Z_GP = NaN(Lg,LD);
    Ca_GP = NaN(Lg,LD);
    Z_WP = NaN(Lg,LD);
    Ca_WP = NaN(Lg,LD);
    Ca_dL = NaN(Lg,LD);
    I_CR = NaN(Lg,LD);
    Cm_CR = NaN(Lg,LD);
    Ca_CR = NaN(Lg,LD);
    A_IK = NaN(Lg,LD);
    Omega_IK = NaN(Lg,LD);
    Ca_IK = NaN(Lg,LD);
    
    qs_E50 = NaN(Lg,LD);
    qs_EF = NaN(Lg,LD);
    qs_IK = NaN(Lg,LD);
    qs_CR = NaN(Lg,LD);
    qs_AF = NaN(Lg,LD);
    qs_GP = NaN(Lg,LD);
    qs_WP = NaN(Lg,LD);
    qs_dL = NaN(Lg,LD);
    qs_vR = NaN(Lg,LD);
    qs_SM = NaN(Lg,LD);
    qs_ML = NaN(Lg,LD);
end

%% Hydro parameters Q_mars = Q_earth
h_guess = 3;
for i = 1:LQ
    for g = 1:Lg
        h(g,i) = h_guess+1;
        while abs(h_guess-h(g,i))>0.0001
            h_guess = h(g,i);
            Rw(g,i) = (h_guess*W)/(h_guess+h_guess+W); %[m] Hydraulic radius
            C(g,i) = 5.75*grav(g)^0.5*log10(12*h_guess/(ks)); %[m0.5/s] Chezy roughness
            u(g,i) = C(g,i)*(Rw(g,i)*S)^0.5; %[m/s] Velocity
            h(g,i) = Q(i)/(W*u(g,i)); %[m] Water depth
        end
        Fr(g,i) = u(g,i)/(grav(g)*h(g,i))^0.5; %[-] Froude number
        Re(g,i) = u(g,i)*h(g,i)/v; %[-] Reynolds number
        tau(g,i) = rho*grav(g)*Rw(g,i)*S; %[N/m2] Bed shear stress
        ust(g,i) = (tau(g,i)/rho)^0.5; %[m/s] Shear velocity
        %lamlyr_Q(g) = 11.63*v/ust_Q(g); %[m] Laminar sublayer thickness
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
        for i = Qselect
            Rest(g,d) = D50(d)*ust(g,i)/v; %[-] Particle Reynolds number (vR84,Leeuw20), Reynolds shear velocity number (Klein05), Boundary Reynolds number (Miedema10)
            shield(g,d) = tau(g,i)/((rhos-rho)*grav(g)*D50(d)); %[-] Shields parameter/Particle mobility parameter
            lambda(g,d) = ust(g,i)/ws(g,d); %[-] Movability number (Liu 1958)
            AdvL(g,d) = u(g,i)*h(g,i)/ws(g,d); %[m] Advection length
        end
    end
end

%% Critical threshold of motion
acc = 1e-5;
for d = 1:LD
    for g = 1:Lg
        %Zanke 2003 iterated
        shieldini = shield(g,d);
        minimise = 1;
        while max(minimise)>acc
            tau_ini = shieldini *((rhos-rho)*grav(g)*D50(d));
            ust_ini = (tau_ini/rho)^0.5;
            Rest_ini(g,d) = D50(d)*ust_ini/v;
            uprmsb_ust = 0.31*Rest_ini(g,d)*exp(-0.1*Rest_ini(g,d))+1.8*exp(-0.88*D50(d)/(shieldini*R*D50(d)/sin(S)))*(1-exp(-0.1*Rest_ini(g,d)));
            Pt = 1 - exp(-0.08*Rest_ini(g,d));
            B = (1-Pt)*(2.5.*log(Rest_ini(g,d))+5.25)+8.5*Pt;
            uy_ust = ((1-Pt)/Rest_ini(g,d)^2+Pt/(2.5*log(1)+B^2))^-0.5;
            ub_ust = 0.8+0.9*uy_ust;
            uprmsb_ub = uprmsb_ust / ub_ust;
            K = 1 + 3e-8/((rhos-rho)*D50(d)^2);
            shield1 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
            minimise = abs(shield1 - shieldini);
            shieldini = shield1;
        end
        shield_cr(g,d) = shieldini;
        tau_cr(g,d) = shield_cr(g,d) * (rhos-rho)*grav(g)*D50(d);
        ust_cr(g,d) = (tau_cr(g,d) / rho)^0.5;
        lambda_cr(g,d) = ust_cr(g,d) / ws(g,d);       
    end
end
clear tau_Q_ini ust_Q_ini Rest_ini uprmsb_ust Pt B uy_ust ub_ust uprmsb_ub shield1 minimise shieldini K

%% Bedload transport
for i = Qselect
for g = 1:Lg
    for d = 1:LD
        S0(g,d) = max(((shield(g,d)-shield_cr(g,d))/shield_cr(g,d)),0);
        T0(g,d) = max(((ust(g,i)^2-ust_cr(g,d)^2)/ust_cr(g,d)^2),0);

        %A(theta-theta_cr)^B
        Qb_MPM(g,d) = 8*(max(shield(g,d)-shield_cr(g,d),0))^1.5;        	%[-] Meyer-Peter and Muller 1948 (also in D3D and Carrillo 2021) | theta_cr=0.047
        qb_MPM(g,d) = Qb_MPM(g,d)*(grav(g)*R*D50(d)^3)^0.5;                 %[m3/ms]
        Qb_MP(g,d) = (4*max(shield(g,d)-0.188,0))^1.5;                     	%[-] Meyer-Peter 1949/1951 as in Carrillo 2021, 1.25<s<4.2
        qb_MP(g,d) = Qb_MP(g,d)*(grav(g)*R*D50(d)^3)^0.5;                   %[m3/ms]
        Qb_E50(g,d) = 3.97*(max(shield(g,d)-shield_cr(g,d),0))^1.5;         %[-] Einstein 1950 as in de Leeuw 2020 (I have the paper, but am unable to rewrite)
        qb_E50(g,d) = Qb_E50(g,d)*(grav(g)*R*D50(d)^3)^0.5;                 %[m3/ms]
        Qb_FB(g,d) = 5.7*(max(shield(g,d)-shield_cr(g,d),0))^1.5;           %[-] Fernandez Luque and van Beek 1976 as in de Leeuw 2020 (payed download)
        qb_FB(g,d) = Qb_FB(g,d)*(grav(g)*R*D50(d)^3)^0.5;                   %[m3/ms]
        Qb_R(g,d) = 11*(max(shield(g,d)-shield_cr(g,d),0))^1.65;        	%[-] Ribberink 1998
        qb_R(g,d) = Qb_R(g,d)*(grav(g)*R*D50(d)^3)^0.5;                     %[m3/ms]
        Qb_HJ(g,d) = 5*(max(shield(g,d)-shield_cr(g,d),0))^1.5;             %[-] Hunziker and Jaeggi 2002 as in Wong and Parker 2006 (payed download) | theta_cr=0.05
        qb_HJ(g,d) = Qb_HJ(g,d)*(grav(g)*R*D50(d)^3)^0.5;                   %[m3/ms]
        Qb_WP1(g,d) = 4.93*(max(shield(g,d)-shield_cr(g,d),0))^1.6;         %[-] Wong and Parker 2006, s=2.55 | theta_cr=0.047
        Qb_WP2(g,d) = 3.97*(max(shield(g,d)-shield_cr(g,d),0))^1.5;         %[-] Wong and Parker 2006, s=2.55 | theta_cr=0.0495
        qb_WP1(g,d) = Qb_WP1(g,d)*(grav(g)*R*D50(d)^3)^0.5;                 %[m3/ms]
        qb_WP2(g,d) = Qb_WP2(g,d)*(grav(g)*R*D50(d)^3)^0.5;                 %[m3/ms]
        
        % A theta^1.5 exp(B theta_cr/theta)
        Qb_Ca(g,d) = 12*shield(g,d)^1.5*exp(-4.5*(shield_cr(g,d)/shield(g,d))); %[-] Camenen 2005
        qb_Ca(g,d) = Qb_Ca(g,d)*(grav(g)*R*D50(d)^3)^0.5;                   %[m3/ms]
        Qb_C(g,d) = 13*shield(g,d)^1.5*exp(-shield_cr(g,d)/shield(g,d)^1.5);          %[-] Cheng 2002 (also in Carrillo 2021), 2.53<s<2.69, 0.73<S<1.2%, 0.068<D<1.27ft, 0.093<Q<1.118ft3/s | theta_cr=0.05
        qb_C(g,d) = Qb_C(g,d)*(grav(g)*R*D50(d)^3)^0.5;                     %[m3/ms]
        
        %A theta^0.5*(theta-theta_cr)
        Qb_N(g,d) = 12*shield(g,d)^0.5*(max(shield(g,d)-shield_cr(g,d),0)); %[-] Nielsen 1992
        qb_N(g,d) = Qb_N(g,d)*(grav(g)*R*D50(d)^3)^0.5;                     %[m3/ms]
        %A(f) theta^0.5*(theta-theta_cr)
        Qb_S(g,d) = 4.2*S^0.6*(u(g,i)/ust(g,i))*shield(g,d)^0.5*(max(shield(g,d)-shield_cr(g,d),0)); %[-] Smart 1984
        qb_S(g,d) = Qb_S(g,d)*(grav(g)*R*D50(d)^3)^0.5;                     %[m3/ms]
        %A (theta^0.5-theta_cr^0.5)*(theta-theta_cr)
        Qb_AM(g,d) = 17*(max(shield(g,d)-shield_cr(g,d),0))*(max(shield(g,d)^0.5-shield_cr(g,d)^0.5,0)); %[-] Ashida and Michiue 1972 as in D3D and Carrillo 2021 (I have the paper, but it is in Japanese) | theta_cr=0.05
        qb_AM(g,d) = Qb_AM(g,d)*(grav(g)*R*D50(d)^3)^0.5;                   %[m3/ms]
        
        %other
        p_EF(g,d) = min((1+((pi/6)*1/(max(shield(g,d)-shield_cr(g,d),0)))^4)^(-0.25),1); %[-] probability | beta=1 Engelund and Fredsoe 1982 as in Garcia and Parker 1991
        Qb_EF(g,d) = 5*p_EF(g,d)*(max(shield(g,d)^0.5-0.7*shield_cr(g,d)^0.5,0)); %[-] Engelund and Fredsoe 1976 (9.3*(pi/6)=~5) | theta_cr=0.05,0.06
        qb_EF(g,d) = Qb_EF(g,d)*(grav(g)*R*D50(d)^3)^0.5;                   %[m3/ms]
        Qb_P(g,d) = 11.2*(max(shield(g,d)-shield_cr(g,d),0))^4.5/shield(g,d)^3; %[-] Parker 1979 as in Carrillo 2021 or Parker 1982 as in Kleinhans 2005 (payed download) | theta_cr=0.03                Qb_P(g,d) = 11.2*(shield(g,d)-shield_cr_nc(g,d))^4.5/shield(g,d)^3; %[-] Parker 1979 as in Carrillo 2021 or Parker 1982 as in Kleinhans 2005 (payed download) | theta_cr=0.03
        qb_P(g,d) = Qb_P(g,d)*(grav(g)*R*D50(d)^3)^0.5;                     %[m3/ms]
        Qb_vR(g,d) = 0.053*Dst(g,d)^-0.3*T0(g,d)^2.1;                       %[-] van Rijn 1984, 200<D<2000mum
        qb_vR(g,d) = Qb_vR(g,d)*(grav(g)*R*D50(d)^3)^0.5;                   %[m3/ms]
        Qb_vR2(g,d) = 0.1*Dst(g,d)^(-0.3)*S0(g,d)^1.5;                      %[-] van Rijn 1984 as in Kleinhans 2005
        qb_vR2(g,d) = Qb_vR2(g,d)*(grav(g)*R*D50(d)^3)^0.5;                 %[m3/ms]
        
        %no theta_cr
        qb_EH(g,d) = 0.05*u(g,i)^5/(grav(g)^0.5*C(g,i)^3*R^2*D50(d));       %[m3/ms]
        Qb_EH(g,d) = qb_EH(g,d)/(grav(g)*R*D50(d)^3)^0.5;                   %[-] Engelund and Hansen as in D3D
        Qb_W(g,d) = 12*shield(g,d)^1.5;                                     %[-] Wilson 1966 as in Soulsby 2005
        qb_W(g,d) = Qb_W(g,d)*(grav(g)*R*D50(d)^3)^0.5;                     %[m3/ms]
        Qb_E42(g,d) = 2.1*exp(-0.391/shield(g,d));                          %[-] Einstein 1942 as in Carrillo 2021, 1.25<s<4.25, 0.315<D<28.6mm, Qb<0.4
        qb_E42(g,d) = Qb_E42(g,d)*(grav(g)*R*D50(d)^3)^0.5;                 %[m3/ms]
        
        G2(g,d) = rhos*D50(d)^2*ust(g,i)^2/(rho*14*v^2);                    %[?] Grain flow number (Bagnold 1966), 150<G2<6000:tan(phi)=-0.236, G2<150:tan(phi)=0.7, G2>6000:tan(phi)=0.374
        if D50(d)>0.000015 && D50(d)<0.00006
            eb(g,d) = -0.012*log(3.28*u(g,i))+0.15;
        elseif D50(d)>=0.00006 && D50(d)<0.0002
            eb(g,d) = -0.013*log(3.28*u(g,i))+0.145;
        elseif D50(d)>=0.0002 && D50(d)<0.0007
            eb(g,d) = -0.016*log(3.28*u(g,i))+0.139;
        elseif D50(d)>=0.0007
            eb(g,d) = -0.028*log(3.28*u(g,i))+0.135;
        else
            eb(g,d) = NaN;
        end
        if G2(g,d)>150 && G2(g,d)<6000
            qb_B(g,d) = eb(g,d)*u(g,i)*tau(g,i)/((rhos-rho)*grav(g)*cos(S*(-0.236-tan(S))));
        elseif G2(g,d)<=150
            qb_B(g,d) = eb(g,d)*u(g,i)*tau(g,i)/((rhos-rho)*grav(g)*cos(S*(0.7-tan(S))));
        elseif G2(g,d)>=6000
            qb_B(g,d) = eb(g,d)*u(g,i)*tau(g,i)/((rhos-rho)*grav(g)*cos(S*(0.374-tan(S))));
        else
            qb_B(g,d) = NaN;
        end
        Qb_B(g,d) = qb_B(g,d)/(grav(g)*R*D50(d)^3)^0.5;                     %[-] Bagnold 1966 as in Kleinhans 2005
        
    end
end
end

%% Suspended transport
for g=1:Lg
    for d=1:LD
        %reference height
        a_E50(d) = 2*D50(d);
        a_EF(d) = 2*D50(d);
        a_IK(g) = 0.05*h(g);
        a_CR(g) = 0.05*h(g);
        a_AF(g) = 0.05*h(g);
        a_GP(g) = 0.05*h(g);
        a_WP(g) = 0.05*h(g);
        a_dL(g) = 0.1*h(g);
        a_vR(g) = max(0.01*h(g),ks);
        a_SM(g,d) = 26.3*(max(shield(g,d)-shield_cr(g,d),0))*D50(d)+ks;
        a_ML(g,d) = 0.68*(tau(g)/tau_cr(g,d))*D50(d)/(1+(0.0204*log(D50(d)*100)^2+0.022*log(D50(d)*100)+0.0709)*(tau(g)/tau_cr(g,d)));
        
        %Rouse number
        Rouse(g,d) = ws(g,d)/(Kappa*ust(g));
        Rouse_dL(g,d) = (ust(g)/ws(g,d))^-0.45;
        
        %Ca / Es
        Ca_E50(g,d) = (1/32.2)*Qb_E50(g,d)/shield(g,d)^0.5;
        labda_EF(g,d) = (max((shield(g,d)-shield_cr(g,d)-(1*p_EF(g,d)*pi/6))/(0.027*(R+1)*shield(g,d)),0))^0.5;
        Ca_EF(g,d) = 0.65/(1+labda_EF(g,d)^-1)^3;
        Ca_SM(g,d) = 0.65*0.0024*S0(g,d)/(1+0.0024*S0(g,d));
        Ca_ML(g,d) = 0.065*0.004*S0(g,d)/(1+0.004*S0(g,d));
        Ca_vR(g,d) = 0.015*(D50(d)/a_vR(g))*S0(g,d)^1.5/Dst(g,d)^0.3;
        Z_AF(g,d) = min(Rep(g,d)^0.5*ust(g)/ws(g,d),13.2); %max susp at 13.2, no susp<5
        Ca_AF(g,d) = 0.000000000003*Z_AF(g,d)^10*(1-(5/Z_AF(g,d)));
        Z_GP(g,d) = Rep(g,d)^0.6*ust(g)/ws(g,d);
        Ca_GP(g,d) = 0.00000013*Z_GP(g,d)^5/(1+(0.00000013/0.3)*Z_GP(g,d)^5);
        Z_WP(g,d) = Rep(g,d)^0.6*ust(g)/ws(g,d);
        Ca_WP(g,d) = 0.00000078*Z_WP(g,d)^5/(1+(0.00000078/0.3)*Z_WP(g,d)^5);
        Ca_dL(g,d) = 0.000474*(ust(g)/ws(g,d))^1.77*Fr(g)^1.18;
        fun = @(eta) (((1-eta)./eta).*(0.05/(1-0.05))).^Rouse(g,d);
        I_CR(g,d) = integral(fun,0.05,1);
        Cm_CR(g,d) = 0.034*(1-(ks/h(g))^0.06)*ust(g)^2*u(g)/(grav(g)*R*h(g)*ws(g,d));
        Ca_CR(g,d) = 1.13*Cm_CR(g,d)/I_CR(g,d);
        A_IK(g,d) = 0.143/shield(g,d)-2;
        fun2 = @(eta2) exp(-1.*eta2.^2);
        Omega_IK(g,d) = (shield(g,d)/0.143)*(2+((exp(-1*A_IK(g,d)^2))/(integral(fun2,A_IK(g,d),Inf))))-1;
        Ca_IK(g,d) = max(0.008*(0.14*ust(g)*Omega_IK(g,d)/(ws(g,d)*shield(g,d))-1),0); %goes negative without max(CA_IK,0)...
        
        %Profile
        z(g,:) = 0:0.01:h(g); %[m] height, vector over water depth for Rouse profile
        Lz = size(z,2);
        for j=1:Lz
            Cz_E50(g,d,j) = Ca_E50(g,d) * ( ((h(g)-z(g,j))/z(g,j)) * (a_E50(d)/(h(g)-a_E50(d))) )^Rouse(g,d);
            Cz_EF(g,d,j) = Ca_EF(g,d) * ( ((h(g)-z(g,j))/z(g,j)) * (a_EF(d)/(h(g)-a_EF(d))) )^Rouse(g,d);
            Cz_IK(g,d,j) = Ca_IK(g,d) * ( ((h(g)-z(g,j))/z(g,j)) * (a_IK(g)/(h(g)-a_IK(g))) )^Rouse(g,d);
            Cz_CR(g,d,j) = Ca_CR(g,d) * ( ((h(g)-z(g,j))/z(g,j)) * (a_CR(g)/(h(g)-a_CR(g))) )^Rouse(g,d);
            Cz_AF(g,d,j) = Ca_AF(g,d) * ( ((h(g)-z(g,j))/z(g,j)) * (a_AF(g)/(h(g)-a_AF(g))) )^Rouse(g,d);
            Cz_GP(g,d,j) = Ca_GP(g,d) * ( ((h(g)-z(g,j))/z(g,j)) * (a_GP(g)/(h(g)-a_GP(g))) )^Rouse(g,d);
            Cz_WP(g,d,j) = Ca_WP(g,d) * ( ((h(g)-z(g,j))/z(g,j)) * (a_WP(g)/(h(g)-a_WP(g))) )^Rouse(g,d);
            Cz_dL(g,d,j) = Ca_dL(g,d) * ( ((h(g)-z(g,j))/z(g,j)) * (a_dL(g)/(h(g)-a_dL(g))) )^Rouse(g,d);
            Cz_vR(g,d,j) = Ca_vR(g,d) * ( ((h(g)-z(g,j))/z(g,j)) * (a_vR(g)/(h(g)-a_vR(g))) )^Rouse(g,d);
            Cz_SM(g,d,j) = Ca_SM(g,d) * ( ((h(g)-z(g,j))/z(g,j)) * (a_SM(g,d)/(h(g)-a_SM(g,d))) )^Rouse(g,d);
            Cz_ML(g,d,j) = Ca_ML(g,d) * ( ((h(g)-z(g,j))/z(g,j)) * (a_ML(g,d)/(h(g)-a_ML(g,d))) )^Rouse(g,d);
        end
        clear z
        
        %Integral
        fun = @(z) Ca_E50(g,d) * ( ((h(g)-z)./z) * (a_E50(d)/(h(g)-a_E50(d))) ).^Rouse(g,d);
        qs_E50(g,d) = integral(fun,a_E50(d),h(g));                          %Einstein 1950
        fun = @(z) Ca_EF(g,d) * ( ((h(g)-z)./z) * (a_EF(d)/(h(g)-a_EF(d))) ).^Rouse(g,d);
        qs_EF(g,d) = integral(fun,a_EF(d),h(g));                            %Engelund and Fredsoe 1976
        fun = @(z) Ca_IK(g,d) * ( ((h(g)-z)./z) * (a_IK(g)/(h(g)-a_IK(g))) ).^Rouse(g,d);
        qs_IK(g,d) = integral(fun,a_IK(g),h(g));                            %Itakura and Kishi 1980 as in Garcia and Parker 1991
        fun = @(z) Ca_CR(g,d) * ( ((h(g)-z)./z) * (a_CR(g)/(h(g)-a_CR(g))) ).^Rouse(g,d);
        qs_CR(g,d) = integral(fun,a_CR(g),h(g));                            %Celik and Rodi 1984 as in Garcia and Parker 1991
        fun = @(z) Ca_AF(g,d) * ( ((h(g)-z)./z) * (a_AF(g)/(h(g)-a_AF(g))) ).^Rouse(g,d);
        qs_AF(g,d) = integral(fun,a_AF(g),h(g));                            %Akiyama and Fukushima 1986 as in Garcia and Parker 1991
        fun = @(z) Ca_GP(g,d) * ( ((h(g)-z)./z) * (a_GP(g)/(h(g)-a_GP(g))) ).^Rouse(g,d);
        qs_GP(g,d) = integral(fun,a_GP(g),h(g));                            %Garcia and Parker 1991
        fun = @(z) Ca_WP(g,d) * ( ((h(g)-z)./z) * (a_WP(g)/(h(g)-a_WP(g))) ).^Rouse(g,d);
        qs_WP(g,d) = integral(fun,a_WP(g),h(g));                            %Wright and Parker 2004
        fun = @(z) Ca_dL(g,d) * ( ((h(g)-z)./z) * (a_dL(g)/(h(g)-a_dL(g))) ).^Rouse(g,d);
        qs_dL(g,d) = integral(fun,a_dL(g),h(g));                            %de Leeuw 2020
        fun = @(z) Ca_vR(g,d) * ( ((h(g)-z)./z) * (a_vR(g)/(h(g)-a_vR(g))) ).^Rouse(g,d);
        qs_vR(g,d) = integral(fun,a_vR(g),h(g));                            %van Rijn 1984
        fun = @(z) Ca_SM(g,d) * ( ((h(g)-z)./z) * (a_SM(g,d)/(h(g)-a_SM(g,d))) ).^Rouse(g,d);
        qs_SM(g,d) = integral(fun,a_SM(g,d),h(g));                          %Smith and McLean 1977
        fun = @(z) Ca_ML(g,d) * ( ((h(g)-z)./z) * (a_ML(g,d)/(h(g)-a_ML(g,d))) ).^Rouse(g,d);
        qs_ML(g,d) = integral(fun,a_ML(g,d),h(g));                          %McLean 1992
    end
end

%% Ratios
for a = 1
qt_vR = qb_vR + qs_vR;
bedload_percent_vR = qb_vR./qt_vR*100;
suspended_percent_vR = qs_vR./qt_vR*100;
ratio_vR = qb_vR./qs_vR;

qt_E50 = qb_E50 + qs_E50;
bedload_percent_E50 = qb_E50./qt_E50*100;
suspended_percent_E50 = qs_E50./qt_E50*100;
ratio_E50 = qb_E50./qs_E50;

qt_EF = qb_EF + qs_EF;
bedload_percent_EF = qb_EF./qt_EF*100;
suspended_percent_EF = qs_EF./qt_EF*100;
ratio_EF = qb_EF./qs_EF;

qt_dL = qb_FB + qs_dL;
bedload_percent_dL = qb_FB./qt_dL*100;
suspended_percent_dL = qs_dL./qt_dL*100;
ratio_dL = qb_FB./qs_dL;
end

%% Testplot Rouse profiles
for a = 1
    % i = 10;
    % j = 2;
    % figure()
    % hold on
    % plot(squeeze(Cz_E50(i,j,:)),z)
    % plot(squeeze(Cz_EF(i,j,:)),z)
    % plot(squeeze(Cz_CR(i,j,:)),z)
    % plot(squeeze(Cz_AF(i,j,:)),z)
    % plot(squeeze(Cz_GP(i,j,:)),z)
    % plot(squeeze(Cz_WP(i,j,:)),z)
    % plot(squeeze(Cz_dL(i,j,:)),z)
    % plot(squeeze(Cz_vR(i,j,:)),z)
    % plot(squeeze(Cz_SM(i,j,:)),z)
    % plot(squeeze(Cz_ML(i,j,:)),z)
    % set(gca,'Xscale','log');
end

%% Plot settings
for a = 1
    pl.width = 18;
    pl.height = 18;
    pl.line = 1;
    pl.line_ax = 0.75;
    pl.fsz = 7;
    pl.fsz2 = 8.5;
    CMAPOBJ = clrmap('read','Gravity.clrmap');
    clr1 = clrmap(CMAPOBJ,Lg);
    CMAPOBJ = clrmap('read','Discharge.clrmap');
    %clr2 = clrmap(CMAPOBJ,Lh);
    clr3 = clrmap(CMAPOBJ,LQ);
    clrM = [1 0.804 0.804];
    clrE = [0.8 0.8 0.9804];
    pl.vertical = 3;
    pl.horizantal = 2;
    pl.lmarge = 0.065;
    pl.rmarge = 0.025;
    pl.bmarge = 0.055;
    pl.tmarge = 0.01;
    pl.intv = 0.065;
    pl.inth = 0.065;
    pl.wt = (1 - pl.lmarge - pl.rmarge - ((pl.horizantal-1)*pl.inth)) / pl.horizantal;
    pl.ht = (1 - pl.bmarge - pl.tmarge - ((pl.vertical-1)*pl.intv)) / pl.vertical;
    y = 0.94;
    y2 = 0.97;
    D50L = D50(1);
    D50R = D50(end);
    xd = 10^(log10(D50L)+( log10(D50R)-log10(D50L) )*0.01);
end

%% Plot
for a = 1
    close all; clc;
    f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[130 0 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');
    
    s1(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(D50,qb_MPM(g,:),'-','color',clr1(g,:),'linewidth',pl.line); 
    end
    for g = 1:Lg
        plot(D50,qb_E42(g,:),':','color',[clr1(g,:) 0.2],'linewidth',pl.line); %-no crit, -down small grains
        plot(D50,qb_MP(g,:),':','color',[clr1(g,:) 0.2],'linewidth',pl.line); %-no crit, critical cutoff without theta_cr
        plot(D50,qb_E50(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,qb_B(g,:),':','color',[clr1(g,:) 0.2],'linewidth',pl.line); %-no crit, low straight line
        plot(D50,qb_W(g,:),':','color',[clr1(g,:) 0.2],'linewidth',pl.line); %-no crit, high straight line
        plot(D50,qb_AM(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,qb_FB(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,qb_EF(g,:),'-','color',clr1(g,:),'linewidth',pl.line); %-down small grains
        plot(D50,qb_P(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,qb_S(g,:),'-','color',clr1(g,:),'linewidth',pl.line); %1 orde lager dan andere
        plot(D50,qb_vR(g,:),'-','color',clr1(g,:),'linewidth',pl.line); %Rare vorm
        plot(D50,qb_vR2(g,:),'-','color',clr1(g,:),'linewidth',pl.line); %Rare vorm, lager dan vR
        plot(D50,qb_N(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,qb_R(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,qb_HJ(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,qb_C(g,:),'-','color',clr1(g,:),'linewidth',pl.line); %no critical cutoff
        plot(D50,qb_Ca(g,:),'-','color',clr1(g,:),'linewidth',pl.line); %no critical cutoff
        plot(D50,qb_WP1(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,qb_WP2(g,:),'-','color',clr1(g,:),'linewidth',pl.line);  
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    %xlim([D50(1) D50(end)]);
    xlim([D50L D50R]);
    ylim([10^-6.5 10^1]);
    l(1) = ylabel('Bedload transport (q_b) [m^3/ms]');
    xlabel('Grain size (D_{50}) [m]');
    t(1) = text(xd,10^(-6.5+(1--6.5)*y2),' a');

    s1(2) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(D50,qs_dL(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
    end
    for g = 1:Lg    
        plot(D50,qs_IK(g,:),':','color',[clr1(g,:) 0.2],'linewidth',pl.line); %small range
        plot(D50,qs_EF(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,qs_E50(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,qs_vR(g,:),'-','color',clr1(g,:),'linewidth',pl.line); 
        plot(D50,qs_CR(g,:),':','color',[clr1(g,:) 0.2],'linewidth',pl.line); %too much transport big grain sizes
        plot(D50,qs_AF(g,:),':','color',[clr1(g,:) 0.2],'linewidth',pl.line); %too much transport big grain sizes
        plot(D50,qs_GP(g,:),':','color',[clr1(g,:) 0.2],'linewidth',pl.line); %too much transport big grain sizes
        plot(D50,qs_WP(g,:),':','color',[clr1(g,:) 0.2],'linewidth',pl.line); %too much transport big grain sizes
        plot(D50,qs_SM(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
        plot(D50,qs_ML(g,:),'-','color',clr1(g,:),'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    %xlim([D50(1) D50(end)]);
    xlim([D50L D50R]);
    ylim([10^-6.5 10^1]);
    l(2) = ylabel('Suspended transport (q_s) [m^3/ms]');
    xlabel('Grain size (D_{50}) [m]');
    t(2) = text(xd,10^(-6.5+(1--6.5)*y2),' b');

    s1(3) = axes('Position',[pl.lmarge 1-pl.tmarge-2*pl.ht-pl.intv pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(D50,qt_E50(g,:),'-','color',clr1(g,:),'linewidth',pl.line)
        plot(D50,qt_EF(g,:),'--','color',clr1(g,:),'linewidth',pl.line)   
        plot(D50,qt_dL(g,:),'-.','color',clr1(g,:),'linewidth',pl.line)   
        plot(D50,qt_vR(g,:),':','color',clr1(g,:),'linewidth',pl.line)
        plot(D50,qb_EH(g,:),'-','color',[clr1(g,:) 0.2],'linewidth',pl.line)   
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    %xlim([D50(1) D50(end)]);
    xlim([D50L D50R]);
    ylim([10^-4 10^1]);
    l(3) = ylabel('Total transport (q_t) [m^3/ms]');
    xlabel('Grain size (D_{50}) [m]');
    t(3) = text(xd,10^(-4+(1--4)*y2),' c');
    
    s1(4) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-2*pl.ht-pl.intv pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(D50,suspended_percent_E50(g,:),'-','color',clr1(g,:),'linewidth',pl.line)
        plot(D50,suspended_percent_EF(g,:),'--','color',clr1(g,:),'linewidth',pl.line)
        plot(D50,suspended_percent_dL(g,:),'-.','color',clr1(g,:),'linewidth',pl.line)
        plot(D50,suspended_percent_vR(g,:),':','color',clr1(g,:),'linewidth',pl.line)
    end
    set(gca,'Xscale','log');
    %xlim([D50(1) D50(end)]);
    xlim([D50L D50R]);
    ylim([0 100]);
    l(3) = ylabel('Suspended transport of total [%]');
    xlabel('Grain size (D_{50}) [m]');
    t(3) = text(xd,100*y2,' d');
    
    s1(5) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    plot(D50,qt_E50(1,:)./qt_E50(2,:),'-','color',[0 0 0],'linewidth',pl.line)
    plot(D50,qt_EF(1,:)./qt_EF(2,:),'--','color',[0 0 0],'linewidth',pl.line)   
    plot(D50,qt_dL(1,:)./qt_dL(2,:),'-.','color',[0 0 0],'linewidth',pl.line)   
    plot(D50,qt_vR(1,:)./qt_vR(2,:),':','color',[0 0 0],'linewidth',pl.line)
    plot(D50,qb_EH(1,:)./qb_EH(2,:),'-','color',[0 0 0 0.2],'linewidth',pl.line)
    plot([D50(1) D50(end)], [1 1],'-g')
    set(gca,'Xscale','log');
    %xlim([D50(1) D50(end)]);
    xlim([D50L D50R]);
    ylim([0 10]);
    l(3) = ylabel('Mars/Earth ratio total transport [-]');
    xlabel('Grain size (D_{50}) [m]');
    t(3) = text(xd,10*y2,' e');

    leg = legend(s1(1),['Mars Q=' num2str(Q(Qselect)) ' m^3/s'],['Earth Q=' num2str(Q(Qselect)) ' m^3/s'],['Unsuitable predictor']);
    set(leg,'location','northeast');
    leg2 = legend(s1(2),['Mars Q=' num2str(Q(Qselect)) ' m^3/s'],['Earth Q=' num2str(Q(Qselect)) ' m^3/s'],['Unsuitable predictor']);
    set(leg2,'location','northeast');
    leg3 = legend(s1(3),'Einstein 1950','Engelund and Fredsoe 1976','de Leeuw 2020','van Rijn 1984','Engelund and Hansen 1967','location','northeast');
    leg4 = legend(s1(4),'Einstein 1950','Engelund and Fredsoe 1976','de Leeuw 2020','van Rijn 1984','location','northeast');
    leg5 = legend(s1(5),'Einstein 1950','Engelund and Fredsoe 1976','de Leeuw 2020','van Rijn 1984','Engelund and Hansen 1967','location','northeast');

    xtick = logspace(-6,1,8);
    set(s1,'box','on','Layer','top', ...
        'XMinorTick','on','YMinorTick','on', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5], ...
        'XTick',xtick(1:1:end));
    
    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);
    
    set(gcf,'renderer','painters');
    print(f1,'-dpng',[output '/FlowOnMars4_Transport_DQ'],'-r400');
    print(f1,'-dpdf',[output '/FlowOnMars4_Transport_DQ'],'-r400');   
end
