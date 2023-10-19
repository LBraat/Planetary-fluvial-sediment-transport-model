%% Flow on Mars 4 - Transport Q (1 Q equal for Mars and Earth)
% Author: Lisanne Braat (lisannebraat@gmail.com)
% Last update: 2023-08-01

%% Initialize
clear variables
close all
clc

output = 'FlowOnMars_exportfig_v2';
addpath(genpath('Checkout'));

%% Input parameters
for a = 1
    h = 2:0.5:15; %[m] Water depth
    hselect = 7; %[1 27];
    W = 200; %[m] Channel width
    grav = [3.7 9.8]; %[m/s2] Gravitational acceleration
    S = 0.001; %[m/m] Channel slope
    TC = 4; %[degrees C] Temperature
    rho = 1000; %[kg/m3] Water density

    % Sediment parameters
    rhos = 2900; %[kg/m3] Sediment density
    D = logspace(log10(63*10^-6),-1,50); %[m] Grain size vector
    CDryB = 1750; %[kg/m3] Dry bed density
    C1 = 20; %18-24, smooth to rough sphere, affects small diameters, Ferguson and Church 2004
    C2 = 1; %0.4-1.2, smooth to rough sphere, affects big diameters, Ferguson and Church 2004
    phi = 30; %angle of repose in degrees
    Kappa = 0.4; %Kappa for Rouse

    % Derived input parameters
    TK = TC+273.15; %[degrees K] Temperature
    mu = 1.856e-14*exp(4209/TK+0.04527*TK-3.376e-5*TK^2); %[Pa | Ns/m2] Dynamic viscosity
    v = mu/rho; %[m2/s] Viscosity
    ks = 2.5*D; %[m] Nikurandse roughness height

    R = (rhos-rho)/rho; %[-] Relative density
    n = 1-CDryB/rhos; %[-] Porosity

    Lh = length(h);
    Lg = length(grav);
    LD = length(D);
end
clear TC TK mu CDryB

%% allocate
for a = 1
    Rw = NaN(1,Lh);
    f = NaN(Lg,Lh,LD);
    u = NaN(Lg,Lh,LD);
    Q = NaN(Lg,Lh,LD);
    Fr = NaN(Lg,Lh,LD);
    Re = NaN(Lg,Lh,LD);
    tau = NaN(Lg,Lh);
    ust = NaN(Lg,Lh);

    ws = NaN(Lg,LD);
    Rep = NaN(Lg,LD);
    Dst = NaN(Lg,LD);
    Rest = NaN(Lg,Lh,LD);
    shield = NaN(Lg,Lh,LD);
    lambda = NaN(Lg,Lh,LD);
    AdvL = NaN(Lg,Lh,LD);

    shield_cr = NaN(Lg,LD);
    tau_cr = NaN(Lg,LD);
    ust_cr = NaN(Lg,LD);
    lambda_cr = NaN(Lg,LD);
    S0 = NaN(Lg,LD);
    T0 = NaN(Lg,LD);
    Rest_ini = NaN(Lg,LD);

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
    qt_EH = NaN(Lg,LD);
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

    ub = NaN(Lg,LD);
    hb = NaN(Lg,LD);
    Cb = NaN(Lg,LD);
    Ca_dL2 = NaN(Lg,LD);
    qs_dL2_a = NaN(Lg,LD);
    qs_dL2_b = NaN(Lg,LD);
    qs_dL2 = NaN(Lg,LD);
end

%% Hydro parameters h_mars = h_earth
for i = 1:Lh
    Rw(i) = (h(i)*W)/(h(i)+h(i)+W); %[m] Hydraulic radius
    for g = 1:Lg
        for d = 1:LD
            f(g,i,d) = 8/(5.74*log10(12.2*h(i)/ks(d)))^2; %[-] friction factor
            u(g,i,d) = (8*grav(g)*Rw(i)*S/f(g,i,d))^0.5; %[m/s] Velocity
            Q(g,i,d) = h(i)*W*u(g,i,d); %[m3/s] Discharge
            Fr(g,i,d) = u(g,i,d)/(grav(g)*h(i))^0.5; %[-] Froude number
            Re(g,i,d) = u(g,i,d)*h(i)/v; %[-] Reynolds number
            tau(g,i) = rho*grav(g)*Rw(i)*S; %[N/m2] Bed shear stress
            ust(g,i) = (tau(g,i)/rho)^0.5; %[m/s] Shear velocity
        end
    end
end

%% Sediment parameters
for d = 1:LD
    for g = 1:Lg
        ws(g,d) = R*grav(g)*D(d)^2/(C1*v+(0.75*C2*R*grav(g)*D(d)^3)^0.5); %[m/s] Settling velocity (Ferguson and Church 2004)
        Rep(g,d) = D(d)^1.5*(R*grav(g))^0.5/v; %[-] Particle Reynolds number (Klein05,Leeuw20,Miedema), >70 turbulent <3.5 laminar
        Dst(g,d) = D(d)*(R*grav(g)/v^2)^(1/3); %[-] Dimensionless particle parameter/Bonnefille number

        % h_mars = h_earth
        for i = 1:Lh
            Rest(g,i,d) = D(d)*ust(g,i)/v; %[-] Particle Reynolds number (vR84,Leeuw20,nino03), Reynolds shear velocity number (Klein05), Boundary Reynolds number (Miedema10), fluid reynolds number (Bagnold1966)
            shield(g,i,d) = tau(g,i)/((rhos-rho)*grav(g)*D(d)); %[-] Shields parameter/Particle mobility parameter
            lambda(g,i,d) = ust(g,i)/ws(g,d); %[-] Movability number (Liu 1958)
            AdvL(g,i,d) = u(g,i,d)*h(i)/ws(g,d); %[m] Advection length
        end
    end
end
clear C1 C2

%% Critical threshold of motion h_mars = h_earth
i = hselect;
acc = 1e-5;
for d = 1:LD
    for g = 1:Lg
        %Zanke 2003 iterated
        shieldini = shield(g,i,d);
        minimise = 1;
        while max(minimise)>acc
            tau_Q_ini = shieldini *((rhos-rho)*grav(g)*D(d));
            ust_Q_ini = (tau_Q_ini/rho)^0.5;
            Rest_ini(g,d) = D(d)*ust_Q_ini/v;
            uprmsb_ust = 0.31*Rest_ini(g,d)*exp(-0.1*Rest_ini(g,d))+1.8*exp(-0.88*D(d)/(shieldini*R*D(d)/sin(S)))*(1-exp(-0.1*Rest_ini(g,d)));
            Pt = 1 - exp(-0.08*Rest_ini(g,d));
            B = (1-Pt)*(2.5.*log(Rest_ini(g,d))+5.25)+8.5*Pt;
            uy_ust = ((1-Pt)/Rest_ini(g,d)^2+Pt/(2.5*log(1)+B^2))^-0.5;
            ub_ust = 0.8+0.9*uy_ust;
            uprmsb_ub = uprmsb_ust / ub_ust;
            K = 1 + 3e-8/((rhos-rho)*D(d)^2);
            shield1 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
            minimise = abs(shield1 - shieldini);
            shieldini = shield1;
        end
        shield_cr(g,d) = shieldini;
        tau_cr(g,d) = shield_cr(g,d) * (rhos-rho)*grav(g)*D(d);
        ust_cr(g,d) = (tau_cr(g,d) / rho)^0.5;
        lambda_cr(g,d) = ust_cr(g,d) / ws(g,d);
    end
end
clear i acc shieldini minimise tau_Q_ini ust_Q_ini Rest_ini uprmsb_ust Pt B uy_ust ub_ust uprmsb_ub K shield1 phi

%% Bedload transport
for i = hselect
    for g = 1:Lg
        for d = 1:LD
            S0(g,d) = max(((shield(g,i,d)-shield_cr(g,d))/shield_cr(g,d)),0);
            T0(g,d) = max(((ust(g,i)^2-ust_cr(g,d)^2)/ust_cr(g,d)^2),0);

            %A(theta-theta_cr)^B
            Qb_MPM(g,d) = 8*(max(shield(g,i,d)-shield_cr(g,d),0))^1.5; %[-] Meyer-Peter and Muller 1948 (also in D3D and Carrillo 2021) | theta_cr=0.047
            qb_MPM(g,d) = Qb_MPM(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            Qb_MP(g,d) = (4*max(shield(g,i,d)-0.188,0))^1.5; %[-] Meyer-Peter 1949/1951 as in Carrillo 2021, 1.25<s<4.2
            qb_MP(g,d) = Qb_MP(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            Qb_E50(g,d) = 3.97*(max(shield(g,i,d)-shield_cr(g,d),0))^1.5; %[-] Einstein 1950 as in de Leeuw 2020 (I have the paper, but am unable to rewrite)
            qb_E50(g,d) = Qb_E50(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            Qb_FB(g,d) = 5.7*(max(shield(g,i,d)-shield_cr(g,d),0))^1.5; %[-] Fernandez Luque and van Beek 1976 (also in de Leeuw 2020)
            qb_FB(g,d) = Qb_FB(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            Qb_R(g,d) = 11*(max(shield(g,i,d)-shield_cr(g,d),0))^1.65; %[-] Ribberink 1998
            qb_R(g,d) = Qb_R(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            Qb_HJ(g,d) = 5*(max(shield(g,i,d)-shield_cr(g,d),0))^1.5; %[-] Hunziker and Jaeggi 2002 as in Wong and Parker 2006 (payed download) | theta_cr=0.05
            qb_HJ(g,d) = Qb_HJ(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            Qb_WP1(g,d) = 4.93*(max(shield(g,i,d)-shield_cr(g,d),0))^1.6; %[-] Wong and Parker 2006, s=2.55 | theta_cr=0.047
            Qb_WP2(g,d) = 3.97*(max(shield(g,i,d)-shield_cr(g,d),0))^1.5; %[-] Wong and Parker 2006, s=2.55 | theta_cr=0.0495
            qb_WP1(g,d) = Qb_WP1(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            qb_WP2(g,d) = Qb_WP2(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]

            % A theta^1.5 exp(B theta_cr/theta)
            Qb_Ca(g,d) = 12*shield(g,i,d)^1.5*exp(-4.5*shield_cr(g,d)/shield(g,i,d)); %[-] Camenen 2005
            qb_Ca(g,d) = Qb_Ca(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            Qb_C(g,d) = 13*shield(g,i,d)^1.5*exp(-1*shield_cr(g,d)/shield(g,i,d)^1.5); %[-] Cheng 2002 (also in Carrillo 2021), 2.53<s<2.69, 0.73<S<1.2%, 0.068<D<1.27ft, 0.093<Q<1.118ft3/s | theta_cr=0.05
            qb_C(g,d) = Qb_C(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]

            %A theta^0.5*(theta-theta_cr)
            Qb_N(g,d) = 12*shield(g,i,d)^0.5*(max(shield(g,i,d)-shield_cr(g,d),0)); %[-] Nielsen 1992 (also in Soulsby 2005 and Camenen 2005) (different in Carrillo 2021) (payed book download)
            qb_N(g,d) = Qb_N(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            %A(f) theta^0.5*(theta-theta_cr)
            Qb_S(g,d) = 4.2*S^0.6*(u(g,i,d)/ust(g,i))*shield(g,i,d)^0.5*(max(shield(g,i,d)-shield_cr(g,d),0)); %[-] Smart 1984
            qb_S(g,d) = Qb_S(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            %A (theta^0.5-theta_cr^0.5)*(theta-theta_cr)
            Qb_AM(g,d) = 17*(max(shield(g,i,d)-shield_cr(g,d),0))*(max(shield(g,i,d)^0.5-shield_cr(g,d)^0.5,0)); %[-] Ashida and Michiue 1972 as in D3D and Carrillo 2021 (I have the paper, but it is in Japanese) | theta_cr=0.05
            qb_AM(g,d) = Qb_AM(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]

            %other
            p_EF(g,d) = min((1+((pi/6)*1/(max(shield(g,i,d)-shield_cr(g,d),0)))^4)^(-0.25),1); %[-] probability | beta=1 Engelund and Fredsoe 1982 as in Garcia and Parker 1991
            Qb_EF(g,d) = 5*p_EF(g,d)*(max(shield(g,i,d)^0.5-0.7*shield_cr(g,d)^0.5,0)); %[-] Engelund and Fredsoe 1976 (9.3*(pi/6)=~5) | theta_cr=0.05,0.06
            qb_EF(g,d) = Qb_EF(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            Qb_P(g,d) = 11.2*(max(shield(g,i,d)-shield_cr(g,d),0))^4.5/shield(g,i,d)^3; %[-] Parker 1979 as in Carrillo 2021 or Parker 1982 as in Kleinhans 2005 (payed download) | theta_cr=0.03                Qb_P(g,d) = 11.2*(shield(g,d)-shield_cr_nc(g,d))^4.5/shield(g,d)^3; %[-] Parker 1979 as in Carrillo 2021 or Parker 1982 as in Kleinhans 2005 (payed download) | theta_cr=0.03
            qb_P(g,d) = Qb_P(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            Qb_vR(g,d) = 0.053*Dst(g,d)^-0.3*T0(g,d)^2.1; %[-] van Rijn 1984, 200<D<2000mum
            qb_vR(g,d) = Qb_vR(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            Qb_vR2(g,d) = 0.1*Dst(g,d)^(-0.3)*S0(g,d)^1.5; %[-] van Rijn 1984 as in Kleinhans 2005
            qb_vR2(g,d) = Qb_vR2(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]

            %no theta_cr
            Qb_EH(g,d) = 0.1*shield(g,i,d)^2.5 / (f(g,i,d)/4); %[-] Engelund and Hansen as in D3D
            qt_EH(g,d) = Qb_EH(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            %qb_EH(g,d) = 0.05*u(g,i,d)^5/(grav(g)^0.5*(8*grav(g)/f(g,i,d))^1.5*R^2*D(d)); %[m3/ms]Engelund and Hansen as in D3D
            Qb_W(g,d) = 12*shield(g,i,d)^1.5; %[-] Wilson 1966 as in Soulsby 2005
            qb_W(g,d) = Qb_W(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]
            Qb_E42(g,d) = 2.1*exp(-0.391/shield(g,i,d)); %[-] Einstein 1942 as in Carrillo 2021, 1.25<s<4.25, 0.315<D<28.6mm, Qb<0.4
            qb_E42(g,d) = Qb_E42(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]

            G2(g,d) = rhos*D(d)^2*ust(g,i)^2/(rho*14*v^2);  %[?] Grain flow number (Bagnold 1966), 150<G2<6000:tan(phi)=-0.236, G2<150:tan(phi)=0.7, G2>6000:tan(phi)=0.374
            if D(d)>0.000015 && D(d)<0.00006
                eb(g,d) = -0.012*log(3.28*u(g,i,d))+0.15;
            elseif D(d)>=0.00006 && D(d)<0.0002
                eb(g,d) = -0.013*log(3.28*u(g,i,d))+0.145;
            elseif D(d)>=0.0002 && D(d)<0.0007
                eb(g,d) = -0.016*log(3.28*u(g,i,d))+0.139;
            elseif D(d)>=0.0007
                eb(g,d) = -0.028*log(3.28*u(g,i,d))+0.135;
            else
                eb(g,d) = NaN;
            end
            if G2(g,d)>150 && G2(g,d)<6000
                qb_B(g,d) = eb(g,d)*u(g,i,d)*tau(g,i)/((rhos-rho)*grav(g)*cos(S*(-0.236-tan(S)))); %[m3/ms]
            elseif G2(g,d)<=150
                qb_B(g,d) = eb(g,d)*u(g,i,d)*tau(g,i)/((rhos-rho)*grav(g)*cos(S*(0.7-tan(S))));
            elseif G2(g,d)>=6000
                qb_B(g,d) = eb(g,d)*u(g,i,d)*tau(g,i)/((rhos-rho)*grav(g)*cos(S*(0.374-tan(S))));
            else
                qb_B(g,d) = NaN;
            end
            Qb_B(g,d) = qb_B(g,d)/(grav(g)*R*D(d)^3)^0.5; %[-] Bagnold 1966 as in Kleinhans 2005

        end
    end
end
% More equations in: Cheng 2002
% Check: Yalin 1963, Carrillo 2021

%% Bedload calculations
for a = 1
    qb(:,:,1) = qb_E42;
    qb(:,:,2) = qb_MPM ;
    qb(:,:,3) = qb_MP;
    qb(:,:,4) = qb_E50;
    qb(:,:,5) = qb_B;
    qb(:,:,6) = qb_W;
    qb(:,:,7) = qb_AM;
    qb(:,:,8) = qb_FB;
    qb(:,:,9) = qb_EF;
    qb(:,:,10) = qb_P;
    qb(:,:,11) = qb_S;
    qb(:,:,12) = qb_vR;
    qb(:,:,13) = qb_vR2;
    qb(:,:,14) = qb_N;
    qb(:,:,15) = qb_R;
    qb(:,:,16) = qb_HJ;
    qb(:,:,17) = qb_C;
    qb(:,:,18) = qb_Ca;
    qb(:,:,19) = qb_WP1;
    qb(:,:,20) = qb_WP2;
    %qb_E42 qb_MPM qb_MP qb_E50 qb_B qb_W qb_AM qb_FB qb_EF qb_P qb_S qb_vR qb_vR2 qb_N qb_R qb_HJ qb_C qb_Ca qb_WP1 qb_WP2
    clear qb_E42 qb_MPM qb_MP qb_E50 qb_B qb_W qb_AM qb_FB qb_EF qb_P qb_S qb_vR qb_vR2 qb_N qb_R qb_HJ qb_C qb_Ca qb_WP1 qb_WP2
    clear Qb_MPM Qb_MP  Qb_FB Qb_R Qb_HJ Qb_WP1 Qb_WP2 Qb_Ca Qb_C Qb_N Qb_S Qb_AM Qb_EF Qb_P Qb_vR Qb_vR2 Qb_W Qb_E42 Qb_B Qb_EH
    qb_ratio = squeeze(qb(1,:,:)./qb(2,:,:));
end

%% Suspended transport
for i = hselect
    for g=1:Lg
        for d=1:LD
            %reference height
            a_E50(d) = 2*D(d); %Einstein 1950
            a_EF(d) = 2*D(d); %Engelund and Fredsoe 1976
            a_IK = 0.05*h(i); %Itakura and Kishi 1980 as in Garcia and Parker 1991
            a_CR = 0.05*h(i); %Celik and Rodi 1984 as in Garcia and Parker 1991
            a_AF = 0.05*h(i); %Akiyama and Fukushima 1986 as in Garcia and Parker 1991
            a_GP = 0.05*h(i); %Garcia and Parker 1991
            a_WP = 0.05*h(i); %Wright and Parker 2004
            a_dL = 0.1*h(i); %de Leeuw 2020
            a_vR(d) = max(0.01*h(i),ks(d)); %van Rijn 1984
            a_SM(g,d) = 26.3*(max(shield(g,i,d)-shield_cr(g,d),0))*D(d)+ks(d); %Smith and McLean 1977
            a_ML(g,d) = 0.68*(tau(g,i)/tau_cr(g,d))*D(d)/(1+(0.0204*log(D(d)*100)^2+0.022*log(D(d)*100)+0.0709)*(tau(g,i)/tau_cr(g,d))); %McLean 1992

            %Rouse number
            Rouse(g,d) = ws(g,d)/(Kappa*ust(g,i)); %Rouse 1937
            %Rouse_dL(g,d) = (ust(g,i)/ws(g,d))^-0.45; %de Leeuw 2020
            %Rouse_dL(g,d) = 0.145 * (ust(g,i,d)/ws(g,d))^-0.46 * (ust(g,i,d)^2/u(g,i,d)^2)^-0.3; %de Leeuw 2020

            %Ca / Es
            Ca_E50(g,d) = (1/23.2)*Qb_E50(g,d)/shield(g,i,d)^0.5; %Einstein 1950
            labda_EF(g,d) = (max((shield(g,i,d)-shield_cr(g,d)-(1*p_EF(g,d)*pi/6))/(0.027*(R+1)*shield(g,i,d)),0))^0.5;
            Ca_EF(g,d) = 0.65/(1+labda_EF(g,d)^-1)^3; %Engelund and Fredsoe 1976
            Ca_SM(g,d) = 0.65*0.0024*S0(g,d)/(1+0.0024*S0(g,d)); %Smith and McLean 1977
            Ca_ML(g,d) = 0.065*0.004*S0(g,d)/(1+0.004*S0(g,d)); %McLean 1992
            Ca_vR(g,d) = 0.015*(D(d)/a_vR(d))*S0(g,d)^1.5/Dst(g,d)^0.3; %van Rijn 1984
            Z_AF(g,d) = Rep(g,d)^0.5*ust(g,i)/ws(g,d);
            if Z_AF(g,d)<5
                Ca_AF(g,d) = 0;
            elseif Z_AF(g,d)>13.2
                Ca_AF(g,d) = 0.3;
            else
                Ca_AF(g,d) = 0.000000000003*Z_AF(g,d)^10*(1-(5/Z_AF(g,d))); %Akiyama and Fukushima 1986 as in Garcia and Parker 1991
            end
            Z_GP(g,d) = Rep(g,d)^0.6*ust(g,i)/ws(g,d);
            Ca_GP(g,d) = 0.00000013*Z_GP(g,d)^5/(1+(0.00000013/0.3)*Z_GP(g,d)^5); %Garcia and Parker 1991
            Z_WP(g,d) = Rep(g,d)^0.6*ust(g,i)/ws(g,d);
            Ca_WP(g,d) = 0.00000078*Z_WP(g,d)^5/(1+(0.00000078/0.3)*Z_WP(g,d)^5); %Wright and Parker 2004
            %Ca_dL(g,d) = 0.000474*(ust(g,i)/ws(g,d))^1.77*Fr(g,i,d)^1.18; %de Leeuw 2020, equation 25
            Ca_dL(g,d) = 4.74*10^-4*((ust(g,i)/ws(g,d))^1.5*Fr(g,i,d)-0.015)^1.18 / (1+3*(4.74*10^-4*((ust(g,i)/ws(g,d))^1.5*Fr(g,i,d)-0.015)^1.18)); %de Leeuw 2020, equation 26b
            fun = @(eta) (((1-eta)./eta).*(0.05/(1-0.05))).^Rouse(g,d);
            I_CR(g,d) = integral(fun,0.05,1);
            Cm_CR(g,d) = 0.034*(1-(ks(d)/h(i))^0.06)*ust(g,i)^2*u(g,i,d)/(grav(g)*R*h(i)*ws(g,d));
            Ca_CR(g,d) = 1.13*Cm_CR(g,d)/I_CR(g,d); %Celik and Rodi 1984
            A_IK(g,d) = 0.143/shield(g,i,d)-2;
            fun2 = @(eta2) exp(-1.*eta2.^2);
            Omega_IK(g,d) = (shield(g,i,d)/0.143)*(2+((exp(-1*A_IK(g,d)^2))/(integral(fun2,A_IK(g,d),Inf))))-1;
            Ca_IK(g,d) = max(0.008*(0.14*ust(g,i)*Omega_IK(g,d)/(ws(g,d)*shield(g,i,d))-1),0); %Itakura and Kishi 1980 as in Garcia and Parker 1991, goes negative without max(CA_IK,0)...

            %Integral
            fun = @(z) (ust(g,i) * Ca_dL(g,d) / Kappa) * ( ((h(i)-z)./z) * (a_dL/(h(i)-a_dL)) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_dL(g,d) = integral(fun,a_dL,h(i)); %de Leeuw 2020
            fun = @(z) (ust(g,i) * Ca_E50(g,d) / Kappa) * ( ((h(i)-z)./z) * (a_E50(d)/(h(i)-a_E50(d))) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_E50(g,d) = integral(fun,a_E50(d),h(i)); %Einstein 1950
            fun = @(z) (ust(g,i) * Ca_EF(g,d) / Kappa) * ( ((h(i)-z)./z) * (a_EF(d)/(h(i)-a_EF(d))) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_EF(g,d) = integral(fun,a_EF(d),h(i)); %Engelund and Fredsoe 1976
            fun = @(z) (ust(g,i) * Ca_IK(g,d) / Kappa) * ( ((h(i)-z)./z) * (a_IK/(h(i)-a_IK)) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_IK(g,d) = integral(fun,a_IK,h(i)); %Itakura and Kishi 1980
            fun = @(z) (ust(g,i) * Ca_CR(g,d) / Kappa) * ( ((h(i)-z)./z) * (a_CR/(h(i)-a_CR)) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_CR(g,d) = integral(fun,a_CR,h(i)); %Celik and Rodi 1984
            fun = @(z) (ust(g,i) * Ca_AF(g,d) / Kappa) * ( ((h(i)-z)./z) * (a_AF/(h(i)-a_AF)) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_AF(g,d) = integral(fun,a_AF,h(i)); %Akiyama and Fukushima 1986
            fun = @(z) (ust(g,i) * Ca_GP(g,d) / Kappa) * ( ((h(i)-z)./z) * (a_GP/(h(i)-a_GP)) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_GP(g,d) = integral(fun,a_GP,h(i)); %Garcia and Parker 1991
            fun = @(z) (ust(g,i) * Ca_WP(g,d) / Kappa) * ( ((h(i)-z)./z) * (a_WP/(h(i)-a_WP)) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_WP(g,d) = integral(fun,a_WP,h(i)); %Wright and Parker 2004
            fun = @(z) (ust(g,i) * Ca_vR(g,d) / Kappa) * ( ((h(i)-z)./z) * (a_vR(d)/(h(i)-a_vR(d))) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_vR(g,d) = integral(fun,a_vR(d),h(i)); %van Rijn 1984
            fun = @(z) (ust(g,i) * Ca_SM(g,d) / Kappa) * ( ((h(i)-z)./z) * (a_SM(g,d)/(h(i)-a_SM(g,d))) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_SM(g,d) = integral(fun,a_SM(g,d),h(i)); %Smith and McLean 1977
            fun = @(z) (ust(g,i) * Ca_ML(g,d) / Kappa) * ( ((h(i)-z)./z) * (a_ML(g,d)/(h(i)-a_ML(g,d))) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_ML(g,d) = integral(fun,a_ML(g,d),h(i)); %McLean 1992

            % de Leeuw et al 2020 for total transport calculations
            ub(g,d) = 0.6 * u(g,i,d); %de Leeuw et al 2020 eq 20
            hb(g,d) = min(0.6 * (Fr(g,i,d) * (D(d)/h(i))^2)^0.3 * h(i), h(i)); %de Leeuw et al 2020 eq 21
            Cb(g,d) = qb(g,d,8) / (hb(g,d) * ub(g,d)); %de Leeuw et al 2020 eq 17
            fun = @(z) (ust(g,i) * Cb(g,d) / Kappa) * ( ((h(i)-z)./z) * (hb(g,d)/(h(i)-hb(g,d))) ).^Rouse(g,d) .* log(z/(0.033*ks(d)));
            qs_dL2(g,d) = integral(fun,hb(g,d),h(i)); 
        end
    end
end
clear labda_EF Z_AF Z_GP Z_WP fun I_CR A_IK fun2 Omega_IK Cm_CR p_EF Qb_E50 ub hb Cb qs_dL2_a qs_dL2_b

%% Suspended calculations
for a = 1
    qs(:,:,1) = qs_E50;
    qs(:,:,2) = qs_EF;
    qs(:,:,3) = qs_SM;
    qs(:,:,4) = qs_IK;
    qs(:,:,5) = qs_CR;
    qs(:,:,6) = qs_vR;
    qs(:,:,7) = qs_AF;
    qs(:,:,8) = qs_GP;
    qs(:,:,9) = qs_ML;
    qs(:,:,10) = qs_WP;
    qs(:,:,11) = qs_dL;
    qs(:,:,12) = qs_dL2;
    %qs_E50 qs_EF qs_SM qs_IK qs_CR qs_vR qs_AF qs_GP qs_ML qs_WP qs_dL qs_dL2
    clear qs_E50 qs_EF qs_SM qs_IK qs_CR qs_vR qs_AF qs_GP qs_ML qs_WP qs_dL qs_dL2
    clear Ca_dL Ca_E50 Ca_EF Ca_IK Ca_CR Ca_AF Ca_GP Ca_WP Ca_vR Ca_SM Ca_ML Ca_dL2
    clear a_dL a_E50 a_EF a_IK a_CR a_AF a_GP a_WP a_vR a_SM a_ML
    qs_ratio = squeeze(qs(1,:,:)./qs(2,:,:));
end

%% Total calculations
for a = 1
    qt_dL2 = squeeze(qb(:,:,8)) + squeeze(qs(:,:,12));

    qt_ratio_dL2 = squeeze( qt_dL2(1,:)./qt_dL2(2,:) );
    qt_ratio_EH = squeeze( qt_EH(1,:)./qt_EH(2,:) );

    suspended_percent_dL2 = squeeze(qs(:,:,12)) ./qt_dL2*100;
end

%% Plot settings
for a = 1
    pl.width = 18; %[cm] plot width
    pl.height = 23; %[cm] plot height
    pl.line = 1; %line thickness
    pl.line_ax = 0.75; %line thickness axes
    pl.fsz = 7; %font size axes
    pl.fsz2 = 8.5; %font size labels
    CMAPOBJ = clrmap('read','Gravity.clrmap');
    clr1 = clrmap(CMAPOBJ,Lg); %colourmap gravity, red-blue
    pl.vertical = 4; %number of vertical subplots
    pl.horizantal = 2; %number of horizontal subplots
    pl.lmarge = 0.065; %left margin
    pl.rmarge = 0.025; %right margin
    pl.bmarge = 0.055; %bottom margin
    pl.tmarge = 0.01; %top margin
    pl.intv = 0.05; %vertical space between graphs
    pl.inth = 0.065; %horizontal space between graphs
    pl.wt = (1 - pl.lmarge - pl.rmarge - ((pl.horizantal-1)*pl.inth)) / pl.horizantal;
    pl.ht = (1 - pl.bmarge - pl.tmarge - ((pl.vertical-1)*pl.intv)) / pl.vertical;
    y = 0.95; %placing subplot letter, y
    D50L = D(1);
    D50R = D(end);
    xd = 10^(log10(D50L)+( log10(D50R)-log10(D50L) )*0.01); %placing subplot letter, x
    alfa_unsuitable = 0.15;
    alfa_predictors = 0.5;
end

%% Plot
for a = 1
    close all; clc;
    f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[100 -10 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    close all; clc;
    f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[100 -10 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s1(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    %qb_E42 qb_MPM qb_MP qb_E50 qb_B qb_W qb_AM qb_FB qb_EF qb_P qb_S qb_vR qb_vR2 qb_N qb_R qb_HJ qb_C qb_Ca qb_WP1 qb_WP2
    for g = 1:Lg
        for j = 2 % acceptable, here for legend
            plot(D,qb(g,:,j),'-','color',[clr1(g,:) alfa_predictors],'linewidth',pl.line);
        end
    end
    for g = 1:Lg
        for j = 8 %used for total
            plot(D,qb(g,:,j),'-','color',clr1(g,:),'linewidth',pl.line*2);
        end
    end
    for g = 1:Lg
        for j = [1 3 5 6] % unsuitable, no critical value
            %plot(D,qb(g,:,j),':','color',[clr1(g,:) alfa_unsuitable],'linewidth',pl.line);
        end
        for j = [9 12 13] % unsuitable, reduction value as small grain sizes
            %plot(D,qb(g,:,j),':','color',[clr1(g,:) alfa_unsuitable],'linewidth',pl.line);
        end
        for j = [2 4 7 8 10 11 14 15 16 17 18 19 20] % acceptable
            plot(D,qb(g,:,j),'-','color',[clr1(g,:) alfa_predictors],'linewidth',pl.line);
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50L D50R]);
    ylim([10^-6.5 10^1]);
    l(1) = ylabel('Bedload transport (q_b) [m^3/ms]');
    xlabel('Grain size (D_{50}) [m]');
    t(1) = text(xd,10^(-6.5+(1--6.5)*y),' a');


    s1(2) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    %qs_E50 qs_EF qs_SM qs_IK qs_CR qs_vR qs_AF qs_GP qs_ML qs_WP qs_dL qs_dL2
    for j = 1 % acceptable, here for legend
        for g = 1:Lg
            plot(D,qs(g,:,j),'-','color',[clr1(g,:) alfa_predictors],'linewidth',pl.line);
        end
    end
    for j = 12 % used for total
        for g = 1:Lg
            plot(D,qs(g,:,j),'-','color',clr1(g,:),'linewidth',pl.line*2);
        end
    end
    for g = 1:Lg
        for j = [4 5 7 8 10] % unsuitable, small range or too much transport big grain sizes
            %plot(D,qs(g,:,j),':','color',[clr1(g,:) alfa_unsuitable],'linewidth',pl.line);
        end
        for j = [1 2 3 6 9 11 12] % acceptable
            plot(D,qs(g,:,j),'-','color',[clr1(g,:) alfa_predictors],'linewidth',pl.line);
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50L D50R]);
    ylim([10^-6.5 10^1]);
    l(2) = ylabel('Suspended transport (q_s) [m^3/ms]');
    xlabel('Grain size (D_{50}) [m]');
    t(2) = text(xd,10^(-6.5+(1--6.5)*y),' b');


    s1(3) = axes('Position',[pl.lmarge 1-pl.tmarge-2*pl.ht-pl.intv pl.wt pl.ht]);
    hold on
    %qb_E42 qb_MPM qb_MP qb_E50 qb_B qb_W qb_AM qb_FB qb_EF qb_P qb_S qb_vR qb_vR2 qb_N qb_R qb_HJ qb_C qb_Ca qb_WP1 qb_WP2
    for j = 2 % acceptable, here for legend
        plot(D,qb_ratio(:,j),'-','color',[0 0 0 alfa_predictors],'linewidth',pl.line);
    end
    for j = 8 % used for total
        plot(D,qb_ratio(:,j),'-','color',[0 0 0],'linewidth',pl.line*2);
    end
    plot([D(1) D(end)], [1 1],'--g')
    for j = [1 3] % unsuitable, no critical value
        %plot(D,qb_ratio(:,j),':','color',[0 0 0 alfa_unsuitable],'linewidth',pl.line);
    end
    for j = [9 12 13] % unsuitbale, reducing value small grian sizes
        %plot(D,qb_ratio(:,j),':','color',[0 0 0 alfa_unsuitable],'linewidth',pl.line);
    end
    for j = [2 4 7 8 10 11 14 15 16 17 18 19 20] % acceptable
        plot(D,qb_ratio(:,j),'-','color',[0 0 0 alfa_predictors],'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50L D50R]);
    ylim([10^-0.4 10^1.1]);
    l(3) = ylabel('Mars/Earth ratio bedload transport [-]');
    xlabel('Grain size (D_{50}) [m]');
    t(3) = text(xd,10^(-0.4+(1.1--0.4)*y),' c');

    s1(4) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-2*pl.ht-pl.intv pl.wt pl.ht]);
    hold on
    %qs_E50 qs_EF qs_SM qs_IK qs_CR qs_vR qs_AF qs_GP qs_ML qs_WP qs_dL
    for j = 1 % acceptable, here for legend
        plot(D,qs_ratio(:,j),'-','color',[0 0 0 alfa_predictors],'linewidth',pl.line);
    end
    for j = 12 % used for total
        plot(D,qs_ratio(:,j),'-','color',[0 0 0],'linewidth',pl.line*2);
    end
    plot([D(1) D(end)], [1 1],'--g')
    for j = [4 5 7 8 10] % unsuitable, small range or too much transport big grain sizes
        %plot(D,qs_ratio(:,j),':','color',[0 0 0 alfa_unsuitable],'linewidth',pl.line);
    end
    for j = [1 2 3 6 9 11 12] % acceptable
        plot(D,qs_ratio(:,j),'-','color',[0 0 0 alfa_predictors],'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50L D50R]);
    ylim([10^-0.4 10^1.1]);
    l(3) = ylabel('Mars/Earth ratio suspended transport [-]');
    xlabel('Grain size (D_{50}) [m]');
    t(4) = text(xd,10^(-0.4+(1.1--0.4)*y),' d');

    s1(5) = axes('Position',[pl.lmarge 1-pl.tmarge-3*pl.ht-2*pl.intv pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(D,qt_dL2(g,:),'-','color',clr1(g,:),'linewidth',pl.line*2)
    end
    for g = 1:Lg
        plot(D,qt_EH(g,:),':','color',[clr1(g,:) 0.5],'linewidth',pl.line*2)
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    %xlim([D50(1) D50(end)]);
    xlim([D50L D50R]);
    ylim([10^-5 10^1]);
    l(3) = ylabel('Total transport (q_t) [m^3/ms]');
    xlabel('Grain size (D_{50}) [m]');
    t(5) = text(xd,10^(-4+(1--4)*y),' e');

    s1(6) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-3*pl.ht-2*pl.intv pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(D,suspended_percent_dL2(g,:),'-','color',clr1(g,:),'linewidth',pl.line*2)
    end
    set(gca,'Xscale','log');
    %xlim([D50(1) D50(end)]);
    xlim([D50L D50R]);
    ylim([0 100]);
    l(3) = ylabel('Suspended transport of total [%]');
    xlabel('Grain size (D_{50}) [m]');
    t(6) = text(xd,100*y,' f');

    s1(7) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    plot(D,qt_ratio_dL2,'-','color',[0 0 0],'linewidth',pl.line*2)
    plot(D,qt_ratio_EH,':','color',[0 0 0 0.5],'linewidth',pl.line*2)
    plot([D(1) D(end)], [1 1],'--g')
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50L D50R]);
    ylim([10^-0.4 10^1.1]);
    l(3) = ylabel('Mars/Earth ratio total transport [-]');
    xlabel('Grain size (D_{50}) [m]');
    t(7) = text(xd,10^(-0.4+(1.1--0.4)*y),' g');
    ytick = [0.1 0.5 1 2 10];
    yticks(ytick);

    %leg1 = legend(s1(1),['Mars h = ' num2str(h(hselect)) ' m'],['Earth h = ' num2str(h(hselect)) ' m'],'Mars - used for total','Earth - used for total','unsuitable','location','northeast');
    %leg2 = legend(s1(2),['Mars h = ' num2str(h(hselect)) ' m'],['Earth h = ' num2str(h(hselect)) ' m'],'Mars - used for total','Earth - used for total','unsuitable','location','northeast');
    leg1 = legend(s1(1),['Mars h = ' num2str(h(hselect)) ' m'],['Earth h = ' num2str(h(hselect)) ' m'],'Mars - used for total','Earth - used for total','location','northeast');
    leg2 = legend(s1(2),['Mars h = ' num2str(h(hselect)) ' m'],['Earth h = ' num2str(h(hselect)) ' m'],'Mars - used for total','Earth - used for total','location','northeast');
    %leg3 = legend(s1(3),'Mars/Earth ratio','Used for total','Mars = Earth','unsuitable','location','nort');
    %leg4 = legend(s1(4),'Mars/Earth ratio','Used for total','Mars = Earth','unsuitable','location','north');
    leg3 = legend(s1(3),'Mars/Earth ratio','Used for total','Mars = Earth','location','north');
    leg4 = legend(s1(4),'Mars/Earth ratio','Used for total','Mars = Earth','location','north');
    leg5 = legend(s1(5),'Mars (de Leeuw et al., 2020)','Earth (de Leeuw et al., 2020)','Mars (Engelund and Hansen, 1976)','Earth (Engelund and Hansen, 1976)','location','northeast');
    leg6 = legend(s1(6),'Mars (de Leeuw et al., 2020)','Earth (de Leeuw et al., 2020)','location','northeast');
    leg7 = legend(s1(7),'Mars/Earth ratio (de Leeuw et al., 2020)','Mars/Earth ration (Engelund and Hansen 1976)','Mars = Earth','location','northeast');

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
