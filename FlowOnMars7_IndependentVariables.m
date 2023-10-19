%% Flow on Mars 7 - Other independent variables
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
    W = 200; %[m] Channel width
    S = 0.001; %[m/m] Channel slope
    TC = 4; %[degrees C] Temperature
    rho = 1000; %[kg/m3] Water density

    grav = [3.7 9.8]; %[m/s2] Gravitational acceleration
    Q = 250:250:15000; %[m3/s] Discharge
    h = 0.5:0.5:15; %[m] Water depth
    u = 0.5:0.1:6.5; %[m/s] Velocity
    tau = 1:1:100; %[N/m2] Bed shear stress

    % Sediment parameters
    rhos = 2900; %[kg/m3] Sediment density
    load('FlowOnMars5_GrainSizeDistribution.mat');
    clear D10 D90;
    CDryB = 1750; %[kg/m3] Dry bed density
    C1 = 20; %18-24, smooth to rough sphere, affects small diameters, Ferguson and Church 2004
    C2 = 1; %0.4-1.2, smooth to rough sphere, affects big diameters, Ferguson and Church 2004
    phi = 30; %angle of repose in degrees
    Kappa = 0.41; %[-] Kappa

    % Derived input parameters
    TK = TC+273.15; %[degrees K] Temperature
    mu = 1.856e-14*exp(4209/TK+0.04527*TK-3.376e-5*TK^2); %[Pa | Ns/m2] Dynamic viscosity
    v = mu/rho; %[m2/s] Viscosity
    ks = 2.5*D50; %[m] Nikurandse roughness height

    R = (rhos-rho)/rho; %[-] Relative density
    n = 1-CDryB/rhos; %[-] Porosity

    Lg = length(grav);
    LD = length(D);
    LQ = length(Q);
    Lh = length(h);
    Lu = length(u);
    Ltau = length(tau);
end
clear a TC TK mu CDryB D10 D90

%% allocate
for a = 1
    h_Q = NaN(Lg,LQ);
    h_tau = NaN(Lg,Ltau);
    h_u = NaN(Lg,Lu);
    f_h = NaN(1,Lh);
    f_Q = NaN(Lg,LQ);
    f_tau = NaN(Lg,Ltau);
    f_u = NaN(Lg,Lu);
    u_h = NaN(Lg,Lh);
    u_Q = NaN(Lg,LQ);
    u_tau = NaN(Lg,Ltau);
    tau_h = NaN(Lg,Lh);
    tau_Q = NaN(Lg,LQ);
    tau_u = NaN(Lg,Lu);
    ust_h = NaN(Lg,Lh);
    ust_Q = NaN(Lg,LQ);
    ust_tau = NaN(1,Ltau);
    ust_u = NaN(Lg,Lu);
    Q_h = NaN(Lg,Lh);
    Q_tau = NaN(Lg,Ltau);
    Q_u = NaN(Lg,Lu);
    Fr_Q = NaN(Lg,LQ);
    Fr_h = NaN(Lg,Lh);
    Fr_u = NaN(Lg,Lu);
    Fr_tau = NaN(Lg,Ltau);

    ws = NaN(Lg,LD);
    shield_Q = NaN(Lg,LQ,LD);
    shield_h = NaN(Lg,Lh,LD);
    shield_u = NaN(Lg,Lu,LD);
    shield_tau = NaN(Lg,Ltau,LD);

    shield_cr_Q = NaN(Lg,LD);
    shield_cr_h = NaN(Lg,LD);
    shield_cr_u = NaN(Lg,LD);
    shield_cr_tau = NaN(Lg,LD);

    Qb_FB_Q = NaN(Lg,LQ,LD);
    Qb_FB_h = NaN(Lg,Lh,LD);
    Qb_FB_u = NaN(Lg,Lu,LD);
    Qb_FB_tau = NaN(Lg,Ltau,LD);
    qb_FB_Q = NaN(Lg,LQ,LD);
    qb_FB_h = NaN(Lg,Lh,LD);
    qb_FB_u = NaN(Lg,Lu,LD);
    qb_FB_tau = NaN(Lg,Ltau,LD);

    a_Q = NaN(Lg,LQ);
    a_h = NaN(1,Lh);
    a_u = NaN(Lg,Lu);
    a_tau = NaN(Lg,Ltau);
    Rouse_Q = NaN(Lg,LQ,LD);
    Rouse_h = NaN(Lg,Lh,LD);
    Rouse_u = NaN(Lg,Lu,LD);
    Rouse_tau = NaN(Lg,Ltau,LD);
    ub_Q = NaN(Lg,LQ);
    ub_h = NaN(Lg,Lh);
    ub_u = NaN(1,Lu);
    ub_tau = NaN(Lg,Ltau);
    hb_Q = NaN(Lg,LQ,LD);
    hb_h = NaN(Lg,Lh,LD);
    hb_u = NaN(Lg,Lu,LD);
    hb_tau = NaN(Lg,Ltau,LD);
    Cb_Q = NaN(Lg,LQ,LD);
    Cb_h = NaN(Lg,Lh,LD);
    Cb_u = NaN(Lg,Lu,LD);
    Cb_tau = NaN(Lg,Ltau,LD);
    Ca_Q = NaN(Lg,LQ,LD);
    Ca_h = NaN(Lg,Lh,LD);
    Ca_u = NaN(Lg,Lu,LD);
    Ca_tau = NaN(Lg,Ltau,LD);
    qs_Q = NaN(Lg,LQ,LD);
    qs_h = NaN(Lg,Lh,LD);
    qs_u = NaN(Lg,Lu,LD);
    qs_tau = NaN(Lg,Ltau,LD);

    qs_Q_a = NaN(Lg,LQ,LD);
    qs_Q_b = NaN(Lg,LQ,LD);
    qs_h_a = NaN(Lg,Lh,LD);
    qs_h_b = NaN(Lg,Lh,LD);
    qs_u_a = NaN(Lg,Lu,LD);
    qs_u_b = NaN(Lg,Lu,LD);
    qs_tau_a = NaN(Lg,Ltau,LD);
    qs_tau_b = NaN(Lg,Ltau,LD);
    
    qb_FB_frac_Q = NaN(Lg,LQ,LD);
    qb_FB_frac_h = NaN(Lg,Lh,LD);
    qb_FB_frac_u = NaN(Lg,Lu,LD);
    qb_FB_frac_tau = NaN(Lg,Ltau,LD);
    qs_dL_frac_Q = NaN(Lg,LQ,LD);
    qs_dL_frac_h = NaN(Lg,Lh,LD);
    qs_dL_frac_u = NaN(Lg,Lu,LD);
    qs_dL_frac_tau = NaN(Lg,Ltau,LD);
end
clear a 

%% Hydro parameters
for a = 1
    h_guess = 3;
    for g = 1:Lg
        for i = 1:LQ
            h_Q(g,i) = h_guess+1;
            while abs(h_guess-h_Q(g,i))>0.0001
                h_guess = h_Q(g,i);
                f_Q(g,i) = 8/(5.75*log10(12*h_guess/ks))^2; %[-] Friction factor
                u_Q(g,i) = (8*grav(g)*h_guess*S/f_Q(g,i))^0.5; %[m/s] Velocity
                h_Q(g,i) = Q(i)/(W*u_Q(g,i)); %[m] Water depth
            end
            tau_Q(g,i) = rho*grav(g)*h_Q(g,i)*S; %[N/m2] Bed shear stress
            ust_Q(g,i) = (tau_Q(g,i)/rho)^0.5; %[m/s] Shear velocity
            Fr_Q(g,i) = u_Q(g,i)/(grav(g)*h_Q(g,i))^0.5; %[-] Froude number
        end
    end

    for g = 1:Lg
        for i = 1:Lh
            f_h(i) = 8/(5.74*log10(12*h(i)/ks))^2; %[-] Friction factor
            u_h(g,i) = (8*grav(g)*h(i)*S/f_h(i))^0.5; %[m/s] Velocity
            Q_h(g,i) = h(i)*W*u_h(g,i); %[m3/s] Discharge
            tau_h(g,i) = rho*grav(g)*h(i)*S; %[N/m2] Bed shear stress
            ust_h(g,i) = (tau_h(g,i)/rho)^0.5; %[m/s] Shear velocity
            Fr_h(g,i) = u_h(g,i)/(grav(g)*h(i))^0.5; %[-] Froude number
        end
    end

    h_guess = 3;
    for g = 1:Lg
        for i = 1:Lu
            h_u(g,i) = h_guess+1;
            while abs(h_guess-h_u(g,i))>0.0001
                h_guess = h_u(g,i);
                f_u(g,i) = 8/(5.75*log10(12*h_guess/ks))^2; %[-] Friction factor
                h_u(g,i) = u(i)^2*f_u(g,i)/(8*grav(g)*S); %[m] Water depth
            end
            Q_u(g,i) = u(i)*W*h_u(g,i); %[m3/s] Discharge
            tau_u(g,i) = rho*grav(g)*h_u(g,i)*S; %[N/m2] Bed shear stress
            ust_u(g,i) = (tau_u(g,i)/rho)^0.5; %[m/s] Shear velocity
            Fr_u(g,i) = u(i)/(grav(g)*h_u(g,i))^0.5; %[-] Froude number
        end
    end

    for g = 1:Lg
        for i = 1:Ltau
            h_tau(g,i) = tau(i)/(rho*grav(g)*S); %[m] Water depth
            f_tau(g,i) = 8/(5.75*log10(12*h_tau(g,i)/ks))^2; %[-] Friction factor
            u_tau(g,i) = (8*grav(g)*h_tau(g,i)*S/f_tau(g,i))^0.5; %[m/s] Velocity
            Q_tau(g,i) = u_tau(g,i)*W*h_tau(g,i); %[m3/s] Discharge
            ust_tau(i) = (tau(i)/rho)^0.5; %[m/s] Shear velocity
            Fr_tau(g,i) = u_tau(g,i)/(grav(g)*h_tau(g,i))^0.5; %[-] Froude number
        end
    end
end
clear a i g d h_guess f_Q f_h f_u f_tau Q_h Q_u Q_tau

%% Sediment parameters
for d = 1:LD
    for g = 1:Lg
        ws(g,d) = R*grav(g)*D(d)^2/(C1*v+(0.75*C2*R*grav(g)*D(d)^3)^0.5); %[m/s] Settling velocity (Ferguson and Church 2004)

        for i = 1:LQ
            shield_Q(g,i,d) = tau_Q(g,i)/((rhos-rho)*grav(g)*D(d)); %[-] Shields parameter/Particle mobility parameter
        end
        for i = 1:Lh
            shield_h(g,i,d) = tau_h(g,i)/((rhos-rho)*grav(g)*D(d));
        end
        for i = 1:Lu
            shield_u(g,i,d) = tau_u(g,i)/((rhos-rho)*grav(g)*D(d));
        end
        for i = 1:Ltau
            shield_tau(g,i,d) = tau(i)/((rhos-rho)*grav(g)*D(d));
        end
    end
end
clear d g i C1 C2 tau_Q tau_h tau_u

%% Critical threshold of motion
for a = 1
    i = 1; 
    d = 1;
    acc = 1e-5;
    for g = 1:Lg
        %Zanke 2003 iterated
        shieldini = shield_Q(g,i,d);
        minimise = 1;
        while max(minimise)>acc
            tau_ini = shieldini *((rhos-rho)*grav(g)*D50);
            ust_ini = (tau_ini/rho)^0.5;
            Rest_ini = D50*ust_ini/v;
            uprmsb_ust = 0.31*Rest_ini*exp(-0.1*Rest_ini)+1.8*exp(-0.88*D50/(shieldini*R*D50/sin(S)))*(1-exp(-0.1*Rest_ini));
            Pt = 1 - exp(-0.08*Rest_ini);
            B = (1-Pt)*(2.5.*log(Rest_ini)+5.25)+8.5*Pt;
            uy_ust = ((1-Pt)/Rest_ini^2+Pt/(2.5*log(1)+B^2))^-0.5;
            ub_ust = 0.8+0.9*uy_ust;
            uprmsb_ub = uprmsb_ust / ub_ust;
            K = 1 + 3e-8/((rhos-rho)*D50^2);
            shield1 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
            minimise = abs(shield1 - shieldini);
            shieldini = shield1;
        end
        shield_cr_Q(g,:) = shieldini*(D/D50).^-0.9;
    end
    clear Rest_ini ust_ini tau_ini shieldini Pt B uy_ust ub_ust uprmsb_ub uprmsb_ust K shield1 minimise

    for g = 1:Lg
        %Zanke 2003 iterated
        shieldini = shield_h(g,i,d);
        minimise = 1;
        while max(minimise)>acc
            tau_ini = shieldini *((rhos-rho)*grav(g)*D50);
            ust_ini = (tau_ini/rho)^0.5;
            Rest_ini = D50*ust_ini/v;
            uprmsb_ust = 0.31*Rest_ini*exp(-0.1*Rest_ini)+1.8*exp(-0.88*D50/(shieldini*R*D50/sin(S)))*(1-exp(-0.1*Rest_ini));
            Pt = 1 - exp(-0.08*Rest_ini);
            B = (1-Pt)*(2.5.*log(Rest_ini)+5.25)+8.5*Pt;
            uy_ust = ((1-Pt)/Rest_ini^2+Pt/(2.5*log(1)+B^2))^-0.5;
            ub_ust = 0.8+0.9*uy_ust;
            uprmsb_ub = uprmsb_ust / ub_ust;
            K = 1 + 3e-8/((rhos-rho)*D50^2);
            shield1 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
            minimise = abs(shield1 - shieldini);
            shieldini = shield1;
        end
        shield_cr_h(g,:) = shieldini*(D/D50).^-0.9;
    end
    clear Rest_ini ust_ini tau_ini shieldini Pt B uy_ust ub_ust uprmsb_ub uprmsb_ust K shield1 minimise

    for g = 1:Lg
        %Zanke 2003 iterated
        shieldini = shield_u(g,i,d);
        minimise = 1;
        while max(minimise)>acc
            tau_ini = shieldini *((rhos-rho)*grav(g)*D50);
            ust_ini = (tau_ini/rho)^0.5;
            Rest_ini = D50*ust_ini/v;
            uprmsb_ust = 0.31*Rest_ini*exp(-0.1*Rest_ini)+1.8*exp(-0.88*D50/(shieldini*R*D50/sin(S)))*(1-exp(-0.1*Rest_ini));
            Pt = 1 - exp(-0.08*Rest_ini);
            B = (1-Pt)*(2.5.*log(Rest_ini)+5.25)+8.5*Pt;
            uy_ust = ((1-Pt)/Rest_ini^2+Pt/(2.5*log(1)+B^2))^-0.5;
            ub_ust = 0.8+0.9*uy_ust;
            uprmsb_ub = uprmsb_ust / ub_ust;
            K = 1 + 3e-8/((rhos-rho)*D50^2);
            shield1 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
            minimise = abs(shield1 - shieldini);
            shieldini = shield1;
        end
        shield_cr_u(g,:) = shieldini*(D/D50).^-0.9;
    end
    clear Rest_ini ust_ini tau_ini shieldini Pt B uy_ust ub_ust uprmsb_ub uprmsb_ust K shield1 minimise

    for g = 1:Lg
        %Zanke 2003 iterated
        shieldini = shield_tau(g,i,d);
        minimise = 1;
        while max(minimise)>acc
            tau_ini = shieldini *((rhos-rho)*grav(g)*D50);
            ust_ini = (tau_ini/rho)^0.5;
            Rest_ini = D50*ust_ini/v;
            uprmsb_ust = 0.31*Rest_ini*exp(-0.1*Rest_ini)+1.8*exp(-0.88*D50/(shieldini*R*D50/sin(S)))*(1-exp(-0.1*Rest_ini));
            Pt = 1 - exp(-0.08*Rest_ini);
            B = (1-Pt)*(2.5.*log(Rest_ini)+5.25)+8.5*Pt;
            uy_ust = ((1-Pt)/Rest_ini^2+Pt/(2.5*log(1)+B^2))^-0.5;
            ub_ust = 0.8+0.9*uy_ust;
            uprmsb_ub = uprmsb_ust / ub_ust;
            K = 1 + 3e-8/((rhos-rho)*D50^2);
            shield1 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
            minimise = abs(shield1 - shieldini);
            shieldini = shield1;
        end
        shield_cr_tau(g,:) = shieldini*(D/D50).^-0.9;
    end
    clear acc Rest_ini ust_ini tau_ini shieldini Pt B uy_ust ub_ust uprmsb_ub uprmsb_ust K shield1 minimise 
end
clear a g d i v S phi n rho rhos D50

%% Bedload transport - de Leeuw et al. 2020
for g = 1:Lg
    for d = 1:LD
        for i = 1:LQ
            Qb_FB_Q(g,i,d) = 5.7*(max(shield_Q(g,i,d)-shield_cr_Q(g,d),0))^1.5;
            qb_FB_Q(g,i,d) = Qb_FB_Q(g,i,d)*(grav(g)*R*D(d)^3)^0.5;
        end
        for i = 1:Lh
            Qb_FB_h(g,i,d) = 5.7*(max(shield_h(g,i,d)-shield_cr_h(g,d),0))^1.5;
            qb_FB_h(g,i,d) = Qb_FB_h(g,i,d)*(grav(g)*R*D(d)^3)^0.5;
        end
        for i = 1:Lu
            Qb_FB_u(g,i,d) = 5.7*(max(shield_u(g,i,d)-shield_cr_u(g,d),0))^1.5;
            qb_FB_u(g,i,d) = Qb_FB_u(g,i,d)*(grav(g)*R*D(d)^3)^0.5;
        end
        for i = 1:Ltau
            Qb_FB_tau(g,i,d) = 5.7*(max(shield_tau(g,i,d)-shield_cr_tau(g,d),0))^1.5;
            qb_FB_tau(g,i,d) = Qb_FB_tau(g,i,d)*(grav(g)*R*D(d)^3)^0.5;
        end
    end
end
clear g d i R shield_cr_Q shield_cr_h shield_cr_u shield_cr_tau
clear Qb_FB_Q Qb_FB_h Qb_FB_u Qb_FB_tau

%% Suspended transport - de Leeuw et al. 2020
for d = 1:LD
    for g = 1:Lg
        for i = 1:LQ
            Rouse_Q(g,i,d) = ws(g,d)/(Kappa*ust_Q(g,i));
            a_Q(g,i) = 0.01*h_Q(g,i);
            ub_Q(g,i) = 0.6 * u_Q(g,i);
            hb_Q(g,i,d) = min(0.6 * (Fr_Q(g,i) * (D(d)/h_Q(g,i))^2)^0.3 * h_Q(g,i),h_Q(g,i));
            Cb_Q(g,i,d) = qb_FB_Q(g,i,d) / (hb_Q(g,i,d) * ub_Q(g,i));
            fun = @(z) (ust_Q(g,i) * Cb_Q(g,i,d) / Kappa) * ( ((h_Q(g,i)-z)./z) * (hb_Q(g,i,d)/(h_Q(g,i)-hb_Q(g,i,d))) ).^Rouse_Q(g,i,d) .* log(z/(0.033*ks)); 
            qs_Q(g,i,d) = integral(fun,hb_Q(g,i,d),h_Q(g,i));
        end
        for i = 1:Lh
            Rouse_h(g,i,d) = ws(g,d)/(Kappa*ust_h(g,i));
            a_h(i) = 0.01*h(i);
            ub_h(g,i) = 0.6 * u_h(g,i);
            hb_h(g,i,d) = min(0.6 * (Fr_h(g,i) * (D(d)/h(i))^2)^0.3 * h(i),h(i));
            Cb_h(g,i,d) = qb_FB_h(g,i,d) / (hb_h(g,i,d) * ub_h(g,i));
                fun = @(z) (ust_h(g,i) * Cb_h(g,i,d) / Kappa) * ( ((h(i)-z)./z) * (hb_h(g,i,d)/(h(i)-hb_h(g,i,d))) ).^Rouse_h(g,i,d) .* log(z/(0.033*ks)); 
                qs_h(g,i,d) = integral(fun,hb_h(g,i,d),h(i));
        end
        for i = 1:Lu
            Rouse_u(g,i,d) = ws(g,d)/(Kappa*ust_u(g,i));
            a_u(g,i) = 0.01*h_u(g,i);
            ub_u(i) = 0.6 * u(i);
            hb_u(g,i,d) = min(0.6 * (Fr_u(g,i) * (D(d)/h_u(g,i))^2)^0.3 * h_u(g,i),h_u(g,i));
            Cb_u(g,i,d) = qb_FB_u(g,i,d) / (hb_u(g,i,d) * ub_u(i));
                fun = @(z) (ust_u(g,i) * Cb_u(g,i,d) / Kappa) * ( ((h_u(g,i)-z)./z) * (hb_u(g,i,d)/(h_u(g,i)-hb_u(g,i,d))) ).^Rouse_u(g,i,d) .* log(z/(0.033*ks)); 
                qs_u(g,i,d) = integral(fun,hb_u(g,i,d),h_u(g,i));
        end
        for i = 1:Ltau
            Rouse_tau(g,i,d) = ws(g,d)/(Kappa*ust_tau(i));
            a_tau(g,i) = 0.01*h_tau(g,i);
            ub_tau(g,i) = 0.6 * u_tau(g,i);
            hb_tau(g,i,d) = min(0.6 * (Fr_tau(g,i) * (D(d)/h_tau(g,i))^2)^0.3 * h_tau(g,i),h_tau(g,i));
            Cb_tau(g,i,d) = qb_FB_tau(g,i,d) / (hb_tau(g,i,d) * ub_tau(g,i));
                fun = @(z) (ust_tau(i) * Cb_tau(g,i,d) / Kappa) * ( ((h_tau(g,i)-z)./z) * (hb_tau(g,i,d)/(h_tau(g,i)-hb_tau(g,i,d))) ).^Rouse_tau(g,i,d) .* log(z/(0.033*ks)); 
                qs_tau(g,i,d) = integral(fun,hb_tau(g,i,d),h_tau(g,i));
        end
    end
end
clear g d i Kappa a_E50 fun ws ks
clear shield_Q shield_h shield_u shield_tau qs_tau_a qs_tau_b qs_u_a qs_u_b qs_h_a qs_h_b qs_Q_a qs_Q_b
clear ust_Q ust_h ust_u ust_tau Rouse_Q Rouse_h Rouse_u Rouse_tau u_Q u_tau u_h FR_h Fr_Q Fr_tau Fr_u h_Q h_tau h_u
clear a_tau a_u a_h a_Q ub_tau ub_u ub_h ub_Q hb_tau hb_u hb_h hb_Q Cb_tau Cb_u Cb_h Cb_Q Ca_tau Ca_u Ca_h Ca_Q

%% multiply sediment classes with bed availability from lognormal sediment curve
for a = 1
    for g = 1:Lg
        for i = 1:LQ
            qb_FB_frac_Q(g,i,:) = squeeze(qb_FB_Q(g,i,:)).*bins_y_perc'/100;
            qs_dL_frac_Q(g,i,:) = squeeze(qs_Q(g,i,:)).*bins_y_perc'/100;
        end
        for i = 1:Lh
            qb_FB_frac_h(g,i,:) = squeeze(qb_FB_h(g,i,:)).*bins_y_perc'/100;
            qs_dL_frac_h(g,i,:) = squeeze(qs_h(g,i,:)).*bins_y_perc'/100;
        end
        for i = 1:Lu
            qb_FB_frac_u(g,i,:) = squeeze(qb_FB_u(g,i,:)).*bins_y_perc'/100;
            qs_dL_frac_u(g,i,:) = squeeze(qs_u(g,i,:)).*bins_y_perc'/100;
        end
        for i = 1:Ltau
            qb_FB_frac_tau(g,i,:) = squeeze(qb_FB_tau(g,i,:)).*bins_y_perc'/100;
            qs_dL_frac_tau(g,i,:) = squeeze(qs_tau(g,i,:)).*bins_y_perc'/100;
        end
    end
    qt_dL_frac_Q = qb_FB_frac_Q + qs_dL_frac_Q;
    qt_dL_frac_h = qb_FB_frac_h + qs_dL_frac_h;
    qt_dL_frac_u = qb_FB_frac_u + qs_dL_frac_u;
    qt_dL_frac_tau = qb_FB_frac_tau + qs_dL_frac_tau;
    qt_dL_sum_Q = sum(qt_dL_frac_Q,3);
    qt_dL_sum_h = sum(qt_dL_frac_h,3);
    qt_dL_sum_u = sum(qt_dL_frac_u,3);
    qt_dL_sum_tau = sum(qt_dL_frac_tau,3);
end

%% Plot settings
for a = 1
    pl.width = 18; %[cm] plot width
    pl.height = 16; %[cm] plot height
    pl.line = 2; %line width
    pl.line_ax = 0.75; %line width axes
    pl.fsz = 7; %font size axes
    pl.fsz2 = 8.5; %font size labels
    CMAPOBJ = clrmap('read','Gravity.clrmap');
    clr1 = clrmap(CMAPOBJ,Lg); %colourmap gravity, red-blue
    pl.vertical = 2; %number of vertical subplots
    pl.horizantal = 2; %number of horizontal subplots
    pl.lmarge = 0.065; %left margin
    pl.rmarge = 0.025; %right margin
    pl.bmarge = 0.06; %bottom margin
    pl.tmarge = 0.01; %top margin
    pl.intv = 0.065; %vertical space between graphs
    pl.inth = 0.065; %horizontal space between graphs
    pl.wt = (1 - pl.lmarge - pl.rmarge - ((pl.horizantal-1)*pl.inth)) / pl.horizantal;
    pl.ht = (1 - pl.bmarge - pl.tmarge - ((pl.vertical-1)*pl.intv)) / pl.vertical;
    y = 0.94; %placing subplot letter, y
    aa = -3; %power y limit 1
    bb = 3.5; %power y limit 2
end

%% Plot
for a = 1
    close all
    clc
    f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[40 1 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s1(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on;
    for g = 1:Lg
        plot(Q,squeeze(qt_dL_sum_Q(g,:)*W),'-','color',clr1(g,:),'linewidth',pl.line)
    end
    plot(Q,(squeeze(qt_dL_sum_Q(1,:)*W)) ./ (squeeze(qt_dL_sum_Q(2,:)*W) ),'-.','color',[0 0 0],'linewidth',pl.line)
    plot([Q(1) Q(end)], [1 1],'--g')
    set(gca,'Yscale','log');
    xlim([Q(1) Q(end)]);
    l(1) = xlabel('Discharge [m^3/s]');
    l(2) = ylabel('Total transport [m^3/s]');
    t(1) = text(Q(1)+(Q(end)-Q(1))*0.04,10^(aa+(bb-aa)*y),' a');

    s1(2) = axes('Position',[pl.lmarge+pl.wt+pl.inth 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on;
    for g = 1:Lg
        plot(h,squeeze(qt_dL_sum_h(g,:)*W),'-','color',clr1(g,:),'linewidth',pl.line)
    end
    plot(h,(squeeze(qt_dL_sum_h(1,:)*W)) ./ (squeeze(qt_dL_sum_h(2,:)*W) ),'-.','color',[0 0 0],'linewidth',pl.line)
    plot([h(1) h(end)], [1 1],'--g')
    set(gca,'Yscale','log');
    xlim([h(1) h(end)]);
    l(3) = xlabel('Water depth [m]');
    l(4) = ylabel('Total transport [m^3/s]');
    t(2) = text(h(1)+(h(end)-h(1))*0.04,10^(aa+(bb-aa)*y),' b');

    s1(3) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on;
    for g = 1:Lg
        plot(u,squeeze(qt_dL_sum_u(g,:)*W),'-','color',clr1(g,:),'linewidth',pl.line)
    end
    plot(u,(squeeze(qt_dL_sum_u(1,:)*W)) ./ (squeeze(qt_dL_sum_u(2,:)*W) ),'-.','color',[0 0 0],'linewidth',pl.line)
    plot([u(1) u(end)], [1 1],'--g')
    set(gca,'Yscale','log');
    xlim([u(1) u(end)]);
    l(5) = xlabel('Velocity [m/s]');
    l(6) = ylabel('Total transport [m^3/s]');
    t(3) = text(u(1)+(u(end)-u(1))*0.04,10^(aa+(bb-aa)*y),' c');

    s1(4) = axes('Position',[pl.lmarge+pl.wt+pl.inth pl.bmarge pl.wt pl.ht]);
    hold on;
    for g = 1:Lg
        plot(tau,squeeze(qt_dL_sum_tau(g,:)*W),'-','color',clr1(g,:),'linewidth',pl.line)
    end
    plot(tau,(squeeze(qt_dL_sum_tau(1,:)*W)) ./ (squeeze(qt_dL_sum_tau(2,:)*W) ),'-.','color',[0 0 0],'linewidth',pl.line)
    plot([tau(1) tau(end)], [1 1],'--g')
    set(gca,'Yscale','log');
    xlim([tau(1) tau(end)]);
    l(7) = xlabel('Shear stress [m^2/s]');
    l(8) = ylabel('Total transport [m^3/s]');
    t(4) = text(tau(1)+(tau(end)-tau(1))*0.04,10^(aa+(bb-aa)*y),' d');

    leg1 = legend(s1(1),'Mars (de Leeuw et al. 2020)','Earth (de Leeuw et al. 2020)','Mars/Earth ratio','Mars = Earth');
    set(leg1,'Position',[0.22 0.615 0.25 0.1])

    t(5) = text(s1(1),min(Q)+(max(Q)-min(Q)),2*10^aa,'W = 200 m, S = 0.001 m/m, Q = 250-15000 m^3/s   ','HorizontalAlignment', 'right');
    t(6) = text(s1(2),min(h)+(max(h)-min(h)),2*10^aa,'W = 200 m, S = 0.001 m/m, h = 0.5-15 m   ','HorizontalAlignment', 'right');
    t(7) = text(s1(3),min(u)+(max(u)-min(u)),2*10^aa,'W = 200 m, S = 0.001 m/m, u = 0.5-6.5 m/s   ','HorizontalAlignment', 'right');
    t(8) = text(s1(4),min(tau)+(max(tau)-min(tau)),2*10^aa,'W = 200 m, S = 0.001 m/m, \tau = 1-100 N/m^2   ','HorizontalAlignment', 'right');

    set(s1,'box','on','Layer','top', ...
        'XMinorTick','on','YMinorTick','on', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);
    ylim(s1,[10^aa 10^bb]);

    set(gcf,'renderer','painters');
    print(f1,'-dpng',[output '/FlowOnMars7_IndependentVariables'],'-r400');
    print(f1,'-dpdf',[output '/FlowOnMars7_IndependentVariables'],'-r400');
end


