%% Flow on Mars 7
% Author: Lisanne Braat (lisannebraat@gmail.com)
% Date: 2022-10-06

%% Initialize
clear variables
close all
clc

output = 'FlowOnMars_exportfig';
addpath(genpath('Checkout'));

%% Input parameters
for a = 1
    W = 200;
    S = 0.001;
    TC = 4;
    rho = 1000;

    grav = [3.7 9.8];
    Q = 250:250:15000;
    h = 0.5:0.5:15;
    u = 0.5:0.1:6.5;
    tau = 1:1:100;

    rhos = 2900;
    load('FlowOnMars5_GrainSizeDistribution.mat');
    Drough = D90;
    clear D10 D50 D90;
    D50 = D; clear D;
    CDryB = 1750;
    C1 = 20;
    C2 = 1;
    phi = 30;
    Kappa = 0.41;

    % Derived input parameters
    TK = TC+273.15; %[degrees K] Temperature
    mu = 1.856e-14*exp(4209/TK+0.04527*TK-3.376e-5*TK^2); %[Pa | Ns/m2] Dynamic viscosity
    v = mu/rho; %[m2/s] Viscosity
    R = (rhos-rho)/rho; %[-] Relative density
    n = 1-CDryB/rhos; %[-] Porosity
    ks = 3*Drough; %[m] Nikuradse
end
clear TC TK mu CDryB Drough

%% allocate
for a=1
    Lg = length(grav);
    LD = length(D50);
    LQ = length(Q);
    Lh = length(h);
    Lu = length(u);
    Ltau = length(tau);

    h_Q = NaN(Lg,LQ);
    C_Q = NaN(Lg,LQ);
    u_Q = NaN(Lg,LQ);
    tau_Q = NaN(Lg,LQ);
    ust_Q = NaN(Lg,LQ);

    Q_h = NaN(Lg,Lh);
    C_h = NaN(Lg,Lh);
    u_h = NaN(Lg,Lh);
    tau_h = NaN(Lg,Lh);
    ust_h = NaN(Lg,Lh);

    h_u = NaN(Lg,Lu);
    C_u = NaN(Lg,Lu);
    Q_u = NaN(Lg,Lu);
    tau_u = NaN(Lg,Lu);
    ust_u = NaN(Lg,Lu);

    h_tau = NaN(Lg,Ltau);
    C_tau = NaN(Lg,Ltau);
    u_tau = NaN(Lg,Ltau);
    Q_tau = NaN(Lg,Ltau);
    ust_tau = NaN(Lg,Ltau);

    ws = NaN(Lg,LD);
    Dst = NaN(Lg,LD);
    shield_Q = NaN(Lg,LD,LQ);
    shield_h = NaN(Lg,LD,Lh);
    shield_u = NaN(Lg,LD,Lu);
    shield_tau = NaN(Lg,LD,Ltau);

    shield_cr_Q = NaN(Lg,LD,LQ);
    shield_cr_h = NaN(Lg,LD,Lh);
    shield_cr_u = NaN(Lg,LD,Lu);
    shield_cr_tau = NaN(Lg,LD,Ltau);

    Qb_E50_Q = NaN(Lg,LD,LQ);
    Qb_E50_h = NaN(Lg,LD,Lh);
    Qb_E50_u = NaN(Lg,LD,Lu);
    Qb_E50_tau = NaN(Lg,LD,Ltau);
    qb_E50_Q = NaN(Lg,LD,LQ);
    qb_E50_h = NaN(Lg,LD,Lh);
    qb_E50_u = NaN(Lg,LD,Lu);
    qb_E50_tau = NaN(Lg,LD,Ltau);

    a_E50 = NaN(1,LD);
    Rouse_Q = NaN(Lg,LD,LQ);
    Rouse_h = NaN(Lg,LD,Lh);
    Rouse_u = NaN(Lg,LD,Lu);
    Rouse_tau = NaN(Lg,LD,Ltau);
    Ca_E50_Q = NaN(Lg,LD,LQ);
    Ca_E50_h = NaN(Lg,LD,Lh);
    Ca_E50_u = NaN(Lg,LD,Lu);
    Ca_E50_tau = NaN(Lg,LD,Ltau);

    qs_E50_Q = NaN(Lg,LD,LQ);
    qs_E50_h = NaN(Lg,LD,Lh);
    qs_E50_u = NaN(Lg,LD,Lu);
    qs_E50_tau = NaN(Lg,LD,Ltau);

    qb_E50_frac_Q = NaN(Lg,LD,LQ);
    qb_E50_frac_h = NaN(Lg,LD,Lh);
    qb_E50_frac_u = NaN(Lg,LD,Lu);
    qb_E50_frac_tau = NaN(Lg,LD,Ltau);

    qs_E50_frac_Q = NaN(Lg,LD,LQ);
    qs_E50_frac_h = NaN(Lg,LD,Lh);
    qs_E50_frac_u = NaN(Lg,LD,Lu);
    qs_E50_frac_tau = NaN(Lg,LD,Ltau);

    qt_E50_frac_Q = NaN(Lg,LD,LQ);
    qt_E50_frac_h = NaN(Lg,LD,Lh);
    qt_E50_frac_u = NaN(Lg,LD,Lu);
    qt_E50_frac_tau = NaN(Lg,LD,Ltau);

    qt_E50_sum_Q = NaN(Lg,LQ);
    qt_E50_sum_h = NaN(Lg,Lh);
    qt_E50_sum_u = NaN(Lg,Lu);
    qt_E50_sum_tau = NaN(Lg,Ltau);
end

%% Hydro parameters
for a = 1
    h_guess = 3;
    for i = 1:LQ
        for g = 1:Lg
            h_Q(g,i) = h_guess+1;
            while abs(h_guess-h_Q(g,i))>0.0001
                h_guess = h_Q(g,i);
                C_Q(g,i) = 5.75*grav(g)^0.5*log10(12*h_guess/(ks));
                u_Q(g,i) = C_Q(g,i)*(h_guess*S)^0.5;
                h_Q(g,i) = Q(i)/(W*u_Q(g,i));
            end
            tau_Q(g,i) = rho*grav(g)*h_Q(g,i)*S;
            ust_Q(g,i) = (tau_Q(g,i)/rho)^0.5;
        end
    end
    clear i g h_guess

    for i = 1:Lh
        for g = 1:Lg
            C_h(g,i) = 5.74*grav(g)^0.5*log10(12*h(i)/(ks));
            u_h(g,i) = C_h(g,i)*(h(i)*S)^0.5;
            Q_h(g,i) = h(i)*W*u_h(g,i);
            tau_h(g,i) = rho*grav(g)*h(i)*S;
            ust_h(g,i) = (tau_h(g,i)/rho)^0.5;
        end
    end
    clear i g

    h_guess = 3;
    for i = 1:Lu
        for g = 1:Lg
            h_u(g,i) = h_guess+1;
            while abs(h_guess-h_u(g,i))>0.0001
                h_guess = h_u(g,i);
                C_u(g,i) = 5.75*grav(g)^0.5*log10(12*h_guess/(ks));
                h_u(g,i) = (u(i)/C_u(g,i))^2/S;
            end
            Q_u(g,i) = u(i)*W*h_u(g,i);
            tau_u(g,i) = rho*grav(g)*h_u(g,i)*S;
            ust_u(g,i) = (tau_u(g,i)/rho)^0.5;
        end
    end
    clear i g h_guess

    for i = 1:Ltau
        for g = 1:Lg
            h_tau(g,i) = tau(i)/(rho*grav(g)*S);
            C_tau(g,i) = 5.75*grav(g)^0.5*log10(12*h_tau(g,i)/(ks));
            u_tau(g,i) = C_tau(g,i)*(h_tau(g,i)*S)^0.5;
            Q_tau(g,i) = u_tau(g,i)*W*h_tau(g,i);
            ust_tau(g,i) = (tau(i)/rho)^0.5;
        end
    end
    clear i g
end
clear C_Q C_h C_u C_tau u_Q u_h u_tau ks Q_h Q_u Q_tau

%% Sediment parameters
for d = 1:LD
    for g = 1:Lg
        ws(g,d) = R*grav(g)*D50(d)^2/(C1*v+(0.75*C2*R*grav(g)*D50(d)^3)^0.5); %[m/s] Settling velocity (Ferguson and Church 2004)
        %Dst(g,d) = D50(d)*(R*grav(g)/v^2)^(1/3); %[-] Dimensionless particle parameter/Bonnefille number

        for i = 1:LQ
            shield_Q(g,d,i) = tau_Q(g,i)/((rhos-rho)*grav(g)*D50(d)); %[-] Shields parameter/Particle mobility parameter
        end
        for i = 1:Lh
            shield_h(g,d,i) = tau_h(g,i)/((rhos-rho)*grav(g)*D50(d));
        end
        for i = 1:Lu
            shield_u(g,d,i) = tau_u(g,i)/((rhos-rho)*grav(g)*D50(d));
        end
        for i = 1:Ltau
            shield_tau(g,d,i) = tau(i)/((rhos-rho)*grav(g)*D50(d));
        end
    end
end
clear C1 C2 d g i tau_Q tau_h tau_u

%% Critical threshold of motion
for a = 1
    acc = 1e-5;
    for d = 1:LD
        for g = 1:Lg
            for i = 1:LQ
                %Zanke 2003 iterated
                shieldini = shield_Q(g,d,i);
                minimise = 1;
                while max(minimise)>acc
                    tau_ini = shieldini *((rhos-rho)*grav(g)*D50(d));
                    ust_ini = (tau_ini/rho)^0.5;
                    Rest_ini = D50(d)*ust_ini/v;
                    uprmsb_ust = 0.31*Rest_ini*exp(-0.1*Rest_ini)+1.8*exp(-0.88*D50(d)/(shieldini*R*D50(d)/sin(S)))*(1-exp(-0.1*Rest_ini));
                    Pt = 1 - exp(-0.08*Rest_ini);
                    B = (1-Pt)*(2.5.*log(Rest_ini)+5.25)+8.5*Pt;
                    uy_ust = ((1-Pt)/Rest_ini^2+Pt/(2.5*log(1)+B^2))^-0.5;
                    ub_ust = 0.8+0.9*uy_ust;
                    uprmsb_ub = uprmsb_ust / ub_ust;
                    K = 1 + 3e-8/((rhos-rho)*D50(d)^2);
                    shield1 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                        ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
                    minimise = abs(shield1 - shieldini);
                    shieldini = shield1;
                end
                shield_cr_Q(g,d,i) = shieldini;
            end
        end
    end
    clear Rest_ini ust_ini tau_ini shieldini Pt B uy_ust ub_ust uprmsb_ub uprmsb_ust K shield1 minimise d g

    for d = 1:LD
        for g = 1:Lg
            for i = 1:Lh
                %Zanke 2003 iterated
                shieldini = shield_h(g,d,i);
                minimise = 1;
                while max(minimise)>acc
                    tau_ini = shieldini *((rhos-rho)*grav(g)*D50(d));
                    ust_ini = (tau_ini/rho)^0.5;
                    Rest_ini = D50(d)*ust_ini/v;
                    uprmsb_ust = 0.31*Rest_ini*exp(-0.1*Rest_ini)+1.8*exp(-0.88*D50(d)/(shieldini*R*D50(d)/sin(S)))*(1-exp(-0.1*Rest_ini));
                    Pt = 1 - exp(-0.08*Rest_ini);
                    B = (1-Pt)*(2.5.*log(Rest_ini)+5.25)+8.5*Pt;
                    uy_ust = ((1-Pt)/Rest_ini^2+Pt/(2.5*log(1)+B^2))^-0.5;
                    ub_ust = 0.8+0.9*uy_ust;
                    uprmsb_ub = uprmsb_ust / ub_ust;
                    K = 1 + 3e-8/((rhos-rho)*D50(d)^2);
                    shield1 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                        ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
                    minimise = abs(shield1 - shieldini);
                    shieldini = shield1;
                end
                shield_cr_h(g,d,i) = shieldini;
            end
        end
    end
    clear Rest_ini ust_ini tau_ini shieldini Pt B uy_ust ub_ust uprmsb_ub uprmsb_ust K shield1 minimise d g

    for d = 1:LD
        for g = 1:Lg
            for i = 1:Lu
                %Zanke 2003 iterated
                shieldini = shield_u(g,d,i);
                minimise = 1;
                while max(minimise)>acc
                    tau_ini = shieldini *((rhos-rho)*grav(g)*D50(d));
                    ust_ini = (tau_ini/rho)^0.5;
                    Rest_ini = D50(d)*ust_ini/v;
                    uprmsb_ust = 0.31*Rest_ini*exp(-0.1*Rest_ini)+1.8*exp(-0.88*D50(d)/(shieldini*R*D50(d)/sin(S)))*(1-exp(-0.1*Rest_ini));
                    Pt = 1 - exp(-0.08*Rest_ini);
                    B = (1-Pt)*(2.5.*log(Rest_ini)+5.25)+8.5*Pt;
                    uy_ust = ((1-Pt)/Rest_ini^2+Pt/(2.5*log(1)+B^2))^-0.5;
                    ub_ust = 0.8+0.9*uy_ust;
                    uprmsb_ub = uprmsb_ust / ub_ust;
                    K = 1 + 3e-8/((rhos-rho)*D50(d)^2);
                    shield1 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                        ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
                    minimise = abs(shield1 - shieldini);
                    shieldini = shield1;
                end
                shield_cr_u(g,d,i) = shieldini;
            end
        end
    end
    clear Rest_ini ust_ini tau_ini shieldini Pt B uy_ust ub_ust uprmsb_ub uprmsb_ust K shield1 minimise d g

    for d = 1:LD
        for g = 1:Lg
            for i = 1:Ltau
                %Zanke 2003 iterated
                shieldini = shield_tau(g,d,i);
                minimise = 1;
                while max(minimise)>acc
                    tau_ini = shieldini *((rhos-rho)*grav(g)*D50(d));
                    ust_ini = (tau_ini/rho)^0.5;
                    Rest_ini = D50(d)*ust_ini/v;
                    uprmsb_ust = 0.31*Rest_ini*exp(-0.1*Rest_ini)+1.8*exp(-0.88*D50(d)/(shieldini*R*D50(d)/sin(S)))*(1-exp(-0.1*Rest_ini));
                    Pt = 1 - exp(-0.08*Rest_ini);
                    B = (1-Pt)*(2.5.*log(Rest_ini)+5.25)+8.5*Pt;
                    uy_ust = ((1-Pt)/Rest_ini^2+Pt/(2.5*log(1)+B^2))^-0.5;
                    ub_ust = 0.8+0.9*uy_ust;
                    uprmsb_ub = uprmsb_ust / ub_ust;
                    K = 1 + 3e-8/((rhos-rho)*D50(d)^2);
                    shield1 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                        ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
                    minimise = abs(shield1 - shieldini);
                    shieldini = shield1;
                end
                shield_cr_tau(g,d,i) = shieldini;
            end
        end
    end
    clear acc Rest_ini ust_ini tau_ini shieldini Pt B uy_ust ub_ust uprmsb_ub uprmsb_ust K shield1 minimise phi n d g rho rhos
end
clear a i d g v S

%% Bedload transport
for g = 1:Lg
    for d = 1:LD
        for i = 1:LQ
            Qb_E50_Q(g,d,i) = 3.97*(max(shield_Q(g,d,i)-shield_cr_Q(g,d,i),0))^1.5;
            qb_E50_Q(g,d,i) = Qb_E50_Q(g,d,i)*(grav(g)*R*D50(d)^3)^0.5;
        end
        for i = 1:Lh
            Qb_E50_h(g,d,i) = 3.97*(max(shield_h(g,d,i)-shield_cr_h(g,d,i),0))^1.5;
            qb_E50_h(g,d,i) = Qb_E50_h(g,d,i)*(grav(g)*R*D50(d)^3)^0.5;
        end
        for i = 1:Lu
            Qb_E50_u(g,d,i) = 3.97*(max(shield_u(g,d,i)-shield_cr_u(g,d,i),0))^1.5;
            qb_E50_u(g,d,i) = Qb_E50_u(g,d,i)*(grav(g)*R*D50(d)^3)^0.5;
        end
        for i = 1:Ltau
            Qb_E50_tau(g,d,i) = 3.97*(max(shield_tau(g,d,i)-shield_cr_tau(g,d,i),0))^1.5;
            qb_E50_tau(g,d,i) = Qb_E50_tau(g,d,i)*(grav(g)*R*D50(d)^3)^0.5;
        end
    end
end
clear R g d i a

%% Suspended transport
for d = 1:LD
    a_E50(d) = 2*D50(d);
    for g = 1:Lg
        for i = 1:LQ
            Rouse_Q(g,d,i) = ws(g,d)/(Kappa*ust_Q(g,i));
            Ca_E50_Q(g,d,i) = (1/32.2)*Qb_E50_Q(g,d,i)/shield_Q(g,d,i)^0.5;
            fun = @(z) Ca_E50_Q(g,d,i) * ( ((h_Q(g,i)-z)./z) * (a_E50(d)/(h_Q(g,i)-a_E50(d))) ).^Rouse_Q(g,d,i);
            qs_E50_Q(g,d,i) = integral(fun,a_E50(d),h_Q(g,i));
        end
        for i = 1:Lh
            Rouse_h(g,d,i) = ws(g,d)/(Kappa*ust_h(g,i));
            Ca_E50_h(g,d,i) = (1/32.2)*Qb_E50_h(g,d,i)/shield_h(g,d,i)^0.5;
            fun = @(z) Ca_E50_h(g,d,i) * ( ((h(i)-z)./z) * (a_E50(d)/(h(i)-a_E50(d))) ).^Rouse_h(g,d,i);
            qs_E50_h(g,d,i) = integral(fun,a_E50(d),h(i));
        end
        for i = 1:Lu
            Rouse_u(g,d,i) = ws(g,d)/(Kappa*ust_u(g,i));
            Ca_E50_u(g,d,i) = (1/32.2)*Qb_E50_u(g,d,i)/shield_u(g,d,i)^0.5;
            fun = @(z) Ca_E50_u(g,d,i) * ( ((h_u(g,i)-z)./z) * (a_E50(d)/(h_u(g,i)-a_E50(d))) ).^Rouse_u(g,d,i);
            qs_E50_u(g,d,i) = integral(fun,a_E50(d),h_u(g,i));
        end
        for i = 1:Ltau
            Rouse_tau(g,d,i) = ws(g,d)/(Kappa*ust_tau(g,i));
            Ca_E50_tau(g,d,i) = (1/32.2)*Qb_E50_tau(g,d,i)/shield_tau(g,d,i)^0.5;
            fun = @(z) Ca_E50_tau(g,d,i) * ( ((h_tau(g,i)-z)./z) * (a_E50(d)/(h_tau(g,i)-a_E50(d))) ).^Rouse_tau(g,d,i);
            qs_E50_tau(g,d,i) = integral(fun,a_E50(d),h_tau(g,i));
        end
    end
end
clear g d i Kappa a_E50 fun ws
clear shield_Q shield_h shield_u shield_tau shield_cr_Q shield_cr_h shield_cr_u shield_cr_tau
clear ust_Q ust_h ust_u ust_tau Rouse_Q Rouse_h Rouse_u Rouse_tau
clear Ca_E50_Q Ca_E50_h Ca_E50_u Ca_E50_tau

%% multiply sediment classes with bed availability from lognormal sediment curve
for a = 1
    for g = 1:Lg
        for i = 1:LQ
            qb_E50_frac_Q(g,:,i) = qb_E50_Q(g,:,i).*bins_y_perc/100;
            qs_E50_frac_Q(g,:,i) = qs_E50_Q(g,:,i).*bins_y_perc/100;
        end
        for i = 1:Lh
            qb_E50_frac_h(g,:,i) = qb_E50_h(g,:,i).*bins_y_perc/100;
            qs_E50_frac_h(g,:,i) = qs_E50_h(g,:,i).*bins_y_perc/100;
        end
        for i = 1:Lu
            qb_E50_frac_u(g,:,i) = qb_E50_u(g,:,i).*bins_y_perc/100;
            qs_E50_frac_u(g,:,i) = qs_E50_u(g,:,i).*bins_y_perc/100;
        end
        for i = 1:Ltau
            qb_E50_frac_tau(g,:,i) = qb_E50_tau(g,:,i).*bins_y_perc/100;
            qs_E50_frac_tau(g,:,i) = qs_E50_tau(g,:,i).*bins_y_perc/100;
        end
    end
    qt_E50_frac_Q = qb_E50_frac_Q + qs_E50_frac_Q;
    qt_E50_frac_h = qb_E50_frac_h + qs_E50_frac_h;
    qt_E50_frac_u = qb_E50_frac_u + qs_E50_frac_u;
    qt_E50_frac_tau = qb_E50_frac_tau + qs_E50_frac_tau;
    qt_E50_sum_Q = sum(qt_E50_frac_Q,2);
    qt_E50_sum_h = sum(qt_E50_frac_h,2);
    qt_E50_sum_u = sum(qt_E50_frac_u,2);
    qt_E50_sum_tau = sum(qt_E50_frac_tau,2);
end
%clear qb_E50 qb_FB qb_EF qb_vR qs_E50 qs_dL qs_EF qs_vR qb_EH

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
    %CMAPOBJ = clrmap('read','Discharge.clrmap');
    %clr2 = clrmap(CMAPOBJ,Lh);
    %clr3 = clrmap(CMAPOBJ,LQ);
    %CMAPOBJ = clrmap('read','GrainSize.clrmap');
    %clr4 = clrmap(CMAPOBJ,LD);
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
    xg = grav(1) + (grav(end)-grav(1))*0.04;
    xq = Q(1) + (Q(end)-Q(1))*0.04;
end

%% Plot
for a = 1
    close all
    clc
    f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[1 1 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s1(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on;
    for g = 1:Lg;
        plot(Q,squeeze(qt_E50_sum_Q(g,:,:)*W),'-','color',clr1(g,:),'linewidth',pl.line)
    end
    set(gca,'Yscale','log');
    xlim([Q(1) Q(end)]);
    ylim([10^1 10^4]);
    l(1) = xlabel('Discharge [m^3/s]');
    l(2) = ylabel('Total transport [m^3/s]');
    t(1) = text(Q(1)+(Q(end)-Q(1))*0.04,10^(1+(4-1)*y),' a');

    s1(2) = axes('Position',[pl.lmarge+pl.wt+pl.inth 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on;
    for g = 1:Lg;
        plot(h,squeeze(qt_E50_sum_h(g,:,:)*W),'-','color',clr1(g,:),'linewidth',pl.line)
    end
    set(gca,'Yscale','log');
    xlim([h(1) h(end)]);
    ylim([10^1 10^4]);
    l(3) = xlabel('Water depth [m]');
    l(4) = ylabel('Total transport [m^3/s]');
    t(2) = text(h(1)+(h(end)-h(1))*0.04,10^(1+(4-1)*y),' b');

    s1(3) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on;
    for g = 1:Lg;
        plot(u,squeeze(qt_E50_sum_u(g,:,:)*W),'-','color',clr1(g,:),'linewidth',pl.line)
    end
    set(gca,'Yscale','log');
    xlim([u(1) u(end)]);
    ylim([10^1 10^4]);
    l(5) = xlabel('Velocity [m/s]');
    l(6) = ylabel('Total transport [m^3/s]');
    t(3) = text(u(1)+(u(end)-u(1))*0.04,10^(1+(4-1)*y),' c');

    s1(4) = axes('Position',[pl.lmarge+pl.wt+pl.inth pl.bmarge pl.wt pl.ht]);
    hold on;
    for g = 1:Lg;
        plot(tau,squeeze(qt_E50_sum_tau(g,:,:)*W),'-','color',clr1(g,:),'linewidth',pl.line)
    end
    set(gca,'Yscale','log');
    xlim([tau(1) tau(end)]);
    ylim([10^1 10^4]);
    l(7) = xlabel('Shear stress [m^2/s]');
    l(8) = ylabel('Total transport [m^3/s]');
    t(4) = text(tau(1)+(tau(end)-tau(1))*0.04,10^(1+(4-1)*y),' d');

    leg1 = legend(s1(1),'Mars','Earth','location','east');

    t(5) = text(s1(1),min(Q)+(max(Q)-min(Q)),1.5*10^1,'W = 200 m + S = 0.001 m/m + Q = 250-15000 m^3/s   ','HorizontalAlignment', 'right');
    t(6) = text(s1(2),min(h)+(max(h)-min(h)),1.5*10^1,'W = 200 m + S = 0.001 m/m + h = 0.5-15 m   ','HorizontalAlignment', 'right');
    t(7) = text(s1(3),min(u)+(max(u)-min(u)),1.5*10^1,'W = 200 m + S = 0.001 m/m + u = 0.5-6.5 m/s   ','HorizontalAlignment', 'right');
    t(8) = text(s1(4),min(tau)+(max(tau)-min(tau)),1.5*10^1,'W = 200 m + S = 0.001 m/m + \tau = 1-100 N/m^2   ','HorizontalAlignment', 'right');

    set(s1,'box','on','Layer','top', ...
        'XMinorTick','on','YMinorTick','on', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);

    set(gcf,'renderer','painters');
    print(f1,'-dpng',[output '/FlowOnMars7'],'-r400');
    print(f1,'-dpdf',[output '/FlowOnMars7'],'-r400');
end


