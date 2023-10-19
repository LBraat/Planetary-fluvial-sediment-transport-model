%% Flow on Mars 6 - Transport mixture
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
    Q = 2000; %[m3/s] Discharge
    W = 200; %[m] Channel width
    grav = 1:0.1:12; %[m/s2] Gravitational acceleration
    gselect = [28 89];
    S = 0.001; %[m/m] Channel slope
    TC = 4; %[degrees C] Temperature
    rho = 1000; %[kg/m3] Water density

    % Sediment parameters
    rhos = 2900; %[kg/m3] Sediment density
    load('FlowOnMars5_GrainSizeDistribution.mat');
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
end
clear a TC TK mu CDryB D10 D90

%% allocate
for a = 1
    h = NaN(1,Lg);
    Rw = NaN(1,Lg);
    f = NaN(1,Lg);
    u = NaN(1,Lg);
    Fr = NaN(1,Lg);
    %     Re = NaN(1,Lg);
    tau = NaN(1,Lg);
    ust = NaN(1,Lg);

    ws = NaN(Lg,LD);
    %     Dst = NaN(Lg,LD);
    shield = NaN(Lg,LD);

    shield_cr = NaN(Lg,LD);
    shield_cr_hide = NaN(Lg,LD);
    shield_cr_d50 = NaN(1,Lg);
    shield_cr_d50_hide = NaN(Lg,LD);

    Qb_FB1 = NaN(Lg,LD);
    Qb_FB2 = NaN(Lg,LD);
    Qb_FB3 = NaN(Lg,LD);
    Qb_FB4 = NaN(Lg,LD);

    Qb_FB = NaN(Lg,LD);
    qb_FB = NaN(Lg,LD);
    Qt_EH = NaN(Lg,LD);
    qt_EH = NaN(Lg,LD);

    a_dL = NaN(1,Lg);
    Rouse = NaN(Lg,LD);
    ub = NaN(1,Lg);
    hb = NaN(Lg,LD);
    Cb = NaN(Lg,LD);
    Ca_dL2 = NaN(Lg,LD);
    qs_dL2 = NaN(Lg,LD);
    qs_dL2_a = NaN(Lg,LD);
    qs_dL2_b = NaN(Lg,LD);

    qb_FB_mix = NaN(Lg,LD);
    qs_dL_mix = NaN(Lg,LD);
    qt_EH_mix = NaN(Lg,LD);
end
clear a

%% Hydro parameters Q_mars = Q_earth
h_guess = 3;
for g = 1:Lg
    h(g) = h_guess+1;
    while abs(h_guess-h(g))>0.0001
        h_guess = h(g);
        Rw(g) = (h_guess*W)/(h_guess+h_guess+W); %[m] Hydraulic radius
        f(g) = 8/(5.75*log10(12*h_guess/ks))^2; %[-] friction factor, use ks of mixture based on D50, so the values are not dependent on D
        u(g) = (8*grav(g)*Rw(g)*S/f(g))^0.5; %[m/s] Velocity
        h(g) = Q/(W*u(g)); %[m] Water depth
    end
    Fr(g) = u(g)/(grav(g)*h(g))^0.5; %[-] Froude number
    %Re(g) = u(g)*h(g)/v; %[-] Reynolds number
    tau(g) = rho*grav(g)*Rw(g)*S; %[N/m2] Bed shear stress
    ust(g) = (tau(g)/rho)^0.5; %[m/s] Shear velocity
end
clear g h_guess Rw Re Q

%% Sediment parameters
for d = 1:LD
    for g = 1:Lg
        ws(g,d) = R*grav(g)*D(d)^2/(C1*v+(0.75*C2*R*grav(g)*D(d)^3)^0.5); %[m/s] Settling velocity (Ferguson and Church 2004)
        %Dst(g,d) = D(d)*(R*grav(g)/v^2)^(1/3); %[-] Dimensionless particle parameter/Bonnefille number
        shield(g,d) = tau(g)/((rhos-rho)*grav(g)*D(d)); %[-] Shields parameter/Particle mobility parameter
    end
end
clear g d C1 C2 tau

%% Critical threshold of motion
acc = 1e-5;
for g = 1:Lg
    for d = 1:LD
        %Zanke 2003 iterated
        shieldini = shield(g,d);
        minimise = 1;
        while max(minimise)>acc
            tau_ini = shieldini *((rhos-rho)*grav(g)*D(d));
            ust_ini = (tau_ini/rho)^0.5;
            Rest_ini = D(d)*ust_ini/v;
            uprmsb_ust = 0.31*Rest_ini*exp(-0.1*Rest_ini)+1.8*exp(-0.88*D(d)/(shieldini*R*D(d)/sin(S)))*(1-exp(-0.1*Rest_ini));
            Pt = 1 - exp(-0.08*Rest_ini);
            B = (1-Pt)*(2.5.*log(Rest_ini)+5.25)+8.5*Pt;
            uy_ust = ((1-Pt)/Rest_ini^2+Pt/(2.5*log(1)+B^2))^-0.5;
            ub_ust = 0.8+0.9*uy_ust;
            uprmsb_ub = uprmsb_ust / ub_ust;
            K = 1 + 3e-8/((rhos-rho)*D(d)^2);
            shield1 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
            minimise = abs(shield1 - shieldini);
            shieldini = shield1;
        end
        shield_cr(g,d) = shieldini;
        shield_cr_hide(g,d) = shield_cr(g,d)*(D(d)/D50).^-0.9;
    end
end

% shield_cr D50, 1 value
for g = 1:Lg
    shield_cr_cal = InterX([D50 D50; 0 10],[D; shield_cr(g,:)]);
    shield_cr_d50(g) = shield_cr_cal(2);
    shield_cr_d50_hide(g,:) = shield_cr_cal(2)*(D/D50).^-0.9;
end
clear g d v S phi n rho rhos
clear acc Rest_ini ust_ini ust_Q_ini tau_ini tau_Q_ini shieldini shield_cr_cal Pt B uy_ust ub_ust uprmsb_ub uprmsb_ust K shield1 minimise 

%% Bedload transport - decide on what critical shear stress to use for a mixture
for a = 0
    if a == 1
        for g = 1:Lg
            for d = 1:LD
                %[-] Fernandez Luque and van Beek 1976 as in de Leeuw 2020 (payed download)
                Qb_FB1(g,d) = 5.7*(max(shield(g,d)-shield_cr(g,d),0))^1.5; %critical based on fractions
                Qb_FB2(g,d) = 5.7*(max(shield(g,d)-shield_cr_hide(g,d),0))^1.5; %critical based on fractions + hiding function
                Qb_FB3(g,d) = 5.7*(max(shield(g,d)-shield_cr_d50(g),0))^1.5; %critical based on mixture
                Qb_FB4(g,d) = 5.7*(max(shield(g,d)-shield_cr_d50_hide(g,d),0))^1.5; %critical based on mixture + hiding function
            end
        end

        close all
        figure()
        hold on
        plot(grav,sum(Qb_FB1,2),'r') %same as b
        plot(grav,sum(Qb_FB2,2),'g') %hiding effect is strongest (reducing transport)
        plot(grav,sum(Qb_FB3,2),'b') %same as r
        plot(grav,sum(Qb_FB4,2),'m') %small hiding effect
        set(gca,'Yscale','log');
        xlabel('g')
    end
end
% Decided to use 1 shield_cr for the whole mixture and add the hiding function.
%clear a g d shield_cr shield_cr_d50 shield_cr_hide

%% Bedload transport
for g = 1:Lg
    for d = 1:LD
        Qb_FB(g,d) = 5.7*(max(shield(g,d)-shield_cr_d50_hide(g,d),0))^(3/2); %[-] Fernandez Luque and van Beek 1976 as in de Leeuw 2020 (payed download)
        qb_FB(g,d) = Qb_FB(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms]

        Qt_EH(g,d) = 0.1*shield(g,d)^(5/2) / (f(g)/4); %[-] Engelund and Hansen - TOTAL LOAD %friction factor of EH is 4 times smaller than our f.
        qt_EH(g,d) = Qt_EH(g,d)*(grav(g)*R*D(d)^3)^0.5; %[m3/ms] - TOTAL LAOD
    end
end
clear g d f R shield shield_cr_d50_hide Qb_FB Qt_EH

%% Suspended transport
for g = 1:Lg
    for d = 1:LD
        %reference height
        a_dL(g) = 0.01*h(g);

        %Rouse number
        %Rouse(g,d) = ws(g,d)/(Kappa*ust(g));
        Rouse(g,d) = (ust(g)/ws(g,d))^-0.45; %0.145*(ust(g)/ws(g,d))^-0.46*(ust(g)^2/u(g)^2)^-0.3;

        %Ca / Es
        %Ca_dL(g,d) = 4.74*10^-4*(ust(g)/ws(g,d))^1.77*Fr(g)^1.18; %equation 25
        %Ca_dL(g,d) = 4.74*10^-4*(max((ust(g)/ws(g,d))^1.5*Fr(g)-0.015,0))^1.18 / (1+3*(4.74*10^-4*((ust(g)/ws(g,d))^1.5*Fr(g)-0.015)^1.18)); %equation 26b
        %fun = @(z) (ust(g) * Ca_dL(g,d) / Kappa) * ( ((h(g)-z)./z) * (a_dL(g)/(h(g)-a_dL(g))) ).^Rouse(g,d) .* log(z/(0.033*ks));
        %qs_dL(g,d) = integral(fun,a_dL(g),h(g)); %de Leeuw 2020

        % de Leeuw et al 2020 for total transport calculations
        ub(g) = 0.6 * u(g); %de Leeuw et al 2020 eq 20
        hb(g,d) = min(0.6 * (Fr(g) * (D(d)/h(g))^2)^0.3 * h(g),h(g)); %de Leeuw et al 2020 eq 21
        Cb(g,d) = qb_FB(g,d) / (hb(g,d) * ub(g)); %de Leeuw et al 2020 eq 17
        fun = @(z) (ust(g) * Cb(g,d) / Kappa) * ( ((h(g)-z)./z) * (hb(g,d)/(h(g)-hb(g,d))) ).^Rouse(g,d) .* log(z/(0.033*ks)); %when hb>a, assume rouse profile above Cb with Cb as reference at hb
        qs_dL2(g,d) = integral(fun,hb(g,d),h(g)); %integrating between Hb and h is the same as integrating between Hb-a and a-h and adding those. So this way of writing is simpler.

    end
end
clear g d ws Kappa ust
%clear ub hb a_dL Cb Ca_dL2 fun qs_dL2_a qs_dL2_b
%clear u Fr Rouse ks h

%% multiply sediment classes with bed availability from lognormal sediment curve
for g = 1:Lg
    for d = 1:LD
        qb_FB_mix(g,d) = qb_FB(g,d).*bins_y_perc(d)/100;
        qs_dL_mix(g,d) = qs_dL2(g,d).*bins_y_perc(d)/100;
        qt_EH_mix(g,d) = qt_EH(g,d).*bins_y_perc(d)/100;
    end
end
qt_dL_mix = qb_FB_mix + qs_dL_mix;

qb_FB_mixsum = sum(qb_FB_mix,2)*W; %[m3/s]
qs_dL_mixsum = sum(qs_dL_mix,2)*W; %[m3/s]
qt_dL_mixsum = sum(qt_dL_mix,2)*W; %[m3/s]
qt_EH_mixsum = sum(qt_EH_mix,2)*W; %[m3/s]

clear g d

%% Plot settings
for a = 1
    pl.width = 18; %[cm] plot width
    pl.height = 9.5; %[cm] plot height
    pl.line = 2; %line thickness
    pl.line_ax = 0.75; %line thickness axes
    pl.fsz = 7; %font size axes
    pl.fsz2 = 8.5; %font size labels
    CMAPOBJ = clrmap('read','Gravity.clrmap');
    clr1 = clrmap(CMAPOBJ,Lg); %colourmap gravity, red-blue
    pl.vertical = 1; %number of vertical subplots
    pl.horizantal = 2; %number of horizontal subplots
    pl.lmarge = 0.065; %left margin
    pl.rmarge = 0.025; %right margin
    pl.bmarge = 0.23; %bottom margin
    pl.tmarge = 0.02; %top margin
    pl.intv = 0.065; %vertical space between graphs
    pl.inth = 0.065; %horizontal space between graphs
    pl.wt = (1 - pl.lmarge - pl.rmarge - ((pl.horizantal-1)*pl.inth)) / pl.horizantal;
    pl.ht = (1 - pl.bmarge - pl.tmarge - ((pl.vertical-1)*pl.intv)) / pl.vertical;
    y = 0.94; %placing subplot letter, y
    xd = 10^(log10(D(1))+( log10(D(end))-log10(D(1)) )*0.04); %placing subplot letter, x grain size
    xg = grav(1) + (grav(end)-grav(1))*0.04;  %placing subplot letter, x gravity
end

%% Plot
for a = 1
    close all
    clc
    f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[120 5 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s1(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot(grav,qt_dL_mixsum,'-','color',[0 0 0],'linewidth',pl.line)
    plot(grav,qt_EH_mixsum,':','color',[0 0 0 0.5],'linewidth',pl.line)
    plot([1.351 1.352],[10^-1 10^6],'--','color',[0 0.348 0.15])
    plot([3.75 3.75],[10^-1 10^6],'--','color',clr1(1,:))
    plot([8.87 8.87],[10^-1 10^6],'--','color',[0.298 0.2 0])
    plot([9.81 9.81],[10^-1 10^6],'--','color',clr1(end,:))
    set(gca,'Yscale','log');
    %set(gca,'Xscale','log');
    xlim([grav(1) grav(end)]);
    ylim([10^0 10^1.5]);
    l(1) = xlabel('Gravity acceleration (g) [m/s^2]');
    l(2) = ylabel('Total transport [m^3/s]');
    t(1) = text(xg,10^(0+(1.5-0)*y),' a');
    t(2) = text(1.35+0.25,10^(0+(1.5-0)*0.025),'Titan','color',[0 0.348 0.15],'HorizontalAlignment', 'left');
    t(3) = text(3.75+0.25,10^(0+(1.5-0)*0.025),'Mars','color',[clr1(1,:)],'HorizontalAlignment', 'left');
    t(4) = text(8.87+0.25,10^(0+(1.5-0)*0.025),'Venus','color',[0.298 0.2 0],'HorizontalAlignment', 'left');
    t(5) = text(9.81+0.25,10^(0+(1.5-0)*0.025),'Earth','color',[clr1(end,:)],'HorizontalAlignment', 'left');
    set(t(2:5),'Rotation',90);
    %t(7) = text(grav(end),1.5*10^1,'W = 200 m + S = 0.001 m/m + Q = 2000 m^3/s   ','HorizontalAlignment', 'right');

    s1(2) = axes('Position',[pl.lmarge+pl.wt+pl.inth 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot(D,W*qt_dL_mix(gselect(1),:)./qt_dL_mixsum(gselect(1)).*100,'-','color',clr1(1,:),'linewidth',pl.line)
    plot(D,W*qt_dL_mix(gselect(2),:)./qt_dL_mixsum(gselect(2)).*100,'-','color',clr1(end,:),'linewidth',pl.line)
    plot(D,W*qt_EH_mix(gselect(1),:)./qt_EH_mixsum(gselect(1)).*100,':','color',[clr1(1,:) 0.5],'linewidth',pl.line)
    plot(D,W*qt_EH_mix(gselect(2),:)./qt_EH_mixsum(gselect(2)).*100,':','color',[clr1(end,:) 0.5],'linewidth',pl.line)

    plot(D, bins_y_perc,'--','color',[0 1 0])
    set(gca,'Xscale','log');
    xlim([D(1) D(end)]);
    ylim([0 60]);
    l(3) = xlabel('Sediment class');
    l(4) = ylabel('% of total');
    t(6) = text(xd,60*y,' b');
    set(gca,'XTick',D,'XtickLabel',bins_txt);

    leg1 = legend(s1(1),'de Leeuw et al. 2020','Engelund and Hansen 1967','location','northeast');
    leg2 = legend(s1(2),'Mars (de Leeuw et al. 2020)','Earth (de Leeuw et al. 2020)','Mars (Engelund and Hansen 1967)','Earth (Engelund and Hansen 1967)','Input sediment distribution','location','northeast');

    xtickangle(s1(2),90);
    set(s1,'box','on','Layer','top', ...
        'XMinorTick','on','YMinorTick','on', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);

    set(gcf,'renderer','painters');
    print(f1,'-dpng',[output '/FlowOnMars6_TransportMixture'],'-r400');
    print(f1,'-dpdf',[output '/FlowOnMars6_TransportMixture'],'-r400');
end

%% Plot Mike
