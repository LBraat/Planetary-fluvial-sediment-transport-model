%% Flow on Mars 3 - Sediment parameters
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
    hselect = 7; %
    Q = 500:500:15000; %[m3/s] Discharge
    Qselect = 4; %
    W = 200; %[m] Channel width
    grav = [3.7 9.8]; %[m/s2] Gravitational acceleration
    S = 0.001; %[m/m] Channel slope
    TC = 4; %[degrees C] Temperature
    rho = 1000; %[kg/m3] Water density

    % Sediment parameters
    rhos = 2900; %[kg/m3] Sediment density
    D = logspace(log10(63*10^-6),0,50); %[m] Grain size vector
    CDryB = 1750; %[kg/m3] Dry bed density
    C1 = 20; %18-24, smooth to rough sphere, affects small diameters, Ferguson and Church 2004
    C2 = 1; %0.4-1.2, smooth to rough sphere, affects big diameters, Ferguson and Church 2004
    phi = 30; %angle of repose in degrees

    % Derived input parameters
    TK = TC+273.15; %[degrees K] Temperature
    mu = 1.856e-14*exp(4209/TK+0.04527*TK-3.376e-5*TK^2); %[Pa | Ns/m2] Dynamic viscosity
    v = mu/rho; %[m2/s] Viscosity
    ks = 2.5*D; %[m] Nikurandse roughness height

    R = (rhos-rho)/rho; %[-] Relative density
    n = 1-CDryB/rhos; %[-] Porosity

    Lh = length(h);
    LQ = length(Q);
    Lg = length(grav);
    LD = length(D);
end
clear TC TK mu CDryB

%% allocate
for a = 1
    h_Q = NaN(Lg,LQ,LD);
    Rw_Q = NaN(Lg,LQ,LD);
    f_Q = NaN(Lg,LQ,LD);
    u_Q = NaN(Lg,LQ,LD);
    Fr_Q = NaN(Lg,LQ,LD);
    Re_Q = NaN(Lg,LQ,LD);
    tau_Q = NaN(Lg,LQ,LD);
    ust_Q = NaN(Lg,LQ,LD);
    %     lamlyr_Q = NaN(Lg,LQ,LD);

    Rw = NaN(1,Lh);
    f = NaN(Lg,Lh,LD);
    u = NaN(Lg,Lh,LD);
    Q_h = NaN(Lg,Lh,LD);
    Fr = NaN(Lg,Lh,LD);
    Re = NaN(Lg,Lh,LD);
    tau = NaN(Lg,Lh);
    ust = NaN(Lg,Lh);
    %     lamlyr = NaN(Lg,Lh);

    ws = NaN(Lg,LD);
    Rep = NaN(Lg,LD);
    Dst = NaN(Lg,LD);
    Rest_Q = NaN(Lg,LQ,LD);
    shield_Q = NaN(Lg,LQ,LD);
    lambda_Q = NaN(Lg,LQ,LD);
    AdvL_Q = NaN(Lg,LQ,LD);
    Rest = NaN(Lg,Lh,LD);
    shield = NaN(Lg,Lh,LD);
    lambda = NaN(Lg,Lh,LD);
    AdvL = NaN(Lg,Lh,LD);

    shield_cr = NaN(Lg,LD);
    tau_cr = NaN(Lg,LD);
    ust_cr = NaN(Lg,LD);
    lambda_cr = NaN(Lg,LD);
    shield_cr_Q = NaN(Lg,LD);
    tau_cr_Q = NaN(Lg,LD);
    ust_cr_Q = NaN(Lg,LD);
    lambda_cr_Q = NaN(Lg,LD);
end

%% Hydro parameters Q_mars = Q_earth
h_guess = 3;
for d = 1:LD
    for i = 1:LQ
        for g = 1:Lg
            h_Q(g,i,d) = h_guess+1;
            while abs(h_guess-h_Q(g,i,d))>0.0001
                h_guess = h_Q(g,i,d);
                Rw_Q(g,i,d) = (h_guess*W)/(h_guess+h_guess+W); %[m] Hydraulic radius
                f_Q(g,i,d) = 8/(5.75*log10(12*h_guess/ks(d)))^2; %[-] friction factor
                u_Q(g,i,d) = (8*grav(g)*Rw_Q(g,i,d)*S/f_Q(g,i,d))^0.5; %[m/s] Velocity
                h_Q(g,i,d) = Q(i)/(W*u_Q(g,i,d)); %[m] Water depth
            end
            tau_Q(g,i,d) = rho*grav(g)*Rw_Q(g,i,d)*S; %[N/m2] Bed shear stress
            ust_Q(g,i,d) = (tau_Q(g,i,d)/rho)^0.5; %[m/s] Shear velocity
            Fr_Q(g,i,d) = u_Q(g,i,d)/(grav(g)*h_Q(g,i,d))^0.5; %[-] Froude number
            Re_Q(g,i,d) = u_Q(g,i,d)*h_Q(g,i,d)/v; %[-] Reynolds number
            %lamlyr_Q(g,i,d) = 11.63*v/ust_Q(g,i,d); %[m] Laminar sublayer thickness
        end
    end
end

%% Hydro parameters h_mars = h_earth
for i = 1:Lh
    Rw(i) = (h(i)*W)/(h(i)+h(i)+W); %[m] Hydraulic radius
    for g = 1:Lg
        for d = 1:LD
            f(g,i,d) = 8/(5.74*log10(12.2*h(i)/ks(d)))^2; %[-] friction factor
            u(g,i,d) = (8*grav(g)*Rw(i)*S/f(g,i,d))^0.5; %[m/s] Velocity
            Q_h(g,i,d) = h(i)*W*u(g,i,d); %[m3/s] Discharge
            Fr(g,i,d) = u(g,i,d)/(grav(g)*h(i))^0.5; %[-] Froude number
            Re(g,i,d) = u(g,i,d)*h(i)/v; %[-] Reynolds number
            tau(g,i) = rho*grav(g)*Rw(i)*S; %[N/m2] Bed shear stress
            ust(g,i) = (tau(g,i)/rho)^0.5; %[m/s] Shear velocity
            %lamlyr(g) = 11.63*v/ust(g); %[m] Laminar sublayer thickness
        end
    end
end

%% Sediment parameters
for d = 1:LD
    for g = 1:Lg
        ws(g,d) = R*grav(g)*D(d)^2/(C1*v+(0.75*C2*R*grav(g)*D(d)^3)^0.5); %[m/s] Settling velocity (Ferguson and Church 2004)
        Rep(g,d) = D(d)^1.5*(R*grav(g))^0.5/v; %[-] Particle Reynolds number (Klein05,Leeuw20,Miedema), >70 turbulent <3.5 laminar
        Dst(g,d) = D(d)*(R*grav(g)/v^2)^(1/3); %[-] Dimensionless particle parameter/Bonnefille number

        % Q_mars = Q_earth
        for i = 1:LQ
            Rest_Q(g,i,d) = D(d)*ust_Q(g,i,d)/v; %[-] Particle Reynolds number (vR84,Leeuw20), Reynolds shear velocity number (Klein05), Boundary Reynolds number (Miedema10)
            shield_Q(g,i,d) = tau_Q(g,i,d)/((rhos-rho)*grav(g)*D(d)); %[-] Shields parameter/Particle mobility parameter
            lambda_Q(g,i,d) = ust_Q(g,i,d)/ws(g,d); %[-] Movability number (Liu 1958)
            AdvL_Q(g,i,d) = u_Q(g,i,d)*h_Q(g,i,d)/ws(g,d); %[m] Advection length
        end
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

%% Critical threshold of motion Q_mars = Q_earth
i = Qselect;
acc = 1e-5;
for d = 1:LD
    for g = 1:Lg
        %Zanke 2003 iterated
        shieldini = shield_Q(g,i,d);
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
clear i shieldini minimise tau_Q_ini ust_Q_ini Rest_ini uprmsb_ust Pt B uy_ust ub_ust uprmsb_ub K shield1

%% Critical threshold of motion h_mars = h_earth
i = hselect;
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
            shield2 = (1-n)*tan(deg2rad(phi/1.5))*K / ...
                ((1+1.8*uprmsb_ub)^2 * (1+0.4*(1.8*uprmsb_ust)^2*tan(deg2rad(phi/1.5))*K));
            minimise = abs(shield2 - shieldini);
            shieldini = shield2;
        end
        shield_cr_Q(g,d) = shieldini;
        tau_cr_Q(g,d) = shield_cr_Q(g,d) * (rhos-rho)*grav(g)*D(d);
        ust_cr_Q(g,d) = (tau_cr_Q(g,d) / rho)^0.5;
        lambda_cr_Q(g,d) = ust_cr_Q(g,d) / ws(g,d);
    end
end
clear i acc shieldini minimise tau_Q_ini ust_Q_ini Rest_ini uprmsb_ust Pt B uy_ust ub_ust uprmsb_ub K shield2 phi

%% Plot settings
for a = 1
    pl.width = 18; %[cm] plot width
    pl.height = 16; %[cm] plot height
    pl.line = 1; %line thickness
    pl.line_ax = 0.75; %line thickness axes
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
    xd = 10^(log10(D(1))+( log10(D(end))-log10(D(1)) )*0.04); %placing subplot letter, x
end

%% plot Q_mars = Q_earth
for a = 1
    close all; clc;
    f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[10 1 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s1(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(D,ws(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D(1) D(end)]);
    ylim([10^-3 10^1]);
    l(1) = ylabel('Settling velocity (w_s) [m/s]');
    xlabel('Grain size (D) [m]');
    t(1) = text(xd,10^(-3+(1--3)*y),' a');

    s1(2) = axes('Position',[pl.lmarge+pl.wt+pl.inth 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        count = 1;
        for i = Qselect
            plot(D,squeeze(shield_Q(g,i,:)),'-','color',[clr1(g,:) 1-(1/length(Qselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
    end
    for g = 1:Lg
        plot(D,shield_cr_Q(g,:),':','color',clr1(g,:),'linewidth',pl.line);
    end
    for g = 1:Lg
        plot(D,(ws(g,:).^2*rho)./((rhos-rho)*grav(g)*D),'--','color',clr1(g,:),'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D(1) D(end)]);
    ylim([10^-3 10^2]);
    l(2) = ylabel('Shields number (\theta) [-]');
    xlabel('Grain size (D) [m]');
    t(2) = text(xd,10^(-3+(2--3)*y),' b');

    s1(3) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:length(grav)
        count = 1;
        for i = Qselect
            plot(D,squeeze(lambda_Q(g,i,:)),'-','color',[clr1(g,:) 1-(1/length(Qselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
        plot(D,lambda_cr_Q(g,:),':','color',clr1(g,:),'linewidth',pl.line);
    end
    plot([D(1) D(end)],[1 1],'--k','linewidth',pl.line)
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D(1) D(end)]);
    ylim([10^-1 10^2]);
    l(3) = ylabel('Movability number (k) [-]');
    l(5) = xlabel('Grain size (D) [m]');
    t(3) = text(xd,10^(-1+(2--1)*y),' c');

    s1(4) = axes('Position',[pl.lmarge+pl.wt+pl.inth pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:length(grav)
        count = 1;
        for i = Qselect
            plot(D,squeeze(AdvL_Q(g,i,:)),'-','color',[clr1(g,:) 1-(1/length(Qselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D(1) D(end)]);
    ylim([10^0 10^4]);
    l(4) = ylabel('Advection length (L_A) [m]');
    l(6) = xlabel('Grain size (D) [m]');
    t(4) = text(xd,10^(0+(4-0)*y),' d');

    if length(Qselect) == 1
        legend(s1(2),'Mars Q=2000 m^3/s','Earth Q=2000 m^3/s','Mars motion threshold','Earth motion threshold','Mars suspension threshold','Earth suspension threshold', ...
            'location','northeast','position',[0.25 0.65 0.2 0.05])
    else
        legend(s1(3),'Mars Q=500 m^3/s','Mars Q=15000 m^3/s','Mars motion threshold','Earth Q=500 m^3/s','Mars Q=15000 m^3/s','Earth motion threshold','Suspension threshold', ...
            'location','northeast')
    end

    set(s1,'box','on','Layer','top', ...
        'XMinorTick','on','YMinorTick','on', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);

    set(gcf,'renderer','painters');
    print(f1,'-dpng',[output '/FlowOnMars3_Sediment_DQ'],'-r400');
    print(f1,'-dpdf',[output '/FlowOnMars3_Sediment_DQ'],'-r400');
end

%% plot h_mars = h_earth
for a = 1
    close all; clc;
    f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[10 1 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s1(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(D,ws(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D(1) D(end)]);
    ylim([10^-3 10^1]);
    l(1) = ylabel('Settling velocity (w_s) [m/s]');
    xlabel('Grain size (D) [m]');
    t(1) = text(xd,10^(-3+(1--3)*y),' a');

    s1(2) = axes('Position',[pl.lmarge+pl.wt+pl.inth 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        count = 1;
        for i = hselect
            plot(D,squeeze(shield(g,i,:)),'-','color',[clr1(g,:) 1-(1/length(hselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
    end
    for g = 1:Lg
        plot(D,shield_cr(g,:),':','color',clr1(g,:),'linewidth',pl.line);
    end
    for g = 1:Lg
        plot(D,(ws(g,:).^2*rho)./((rhos-rho)*grav(g)*D),'--','color',clr1(g,:),'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D(1) D(end)]);
    ylim([10^-3 10^2]);
    l(2) = ylabel('Shields number (\theta) [-]');
    xlabel('Grain size (D) [m]');
    t(2) = text(xd,10^(-3+(2--3)*y),' b');

    s1(3) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:length(grav)
        count = 1;
        for i = hselect
            plot(D,squeeze(lambda(g,i,:)),'-','color',[clr1(g,:) 1-(1/length(hselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
        plot(D,lambda_cr(g,:),':','color',clr1(g,:),'linewidth',pl.line);
    end
    plot([D(1) D(end)],[1 1],'--k','linewidth',pl.line)
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D(1) D(end)]);
    ylim([10^-1 10^2]);
    l(3) = ylabel('Movability number (k) [-]');
    l(5) = xlabel('Grain size (D) [m]');
    t(3) = text(xd,10^(-1+(2--1)*y),' c');

    s1(4) = axes('Position',[pl.lmarge+pl.wt+pl.inth pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:length(grav)
        count = 1;
        for i = hselect
            plot(D,squeeze(AdvL(g,i,:)),'-','color',[clr1(g,:) 1-(1/length(hselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D(1) D(end)]);
    ylim([10^0 10^4]);
    l(4) = ylabel('Advection length (L_A) [m]');
    l(6) = xlabel('Grain size (D) [m]');
    t(4) = text(xd,10^(0+(4-0)*y),' d');

    if length(hselect) == 1
        legend(s1(2),'Mars h=5 m','Earth h=5 m','Mars motion threshold','Earth motion threshold','Mars suspension threshold','Earth suspension threshold', ...
            'location','northeast','position',[0.25 0.65 0.2 0.05])
    else
        legend(s1(3),'Mars h=3 m','Mars h=15 m','Mars motion threshold','Earth h=3 m','Earth h=15 m','Earth motion threshold','Suspension threshold', ...
            'location','northeast')
    end

    set(s1,'box','on','Layer','top', ...
        'XMinorTick','on','YMinorTick','on', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);

    set(gcf,'renderer','painters');
    print(f1,'-dpng',[output '/FlowOnMars3_Sediment_Dh'],'-r400');
    print(f1,'-dpdf',[output '/FlowOnMars3_Sediment_Dh'],'-r400');
end
