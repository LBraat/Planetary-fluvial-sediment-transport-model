%% Flow on Mars 3 - Sediment
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

%% Critical threshold of motion Q_mars = Q_earth
i = 1;
acc = 1e-5;
for d = 1:LD
    for g = 1:Lg
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
clear tau_Q_ini ust_Q_ini Rest_ini uprmsb_ust Pt B uy_ust ub_ust uprmsb_ub K shield1 minimise shieldini

%% Critical threshold of motion h_mars = h_earth
i = 1;
for d = 1:LD
    for g = 1:Lg
        %Zanke 2003 iterated
        shieldini = shield(g,d,i);
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
        shield_cr_Q(g,d) = shieldini;
        tau_cr_Q(g,d) = shield_cr_Q(g,d) * (rhos-rho)*grav(g)*D50(d);
        ust_cr_Q(g,d) = (tau_cr_Q(g,d) / rho)^0.5;
        lambda_cr_Q(g,d) = ust_cr_Q(g,d) / ws(g,d);
    end
end
clear tau_Q_ini ust_Q_ini Rest_ini uprmsb_ust Pt B uy_ust ub_ust uprmsb_ub K shield2 minimise shieldini

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
    for g = 1:Lg
        plot(D50,ws(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-6 10^1]);
    l(1) = ylabel('Settling velocity (w_s) [m/s]');
    %xlabel('Grain size (D_{50}) [m]');
    t(1) = text(xd,10^(-6+(1--6)*y),' a');
    
    s1(2) = axes('Position',[pl.lmarge+pl.wt+pl.inth 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        count = 1;
        for i = Qselect
            plot(D50,shield_Q(g,:,i),'-','color',[clr1(g,:) 1-(1/length(Qselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
    end
    for g = 1:Lg
        plot(D50,shield_cr_Q(g,:),':','color',clr1(g,:),'linewidth',pl.line);
    end
    for g = 1:Lg
        plot(D50,(ws(g,:).^2*rho)./((rhos-rho)*grav(g)*D50),'--','color',clr1(g,:),'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-3 10^3]);
    l(2) = ylabel('Shields parameter (\theta) [-]');
    %xlabel('Grain size (D_{50}) [m]');
    t(2) = text(xd,10^(-3+(3--3)*y),' b');
    
    s1(3) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:length(grav)
        count = 1;
        for i = Qselect
            plot(D50,lambda_Q(g,:,i),'-','color',[clr1(g,:) 1-(1/length(Qselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
        plot(D50,lambda_cr_Q(g,:),':','color',clr1(g,:),'linewidth',pl.line);
    end
    plot([D50(1) D50(end)],[1 1],'--k','linewidth',pl.line)
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-1 10^4]);
    l(3) = ylabel('Movability number (k) [-]');
    l(5) = xlabel('Grain size (D_{50}) [m]');
    t(3) = text(xd,10^(-1+(4--1)*y),' c');
    
    s1(4) = axes('Position',[pl.lmarge+pl.wt+pl.inth pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:length(grav)
        count = 1;
        for i = Qselect
            plot(D50,AdvL_Q(g,:,i),'-','color',[clr1(g,:) 1-(1/length(Qselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^0 10^7]);
    l(4) = ylabel('Advection length (A) [m]');
    l(6) = xlabel('Grain size (D_{50}) [m]');
    t(4) = text(xd,10^(0+(7-0)*y),' d');
    
    if length(Qselect) == 1
    legend(s1(2),'Mars Q=2000 m^3/s','Earth Q=2000 m^3/s','Mars motion threshold','Earth motion threshold','Mars suspension threshold','Earth suspension threshold', ...
        'location','northeast','position',[0.25 0.65 0.2 0.05])
    else
    legend(s1(3),'Mars Q=500 m^3/s','Mars Q=15000 m^3/s','Mars motion threshold','Earth Q=500 m^3/s','Mars Q=15000 m^3/s','Earth motion threshold','Suspension threshold', ...
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
        plot(D50,ws(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-6 10^1]);
    l(1) = ylabel('Settling velocity (w_s) [m/s]');
    %xlabel('Grain size (D_{50}) [m]');
    t(1) = text(xd,10^(-6+(1--6)*y),' a');
    
    s1(2) = axes('Position',[pl.lmarge+pl.wt+pl.inth 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        count = 1;
        for i = hselect
            plot(D50,shield(g,:,i),'-','color',[clr1(g,:) 1-(1/length(hselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
    end
    for g = 1:Lg
        plot(D50,shield_cr(g,:),':','color',clr1(g,:),'linewidth',pl.line);
    end
    for g = 1:Lg
        plot(D50,(ws(g,:).^2*rho)./((rhos-rho)*grav(g)*D50),'--','color',clr1(g,:),'linewidth',pl.line);
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-3 10^3]);
    l(2) = ylabel('Shields parameter (\theta) [-]');
    %xlabel('Grain size (D_{50}) [m]');
    t(2) = text(xd,10^(-3+(3--3)*y),' b');
    
    s1(3) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:length(grav)
        count = 1;
        for i = hselect
            plot(D50,lambda(g,:,i),'-','color',[clr1(g,:) 1-(1/length(hselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
        plot(D50,lambda_cr(g,:),':','color',clr1(g,:),'linewidth',pl.line);
    end
    plot([D50(1) D50(end)],[1 1],'--k','linewidth',pl.line)
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^-1 10^4]);
    l(3) = ylabel('Movability number (k) [-]');
    l(5) = xlabel('Grain size (D_{50}) [m]');
    t(3) = text(xd,10^(-1+(4--1)*y),' c');
    
    s1(4) = axes('Position',[pl.lmarge+pl.wt+pl.inth pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:length(grav)
        count = 1;
        for i = hselect
            plot(D50,AdvL(g,:,i),'-','color',[clr1(g,:) 1-(1/length(hselect))*(count-1)],'linewidth',pl.line);
            count = count + 1;
        end
    end
    set(gca,'Xscale','log');
    set(gca,'Yscale','log');
    xlim([D50(1) D50(end)]);
    ylim([10^0 10^7]);
    l(4) = ylabel('Advection length (A) [m]');
    l(6) = xlabel('Grain size (D_{50}) [m]');
    t(4) = text(xd,10^(0+(7-0)*y),' d');
    
    if length(hselect) == 1
    legend(s1(2),'Mars h=5 m','Earth h=5 m','Mars motion threshold','Earth motion threshold','Mars suspension threshold','Earth suspension threshold', ...
        'location','northeast','position',[0.25 0.65 0.2 0.05])
    else
    legend(s1(3),'Mars h=3 m','Mars h=15 m','Mars motion threshold','Earth h=3 m','Earth h=15 m','Earth motion threshold','Suspension threshold', ...
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
    print(f1,'-dpng',[output '/FlowOnMars3_Sediment_Dh'],'-r400');
    print(f1,'-dpdf',[output '/FlowOnMars3_Sediment_Dh'],'-r400');
end
