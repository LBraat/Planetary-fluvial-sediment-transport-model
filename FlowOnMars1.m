%% Flow on Mars 1 - Hydrological parameters
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
    h = 0.5:0.5:15; %[m] Water depth
    Q = 250:500:15000; %[m3/s] Discharge
    Qselect = 4; %Discharge selection last plot
    W = 200; %[m] Channel width
    grav = 1:0.5:12; %[m/s2] Gravitational acceleration
    S = 0.001; %[m/m] Channel slope
    TC = 4; %[degrees C] Temperature
    rho = 1000; %[kg/m3] Water density
    gE = 9.8; %[m/s2] Gravity Earth
    gM = 3.7; %[m/s2] Gravity Mars
    D = 0.01; %[m] Grain size for bed roughness only

    % Derived input parameters
    TK = TC+273.15; %[degrees K] Temperature
    mu = 1.856e-14*exp(4209/TK+0.04527*TK-3.376e-5*TK^2); %[Pa | Ns/m2] Absolute/Dynamic viscosity
    v = mu/rho; %[m2/s] Kinematic Viscosity
    ks = 2.5*D; %[m] Nikurandse roughness height

    Lh = length(h);
    LQ = length(Q);
    Lg = length(grav);
end
clear TC TK mu

%% allocate
for a = 1
    h_Q = NaN(Lg,LQ);
    Rw_Q = NaN(Lg,LQ);
    f_Q = NaN(Lg,LQ);
    u_Q = NaN(Lg,LQ);
    Fr_Q = NaN(Lg,LQ);
    Re_Q = NaN(Lg,LQ);
    tau_Q = NaN(Lg,LQ);
    ust_Q = NaN(Lg,LQ);
    %lamlyr_Q = NaN(Lg,LQ);
    
    Rw = NaN(1,Lh);
    f = NaN(Lg,Lh);
    u = NaN(Lg,Lh);
    Q_h = NaN(Lg,Lh);
    Fr = NaN(Lg,Lh);
    Re = NaN(Lg,Lh);
    tau = NaN(Lg,Lh);
    ust = NaN(Lg,Lh);
    %lamlyr = NaN(Lg,Lh);
end

%% Hydro parameters Q_mars = Q_earth
h_guess = 3;
for i = 1:LQ
    for g = 1:Lg
        h_Q(g,i) = h_guess+1;
        while abs(h_guess-h_Q(g,i))>0.0001
            h_guess = h_Q(g,i);
            Rw_Q(g,i) = (h_guess*W)/(h_guess+h_guess+W); %[m] Hydraulic radius
            f_Q(g,i) = 8/(5.75*log10(12*h_guess/ks))^2; %[-] friction factor
            u_Q(g,i) = (8*grav(g)*Rw_Q(g,i)*S/f_Q(g,i))^0.5; %[m/s] Velocity
            h_Q(g,i) = Q(i)/(W*u_Q(g,i)); %[m] Water depth
        end
        tau_Q(g,i) = rho*grav(g)*Rw_Q(g,i)*S; %[N/m2] Bed shear stress
        ust_Q(g,i) = (tau_Q(g,i)/rho)^0.5; %[m/s] Shear velocity
        Fr_Q(g,i) = u_Q(g,i)/(grav(g)*h_Q(g,i))^0.5; %[-] Froude number
        Re_Q(g,i) = u_Q(g,i)*h_Q(g,i)/v; %[-] Reynolds number
        %lamlyr_Q(g) = 11.63*v/ust_Q(g); %[m] Laminar sublayer thickness
    end
end

%% Hydro parameters h_mars = h_earth
for i = 1:Lh
    Rw(i) = (h(i)*W)/(h(i)+h(i)+W); %[m] Hydraulic radius
    for g = 1:Lg
        f(g,i) = 8/(5.74*log10(12.2*h(i)/ks))^2; %[-] friction factor
        u(g,i) = (8*grav(g)*Rw(i)*S/f(g,i))^0.5; %[m/s] Velocity
        Q_h(g,i) = h(i)*W*u(g,i); %[m3/s] Discharge
        Fr(g,i) = u(g,i)/(grav(g)*h(i))^0.5; %[-] Froude number
        Re(g,i) = u(g,i)*h(i)/v; %[-] Reynolds number
        tau(g,i) = rho*grav(g)*Rw(i)*S; %[N/m2] Bed shear stress
        ust(g,i) = (tau(g,i)/rho)^0.5; %[m/s] Shear velocity
        %lamlyr(g) = 11.63*v/ust(g); %[m] Laminar sublayer thickness
    end
end

%% Plot settings
for a = 1
    pl.width = 18; %[cm] plot width
    pl.height = 9; %[cm] plot height
    pl.line = 1; %line thickness
    pl.line_ax = 0.75; %line thickness axes
    pl.fsz = 7; %font size axes
    pl.fsz2 = 8.5; %font size labels
    CMAPOBJ = clrmap('read','Gravity.clrmap'); 
    clr1 = clrmap(CMAPOBJ,Lg); %colourmap gravity, red-blue
    CMAPOBJ = clrmap('read','Discharge.clrmap'); 
    clr2 = clrmap(CMAPOBJ,Lh); %colourmap water depth, licht-dark blue
    clr3 = clrmap(CMAPOBJ,LQ); %colourmap discharge, licht-dark blue
    clrM = [1 0.804 0.804]; %red line Mars
    clrE = [0.8 0.8 0.9804]; %blue line Earth
    pl.vertical = 2; %number of vertical subplots
    pl.horizantal = 4; %number of horizontal subplots
    pl.lmarge = 0.07; %left margin
    pl.rmarge = 0.025; %right margin
    pl.bmarge = 0.10; %bottom margin
    pl.tmarge = 0.13; %top margin
    pl.intv = 0.065; %vertical space between graphs
    pl.inth = 0.066; %horizontal space between graphs
    pl.wt = (1 - pl.lmarge - pl.rmarge - ((pl.horizantal-1)*pl.inth)) / pl.horizantal;
    pl.ht = (1 - pl.bmarge - pl.tmarge - ((pl.vertical-1)*pl.intv)) / pl.vertical;
    y = 0.94; %placing subplot letter, y
    xq = Q(1) + (Q(end)-Q(1))*0.01; %placing subplot letter, x-axis discharge
    xh = h(1) + (h(end)-h(1))*0.01; %placing subplot letter x-axis water depth
    xg = grav(1) + (grav(end)-grav(1))*0.01; %placing subplot letter, x-axis gravity
end

%% plot x=h h_mars = h_earth
for a = 1
    if a == 1
    f1 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[10 5 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s1(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h,Q_h(g,:),'color',clr1(g,:),'LineWidth',pl.line);
    end
    l(1) = ylabel('Discharge (Q) [m^3/s]');
    axis([h(1) h(end) 0 15000]);
    t(1) = text(xh,15000*y,' a');

    s1(2) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h,u(g,:),'color',clr1(g,:),'LineWidth',pl.line);
    end
    l(2) = ylabel('Velocity (u) [m/s]');
    l(3) = xlabel('Water depth (h) [m]');
    axis([h(1) h(end) 0 7]);
    t(2) = text(xh,7*y,' e');

    s1(3) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h,Rw,'color',clr1(g,:),'LineWidth',pl.line);
    end
    l(4) = ylabel('Hydraulic radius (R_w) [m]');
    axis([h(1) h(end) 0 15]);
    t(3) = text(xh,15*y,' b');

    s1(4) = axes('Position',[pl.lmarge+pl.inth+pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h,tau(g,:),'color',clr1(g,:),'LineWidth',pl.line);
    end
    l(5) = ylabel('Shear stress (\tau) [N/m^2]');
    l(6) = xlabel('Water depth (h) [m]');
    axis([h(1) h(end) 0 100]);
    t(4) = text(xh,100*y,' f');

    s1(5) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h,f(g,:),'color',clr1(g,:),'LineWidth',pl.line);
    end
    l(7) = ylabel({'Friction factor (f) [-]'});
    axis([h(1) h(end) 0 0.035]);
    t(5) = text(xh,0.035*y,' c');

    s1(6) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h,ust(g,:),'color',clr1(g,:),'LineWidth',pl.line);
    end
    l(8) = xlabel('Water depth (h) [m]');
    l(9) = ylabel('Shear velocity (u_*) [m/s]');
    axis([h(1) h(end) 0 0.3]);
    t(6) = text(xh,0.3*y,' g');

    s1(7) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h,Fr(g,:),'color',clr1(g,:),'LineWidth',pl.line);
    end
    l(10) = ylabel('Froude number (Fr) [-]');
    axis([h(1) h(end) 0 1.2]);
    plot([h(1) h(end)],[1 1],'--','color','k','LineWidth',pl.line);
    t(9) = text(h(end),0.93,'subcritical ','HorizontalAlignment','right');
    t(10) = text(h(end),1.09,'supercritical ','HorizontalAlignment','right');
    t(7) = text(xh,1.2*y,' d');

    s1(8) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h,Re(g,:),'color',clr1(g,:),'LineWidth',pl.line);
    end
    l(11) = xlabel('Water depth (h) [m]');
    l(12) = ylabel('Reynolds number (Re) [-]');
    axis([h(1) h(end) 150 10e7]);
    plot([h(1) h(end)],[500 500],'--','color','k','LineWidth',pl.line);
    t(11) = text(h(end),300,'laminar ','HorizontalAlignment','right');
    t(12) = text(h(end),1100,'tubulent ','HorizontalAlignment','right');
    set(gca,'Yscale','log');
    t(8) = text(xh,10^(log10(150)+(8-log10(150))*y),' h');

    colormap(clr1)
    cb = colorbar('north','position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge+0.015 pl.wt 0.02],'Ticks',[grav(1) grav(end)]);
    cb.Label.String ='Gravity [m/s^2]';
    cb.Label.FontSize = pl.fsz2;
    clim([grav(1) grav(end)]);

    set(s1,'box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);

    set(gcf,'renderer','painters');
    print(f1,'-dpng',[output '/FlowOnMars1_hh'],'-r400');
    print(f1,'-dpdf',[output '/FlowOnMars1_hh'],'-r400');
    end
end

%% plot x=h Q_mars = Q_earth
for a = 1
    if a == 1
    f2 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[10 5 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s2(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h_Q(g,:),Q,'color',clr1(g,:),'linewidth',pl.line);
    end
    l(1) = ylabel('Discharge (Q) [m^3/s]');
    axis([h(1) h(end) 0 15000]);
    t(1) = text(xh,15000*y,' a');

    s2(2) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h_Q(g,:),u_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(2) = xlabel('Water depth (h) [m]');
    l(3) = ylabel('Velocity (u) [m/s]');
    axis([h(1) h(end) 0 7]);
    t(2) = text(xh,7*y,' e');

    s2(3) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h_Q(g,:),Rw_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(4) = ylabel('Hydraulic radius (R_w) [m]');
    axis([h(1) h(end) 0 15]);
    t(3) = text(xh,15*y,' b');

    s2(4) = axes('Position',[pl.lmarge+pl.inth+pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h_Q(g,:),tau_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(5) = ylabel('Shear stress (\tau) [N/m^2]');
    l(6) = xlabel('Water depth (h) [m]');
    axis([h(1) h(end) 0 100]);
    t(4) = text(xh,100*y,' f');

    s2(5) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h_Q(g,:),f_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(7) = ylabel({'Friction factor (f) [-]'});
    axis([h(1) h(end) 0 0.035]);
    t(5) = text(xh,0.035*y,' c');

    s2(6) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h_Q(g,:),ust_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(8) = xlabel('Water depth (h) [m]');
    l(9) = ylabel('Shear velocity (u_*) [m/s]');
    axis([h(1) h(end) 0 0.3]);
    t(6) = text(xh,0.3*y,' g');

    s2(7) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h_Q(g,:),Fr_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(10) = ylabel('Froude number (Fr) [-]');
    axis([h(1) h(end) 0 1.2]);
    plot([h(1) h(end)],[1 1],'--','color','k','linewidth',pl.line);
    t(9) = text(h(end),0.93,'subcritical ','HorizontalAlignment','right');
    t(10) = text(h(end),1.09,'supercritical ','HorizontalAlignment','right');
    t(7) = text(xh,1.2*y,' d');

    s2(8) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(h_Q(g,:),Re_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(11) = xlabel('Water depth (h) [m]');
    l(12) = ylabel('Reynolds number (Re) [-]');
    axis([h(1) h(end) 150 10e7]);
    plot([h(1) h(end)],[500 500],'--','color','k','linewidth',pl.line);
    t(11) = text(h(end),300,'laminar ','HorizontalAlignment','right');
    t(12) = text(h(end),1100,'tubulent ','HorizontalAlignment','right');
    set(gca,'Yscale','log');
    t(8) = text(xh,10^(log10(150)+(8-log10(150))*y),' h');

    colormap(clr1)
    cb = colorbar('north','position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge+0.015 pl.wt 0.02],'Ticks',[grav(1) grav(end)]);
    cb.Label.String ='Gravity [m/s^2]';
    cb.Label.FontSize = pl.fsz2;
    clim([grav(1) grav(end)]);

    set(s2,'box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);

    set(gcf,'renderer','painters');
    print(f2,'-dpng',[output '/FlowOnMars1_hQ'],'-r400');
    print(f2,'-dpdf',[output '/FlowOnMars1_hQ'],'-r400');
    end
end

%% plot x=Q h_mars = h_earth
for a = 1
    if a == 1
    f3 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[10 5 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s3(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q_h(g,:),h,'color',clr1(g,:),'linewidth',pl.line);
    end
    l(1) = ylabel('Water depth (h) [m]');
    axis([Q(1) Q(end) 0 15]);
    t(1) = text(xq,15*y,' a');

    s3(2) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q_h(g,:),u(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(2) = xlabel('Discharge (Q) [m^3/s]');
    l(3) = ylabel('Velocity (u) [m/s]');
    axis([Q(1) Q(end) 0 7]);
    t(2) = text(xq,7*y,' e');

    s3(3) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q_h(g,:),Rw,'color',clr1(g,:),'linewidth',pl.line);
    end
    l(4) = ylabel('Hydraulic radius (R_w) [m]');
    axis([Q(1) Q(end) 0 15]);
    t(3) = text(xq,15*y,' b');

    s3(4) = axes('Position',[pl.lmarge+pl.inth+pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q_h(g,:),tau(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(5) = xlabel('Discharge (Q) [m^3/s]');
    l(6) = ylabel('Shear stress (\tau) [N/m^2]');
    axis([Q(1) Q(end) 0 100]);
    t(4) = text(xq,100*y,' f');

    s3(5) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q_h(g,:),f(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(7) = ylabel([{'Friction factor (f) [-]'}]);
    axis([Q(1) Q(end) 0 0.035]);
    t(5) = text(xh,0.035*y,' c');

    s3(6) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q_h(g,:),ust(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(8) = xlabel('Discharge (Q) [m^3/s]');
    l(9) = ylabel('Shear velocity (u_*) [m/s]');
    axis([Q(1) Q(end) 0 0.3]);
    t(6) = text(xq,0.3*y,' g');

    s3(7) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q_h(g,:),Fr(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(10) = ylabel('Froude number (Fr) [-]');
    axis([Q(1) Q(end) 0 1.2]);
    plot([Q(1) Q(end)],[1 1],'--','color','k','linewidth',pl.line);
    t(9) = text(Q(end),0.93,'subcritical ','HorizontalAlignment','right');
    t(10) = text(Q(end),1.09,'supercritical ','HorizontalAlignment','right');
    t(7) = text(xq,1.2*y,' d');

    s3(8) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q_h(g,:),Re(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(11) = xlabel('Discharge (Q) [m^3/s]');
    l(12) = ylabel('Reynolds number (Re) [-]');
    axis([Q(1) Q(end) 150 10e7]);
    plot([Q(1) Q(end)],[500 500],'--','color','k','linewidth',pl.line);
    t(11) = text(Q(end),300,'laminar ','HorizontalAlignment','right');
    t(12) = text(Q(end),1100,'tubulent ','HorizontalAlignment','right');
    set(gca,'Yscale','log');
    t(8) = text(xq,10^(log10(150)+(8-log10(150))*y),' h');

    colormap(clr1)
    cb = colorbar('north','position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge+0.015 pl.wt 0.02],'Ticks',[grav(1) grav(end)]);
    cb.Label.String ='Gravity [m/s^2]';
    cb.Label.FontSize = pl.fsz2;
    clim([grav(1) grav(end)]);

    set(s3,'box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);

    set(gcf,'renderer','painters');
    print(f3,'-dpng',[output '/FlowOnMars1_Qh'],'-r400');
    print(f3,'-dpdf',[output '/FlowOnMars1_Qh'],'-r400');
    end
end

%% plot x=Q Q_mars = Q_earth
for a = 1
    if a == 1
    f4 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[10 5 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s4(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q,h_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(1) = ylabel('Water depth (h) [m]');
    axis([Q(1) Q(end) 0 15]);
    t(1) = text(xq,15*y,' a');

    s4(2) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q,Rw_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(2) = ylabel('Hydraulic radius (R_w) [m]');
    axis([Q(1) Q(end) 0 15]);
    t(2) = text(xq,15*y,' b');

    s4(3) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q,f_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(3) = ylabel([{'Friction factor (f) [-]'}]);
    axis([Q(1) Q(end) 0 0.035]);
    t(3) = text(xq,0.035*y,' c');

    s4(4) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q,Fr_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(4) = ylabel('Froude number (Fr) [-]');
    axis([Q(1) Q(end) 0 1.2]);
    plot([Q(1) Q(end)],[1 1],'--','color','k','linewidth',pl.line);
    t(4) = text(Q(end),0.93,'subcritical ','HorizontalAlignment','right');
    t(5) = text(Q(end),1.09,'supercritical ','HorizontalAlignment','right');
    t(6) = text(xq,1.2*y,' d');

    s4(5) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q,u_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(5) = xlabel('Discharge (Q) [m^3/s]');
    l(6) = ylabel('Velocity (u) [m/s]');
    axis([Q(1) Q(end) 0 7]);
    t(7) = text(xq,7*y,' e');

    s4(6) = axes('Position',[pl.lmarge+pl.inth+pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q,tau_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(7) = xlabel('Discharge (Q) [m^3/s]');
    l(8) = ylabel('Shear stress (\tau) [N/m^2]');
    axis([Q(1) Q(end) 0 100]);
    t(8) = text(xq,100*y,' f');

    s4(7) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q,ust_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(9) = xlabel('Discharge (Q) [m^3/s]');
    l(10) = ylabel('Shear velocity (u_*) [m/s]');
    axis([Q(1) Q(end) 0 0.3]);
    t(9) = text(xq,0.3*y,' g');

    s4(8) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    for g = 1:Lg
        plot(Q,Re_Q(g,:),'color',clr1(g,:),'linewidth',pl.line);
    end
    l(11) = xlabel('Discharge (Q) [m^3/s]');
    l(12) = ylabel('Reynolds number (Re) [-]');
    axis([Q(1) Q(end) 150 10e7]);
    plot([Q(1) Q(end)],[500 500],'--','color','k','linewidth',pl.line);
    t(10) = text(Q(end),300,'laminar ','HorizontalAlignment','right');
    t(11) = text(Q(end),1100,'tubulent ','HorizontalAlignment','right');
    set(gca,'Yscale','log');
    t(12) = text(xq,10^(log10(150)+(8-log10(150))*y),' h');

    colormap(clr1)
    cb = colorbar('north','position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge+0.015 pl.wt 0.02],'Ticks',[grav(1) grav(end)]);
    cb.Label.String ='Gravity [m/s^2]';
    cb.Label.FontSize = pl.fsz2;
    clim([grav(1) grav(end)]);

    set(s4,'box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);

    set(gcf,'renderer','painters');
    print(f4,'-dpng',[output '/FlowOnMars1_QQ'],'-r400');
    print(f4,'-dpdf',[output '/FlowOnMars1_QQ'],'-r400');
    end
end

%% plot x=g Q_mars = Q_earth
for a = 1
    if a == 1
    f5 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[10 5 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s5(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:LQ
        plot(grav,h_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(1) = ylabel('Water depth (h) [m]');
    axis([grav(1) grav(end) 0 15]);
    t(1) = text(gM-0.5,0,' Mars','color',clrM,'rotation',90);
    t(2) = text(gE-0.5,0,' Earth','color',clrE,'rotation',90);
    t(3) = text(xg,15*y,' a');

    s5(2) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:LQ
        plot(grav,u_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(2) = xlabel('Gravity (g) [m/s^2]');
    l(3) = ylabel('Velocity (u) [m/s]');
    axis([grav(1) grav(end) 0 7]);
    t(4) = text(xg,7*y,' e');

    s5(3) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:LQ
        plot(grav,Rw_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(4) = ylabel('Hydraulic radius (R_w) [m]');
    axis([grav(1) grav(end) 0 15]);
    t(5) = text(xg,15*y,' b');

    s5(4) = axes('Position',[pl.lmarge+pl.inth+pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:LQ
        plot(grav,tau_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(5) = xlabel('Gravity (g) [m/s^2]');
    l(6) = ylabel('Shear stress (\tau) [N/m^2]');
    axis([grav(1) grav(end) 0 100]);
    t(6) = text(xg,100*y,' f');

    s5(5) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:LQ
        plot(grav,f_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(7) = ylabel({'Friction factor (f) [-]'});
    axis([grav(1) grav(end) 0 0.04]);
    t(7) = text(xg,0.04*y,' c');

    s5(6) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:LQ
        plot(grav,ust_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(8) = xlabel('Gravity (g) [m/s^2]');
    l(9) = ylabel('Shear velocity (u_*) [m/s]');
    axis([grav(1) grav(end) 0 0.3]);
    t(8) = text(xg,0.3*y,' g');

    s5(7) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:LQ
        plot(grav,Fr_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(10) = ylabel('Froude number (Fr) [-]');
    axis([grav(1) grav(end) 0 1.20]);
    plot([grav(1) grav(end)],[1 1],'--','color','k','linewidth',pl.line);
    t(9) = text(grav(end),0.93,'subcritical ','HorizontalAlignment','right');
    t(10) = text(grav(end),1.09,'supercritical ','HorizontalAlignment','right');
    t(11) = text(xg,1.2*y,' d');

    s5(8) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[1 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[1 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:LQ
        plot(grav,Re_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(11) = xlabel('Gravity (g) [m/s^2]');
    l(12) = ylabel('Reynolds number (Re) [-]');
    axis([grav(1) grav(end) 150 10e7]);
    plot([grav(1) grav(end)],[500 500],'--','color','k','linewidth',pl.line);
    t(12) = text(grav(end),300,'laminar ','HorizontalAlignment','right');
    t(13) = text(grav(end),1100,'tubulent ','HorizontalAlignment','right');
    set(gca,'Yscale','log');
    t(14) = text(xg,10^(log10(150)+(8-log10(150))*y),' h');

    colormap(clr3)
    cb = colorbar('north','position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge+0.015 pl.wt 0.02],'Ticks',[Q(1) Q(end)]);
    cb.Label.String ='Discharge [m^3/s]';
    cb.Label.FontSize = pl.fsz2;
    clim([Q(1) Q(end)]);
    set(cb,'Ticks',[500 5000 10000 15000]);

    set(s5,'box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);

    set(gcf,'renderer','painters');
    print(f5,'-dpng',[output '/FlowOnMars1_gQ'],'-r400');
    print(f5,'-dpdf',[output '/FlowOnMars1_gQ'],'-r400');
    end
end

%% plot x=g h_mars = h_earth
for a = 1
    if a == 1
    f6 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[10 5 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s6(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:Lh
        plot(grav,Q_h(:,i),'color',clr2(i,:),'linewidth',pl.line);
    end
    l(1) = ylabel('Discharge (Q) [m^3/s]');
    axis([grav(1) grav(end) 0 15000]);
    t(1) = text(gM-0.5,0,' Mars','color',clrM,'rotation',90);
    t(2) = text(gE-0.5,0,' Earth','color',clrE,'rotation',90);
    t(3) = text(xg,15000*y,' a');

    s6(2) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:Lh
        plot(grav,u(:,i),'color',clr2(i,:),'linewidth',pl.line);
    end
    l(2) = xlabel('Gravity (g) [m/s^2]');
    l(3) = ylabel('Velocity (u) [m/s]');
    axis([grav(1) grav(end) 0 7]);
    t(4) = text(xg,7*y,' e');

    s6(3) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:Lh
        plot([grav(1) grav(end)],[Rw(:,i) Rw(:,i)],'color',clr2(i,:),'linewidth',pl.line);
    end
    l(4) = ylabel('Hydraulic radius (R_w) [m]');
    axis([grav(1) grav(end) 0 15]);
    t(5) = text(xg,15*y,' b');

    s6(4) = axes('Position',[pl.lmarge+pl.inth+pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:Lh
        plot(grav,tau(:,i),'color',clr2(i,:),'linewidth',pl.line);
    end
    l(5) = xlabel('Gravity (g) [m/s^2]');
    l(6) = ylabel('Shear stress (\tau) [N/m^2]');
    axis([grav(1) grav(end) 0 100]);
    t(6) = text(xg,100*y,' f');

    s6(5) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:Lh
        plot(grav,f(:,i),'color',clr2(i,:),'linewidth',pl.line);
    end
    l(7) = ylabel({'Friction factor (f) [-]'});
    axis([grav(1) grav(end) 0 0.04]);
    t(7) = text(xg,0.04*y,' c');

    s6(6) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:Lh
        plot(grav,ust(:,i),'color',clr2(i,:),'linewidth',pl.line);
    end
    l(8) = xlabel('Gravity (g) [m/s^2]');
    l(9) = ylabel('Shear velocity (u_*) [m/s]');
    axis([grav(1) grav(end) 0 0.3]);
    t(8) = text(xg,0.3*y,' g');

    s6(7) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:Lh
        plot(grav,Fr(:,i),'color',clr2(i,:),'linewidth',pl.line);
    end
    l(10) = ylabel('Froude number (Fr) [-]');
    axis([grav(1) grav(end) 0 1.20]);
    plot([grav(1) grav(end)],[1 1],'--','color','k','linewidth',pl.line);
    t(9) = text(grav(end),0.93,'subcritical ','HorizontalAlignment','right');
    t(10) = text(grav(end),1.09,'supercritical ','HorizontalAlignment','right');
    t(11) = text(xg,1.2*y,' d');

    s6(8) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[1 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[1 100000000],'color',clrE,'linewidth',pl.line);
    for i = 1:Lh
        plot(grav,Re(:,i),'color',clr2(i,:),'linewidth',pl.line);
    end
    l(11) = xlabel('Gravity (g) [m/s^2]');
    l(12) = ylabel('Reynolds number (Re) [-]');
    axis([grav(1) grav(end) 150 10e7]);
    plot([grav(1) grav(end)],[500 500],'--','color','k','linewidth',pl.line);
    t(12) = text(grav(end),300,'laminar ','HorizontalAlignment','right');
    t(13) = text(grav(end),1100,'tubulent ','HorizontalAlignment','right');
    set(gca,'Yscale','log');
    t(14) = text(xg,10^(log10(150)+(8-log10(150))*y),' h');

    colormap(clr2)
    cb = colorbar('north','position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge+0.015 pl.wt 0.02],'Ticks',[h(1) h(end)]);
    cb.Label.String ='Water depth [m]';
    cb.Label.FontSize = pl.fsz2;
    clim([h(1) h(end)]);
    set(cb,'Ticks',[2 5 10 15]);

    set(s6,'box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);

    set(gcf,'renderer','painters');
    print(f6,'-dpng',[output '/FlowOnMars1_gh'],'-r400');
    print(f6,'-dpdf',[output '/FlowOnMars1_gh'],'-r400');
    end
end

%% plot x=g Q_mars = Q_earth Q=2000
for a = 1
    if a == 1
    f7 = figure('units','centimeters','PaperUnits','centimeters', ...
        'PaperPosition',[0 0 pl.width pl.height],'Position',[10 5 pl.width pl.height], ...
        'PaperSize',[pl.width pl.height],'visible','on');

    s7(1) = axes('Position',[pl.lmarge 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = Qselect
        plot(grav,h_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(1) = ylabel('Water depth (h) [m]');
    axis([grav(1) grav(end) 0 7]);
    t(1) = text(gM-0.5,0,' Mars','color',clrM,'rotation',90);
    t(2) = text(gE-0.5,0,' Earth','color',clrE,'rotation',90);
    t(3) = text(xg,7*y,' a');

    s7(2) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = Qselect
        plot(grav,u_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(2) = xlabel('Gravity (g) [m/s^2]');
    l(3) = ylabel('Velocity (u) [m/s]');
    axis([grav(1) grav(end) 0 4]);
    t(4) = text(xg,4*y,' e');

    s7(3) = axes('Position',[pl.lmarge+pl.inth+pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = Qselect
        plot(grav,Rw_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(4) = ylabel('Hydraulic radius (R_w) [m]');
    axis([grav(1) grav(end) 0 7]);
    t(5) = text(xg,7*y,' b');

    s7(4) = axes('Position',[pl.lmarge+pl.inth+pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = Qselect
        plot(grav,tau_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(5) = xlabel('Gravity (g) [m/s^2]');
    l(6) = ylabel('Shear stress (\tau) [N/m^2]');
    axis([grav(1) grav(end) 0 40]);
    t(6) = text(xg,40*y,' f');

    s7(5) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = Qselect
        plot(grav,f_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(7) = ylabel({'Friction factor (f) [-]'});
    axis([grav(1) grav(end) 0 0.03]);
    t(7) = text(xg,0.03*y,' c');

    s7(6) = axes('Position',[pl.lmarge+2*pl.inth+2*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = Qselect
        plot(grav,ust_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(8) = xlabel('Gravity (g) [m/s^2]');
    l(9) = ylabel('Shear velocity (u_*) [m/s]');
    axis([grav(1) grav(end) 0 0.2]);
    t(8) = text(xg,0.2*y,' g');

    s7(7) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    plot([gM gM],[0 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[0 100000000],'color',clrE,'linewidth',pl.line);
    for i = Qselect
        plot(grav,Fr_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(10) = ylabel('Froude number (Fr) [-]');
    axis([grav(1) grav(end) 0 1.20]);
    plot([grav(1) grav(end)],[1 1],'--','color','k','linewidth',pl.line);
    t(9) = text(grav(end),0.93,'subcritical ','HorizontalAlignment','right');
    t(10) = text(grav(end),1.09,'supercritical ','HorizontalAlignment','right');
    t(11) = text(xg,1.2*y,' d');

    s7(8) = axes('Position',[pl.lmarge+3*pl.inth+3*pl.wt pl.bmarge pl.wt pl.ht]);
    hold on
    plot([gM gM],[1 100000000],'color',clrM,'linewidth',pl.line);
    plot([gE gE],[1 100000000],'color',clrE,'linewidth',pl.line);
    for i = Qselect
        plot(grav,Re_Q(:,i),'color',clr3(i,:),'linewidth',pl.line);
    end
    l(11) = xlabel('Gravity (g) [m/s^2]');
    l(12) = ylabel('Reynolds number (Re) [-]');
    axis([grav(1) grav(end) 150 10e7]);
    plot([grav(1) grav(end)],[500 500],'--','color','k','linewidth',pl.line);
    t(12) = text(grav(end),300,'laminar ','HorizontalAlignment','right');
    t(13) = text(grav(end),1100,'tubulent ','HorizontalAlignment','right');
    set(gca,'Yscale','log');
    t(14) = text(xg,10^(log10(150)+(8-log10(150))*y),' h');

    colormap(clr3)
    cb = colorbar('north','position',[pl.lmarge+3*pl.inth+3*pl.wt 1-pl.tmarge+0.015 pl.wt 0.02],'Ticks',[Q(1) Q(end)]);
    cb.Label.String ='Discharge [m^3/s]';
    cb.Label.FontSize = pl.fsz2;
    clim([Q(1) Q(end)]);
    set(cb,'Ticks',[500 5000 10000 15000]);

    set(s7,'box','on','Layer','top', ...
        'XMinorTick','off','YMinorTick','off', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);

    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);

    set(gcf,'renderer','painters');
    print(f7,'-dpng',[output '/FlowOnMars1_gQ_2000m3s'],'-r400');
    print(f7,'-dpdf',[output '/FlowOnMars1_gQ_2000m3s'],'-r400');
    end
end

%%
close all