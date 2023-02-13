%% Flow on Mars 6 - Sediment mixture flux
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
    Q = 2000; %500:500:15000; %[m3/s] Discharge
    %Qselect = 4;
    W = 200; %[m] Channel width
    grav = 1:0.25:12; %[m/s2] Gravitational acceleration
    gselect = [12 36];
    S = 0.001; %logspace(-4,-2,27);
    %Sselect = 14; %[m/m] Slope
    TC = 4; %[degrees C] Temperature
    rho = 1000; %[kg/m3] Water density
    
    rhos = 2900;                    %[kg/m3] Sediment density
    load('FlowOnMars5_GrainSizeDistribution.mat');
    Drough = D90;                   %[m] Grain size, only used for roughness
    clear D10 D50 D90;
    D50 = D; clear D;
    %Ha = 1e-20; %[J | Nm] Hamaker constant, order of e-20
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
    %kk = ((1-(1-pi/(3*2^0.5)))/(1-n))^(1/3)-1; %[-] ?
    ks = 3*Drough; %[m] Nikuradse
    clear TC TK mu CDryB Drough
end

%% allocate
for a=1
    LQ = length(Q);
    Lg = length(grav);
    LD = length(D50);
    LS = length(S);
    
    h = NaN(Lg,LQ);
    Rw = NaN(Lg,LQ);
    C = NaN(Lg,LQ);
    u = NaN(Lg,LQ,LS);
    Fr = NaN(Lg,LQ,LS);
    tau = NaN(Lg,LQ,LS);
    ust = NaN(Lg,LQ,LS);
    
    ws = NaN(Lg,LD);
    Dst = NaN(Lg,LD);
    shield = NaN(Lg,LD,LQ,LS);
    
    shield_cr = NaN(Lg,LD,LQ,LS);
    S0 = NaN(Lg,LD);
    
    Qb_E50 = NaN(Lg,LD,LQ,LS);
    qb_E50 = NaN(Lg,LD,LQ,LS);
    Qb_FB = NaN(Lg,LD,LQ,LS);
    qb_FB = NaN(Lg,LD,LQ,LS);
    p_EF = NaN(Lg,LD,LQ,LS);
    Qb_EF = NaN(Lg,LD,LQ,LS);
    qb_EF = NaN(Lg,LD,LQ,LS);
    Qb_vR = NaN(Lg,LD,LQ,LS);
    qb_vR = NaN(Lg,LD,LQ,LS);
    qb_EH = NaN(Lg,LD,LQ,LS);
    Qb_EH = NaN(Lg,LD,LQ,LS);
    
    a_E50 = NaN(1,LD);
    a_EF = NaN(1,LD);
    a_dL = NaN(Lg,LQ,LS);
    a_vR = NaN(Lg,LQ,LS);
    
    Rouse = NaN(Lg,LD,LQ,LS);
    
    Ca_E50 = NaN(Lg,LD,LQ,LS);
    labda_EF = NaN(Lg,LD,LQ,LS);
    Ca_EF = NaN(Lg,LD,LQ,LS);
    Ca_vR = NaN(Lg,LD,LQ,LS);
    Ca_dL = NaN(Lg,LD,LQ,LS);
    
    qs_E50 = NaN(Lg,LD,LQ,LS);
    qs_EF = NaN(Lg,LD,LQ,LS);
    qs_dL = NaN(Lg,LD,LQ,LS);
    qs_vR = NaN(Lg,LD,LQ,LS);
    
    qb_E50_frac = NaN(Lg,LD,LQ,LS);
    qb_FB_frac = NaN(Lg,LD,LQ,LS);
    qb_EF_frac = NaN(Lg,LD,LQ,LS);
    qb_vR_frac = NaN(Lg,LD,LQ,LS);
    qs_E50_frac = NaN(Lg,LD,LQ,LS);
    qs_dL_frac = NaN(Lg,LD,LQ,LS);
    qs_EF_frac = NaN(Lg,LD,LQ,LS);
    qs_vR_frac = NaN(Lg,LD,LQ,LS);
    qt_EH_frac = NaN(Lg,LD,LQ,LS);
    
    qt_vR_frac = NaN(Lg,LD,LQ,LS);
    qt_E50_frac = NaN(Lg,LD,LQ,LS);
    qt_EF_frac = NaN(Lg,LD,LQ,LS);
    qt_dL_frac = NaN(Lg,LD,LQ,LS);
    
    qt_vR_sum = NaN(Lg,LQ,LS);
    qt_E50_sum = NaN(Lg,LQ,LS);
    qt_EF_sum = NaN(Lg,LQ,LS);
    qt_dL_sum = NaN(Lg,LQ,LS);
    qt_EH_sum = NaN(Lg,LQ,LS);

end

%% Hydro parameters Q_mars = Q_earth
h_guess = 3;
for j = 1:LS
    for i = 1:LQ
        for g = 1:Lg
            h(g,i,j) = h_guess+1;
            while abs(h_guess-h(g,i,j))>0.0001
                h_guess = h(g,i,j);
                Rw(g,i) = (h_guess*W)/(h_guess+h_guess+W); %[m] Hydraulic radius
                C(g,i) = 5.75*grav(g)^0.5*log10(12*h_guess/(ks)); %[m0.5/s] Chezy roughness
                u(g,i,j) = C(g,i)*(Rw(g,i)*S(j))^0.5; %[m/s] Velocity
                h(g,i,j) = Q(i)/(W*u(g,i,j)); %[m] Water depth
            end
            Fr(g,i,j) = u(g,i,j)/(grav(g)*h(g,i,j))^0.5; %[-] Froude number
            %Re(g,i,j) = u(g,i,j)*h(g,i,j)/v; %[-] Reynolds number
            tau(g,i,j) = rho*grav(g)*Rw(g,i)*S(j); %[N/m2] Bed shear stress
            ust(g,i,j) = (tau(g,i,j)/rho)^0.5; %[m/s] Shear velocity
            %lamlyr_Q(g) = 11.63*v/ust_Q(g); %[m] Laminar sublayer thickness
        end
    end
end
clear i g j h_guess Rw Re

%% Sediment parameters
for d = 1:LD
    for g = 1:Lg
        ws(g,d) = R*grav(g)*D50(d)^2/(C1*v+(0.75*C2*R*grav(g)*D50(d)^3)^0.5); %[m/s] Settling velocity (Ferguson and Church 2004)
        Dst(g,d) = D50(d)*(R*grav(g)/v^2)^(1/3); %[-] Dimensionless particle parameter/Bonnefille number
        
        % Q_mars = Q_earth
        for i = 1:LQ
            for j = 1:LS
                shield(g,d,i,j) = tau(g,i,j)/((rhos-rho)*grav(g)*D50(d)); %[-] Shields parameter/Particle mobility parameter
            end
        end
    end
end
clear C1 C2 d g i

%% Critical threshold of motion
acc = 1e-5;
for d = 1:LD
    for g = 1:Lg
        for i = 1:LQ
            for j = 1:LS
                %Zanke 2003 iterated
                shieldini = shield(g,d,i,j);
                minimise = 1;
                while max(minimise)>acc
                    tau_ini = shieldini *((rhos-rho)*grav(g)*D50(d));
                    ust_ini = (tau_ini/rho)^0.5;
                    Rest_ini = D50(d)*ust_ini/v;
                    uprmsb_ust = 0.31*Rest_ini*exp(-0.1*Rest_ini)+1.8*exp(-0.88*D50(d)/(shieldini*R*D50(d)/sin(S(j))))*(1-exp(-0.1*Rest_ini));
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
                shield_cr(g,d,i,j) = shieldini;
            end
        end
    end
end
clear acc Rest_ini ust_ini tau_ini shieldini Pt B uy_ust ub_ust uprmsb_ub uprmsb_ust K shield1 minimise phi n d g

%% get 1 value for shield_cr
% for i = 1:LQ
%     for g = 1:Lg
%         for j = 1:LS
%             shield_cr_cal = InterX([D50; shield_cr(g,:,i,j)],[D50; shield(g,:,i,j)]);
%             shield_cr_val(g,i,j) = shield_cr_cal(2);
%         end
%     end
% end
% clear shield_cr
% shield_cr = shield_cr_val;
% clear shield_cr_calc shield_cr_val g

%% Bedload transport
for i = 1:LQ
    for g = 1:Lg
        for d = 1:LD
            for j = 1:LS
                S0(g,d,i,j) = max(((shield(g,d,i,j)-shield_cr(g,d,i,j))/shield_cr(g,d,i,j)),0);
                
                %A(theta-theta_cr)^B
                Qb_E50(g,d,i,j) = 3.97*(max(shield(g,d,i,j)-shield_cr(g,d,i,j),0))^1.5; %[-] Einstein 1950 as in de Leeuw 2020 (I have the paper, but am unable to rewrite)
                qb_E50(g,d,i,j) = Qb_E50(g,d,i,j)*(grav(g)*R*D50(d)^3)^0.5; %[m3/ms]
%                 Qb_FB(g,d,i,j) = 5.7*(max(shield(g,d,i,j)-shield_cr(g,d,i,j),0))^1.5; %[-] Fernandez Luque and van Beek 1976 as in de Leeuw 2020 (payed download)
%                 qb_FB(g,d,i,j) = Qb_FB(g,d,i,j)*(grav(g)*R*D50(d)^3)^0.5; %[m3/ms]
%                 
%                 %other
%                 p_EF(g,d,i,j) = min((1+((pi/6)*1/(max(shield(g,d,i,j)-shield_cr(g,d,i,j),0)))^4)^(-0.25),1); %[-] probability | beta=1 Engelund and Fredsoe 1982 as in Garcia and Parker 1991
%                 Qb_EF(g,d,i,j) = 5*p_EF(g,d,i,j)*(max(shield(g,d,i,j)^0.5-0.7*shield_cr(g,d,i,j)^0.5,0)); %[-] Engelund and Fredsoe 1976 (9.3*(pi/6)=~5) | theta_cr=0.05,0.06
%                 qb_EF(g,d,i,j) = Qb_EF(g,d,i,j)*(grav(g)*R*D50(d)^3)^0.5; %[m3/ms]
%                 Qb_vR(g,d,i,j) = 0.053*Dst(g,d)^-0.3*S0(g,d,i,j)^2.1; %[-] van Rijn 1984, 200<D<2000mum
%                 qb_vR(g,d,i,j) = Qb_vR(g,d,i,j)*(grav(g)*R*D50(d)^3)^0.5; %[m3/ms]
%                 
%                 %no theta_cr
%                 qb_EH(g,d,i,j) = 0.05*u(g,i,j)^5/(grav(g)^0.5*C(g,i)^3*R^2*D50(d)); %[m3/ms]
%                 Qb_EH(g,d,i,j) = qb_EH(g,d,i,j)/(grav(g)*R*D50(d)^3)^0.5; %[-] Engelund and Hansen as in D3D
            end
        end
    end
end

%% Suspended transport
tic
for i = 1:LQ
    for g=1:Lg
        for d=1:LD
            for j = 1:LS
                %reference height
                a_E50(d) = 2*D50(d);
%                 a_EF(d) = 2*D50(d);
%                 a_dL(g,i,j) = 0.1*h(g,i,j);
%                 a_vR(g,i,j) = max(0.01*h(g,i,j),ks);
                
                %Rouse number
                Rouse(g,d,i,j) = ws(g,d)/(Kappa*ust(g,i,j));
                
                %Ca / Es
                Ca_E50(g,d,i,j) = (1/32.2)*Qb_E50(g,d,i,j)/shield(g,d,i,j)^0.5;
%                 labda_EF(g,d,i,j) = (max((shield(g,d,i,j)-shield_cr(g,d,i,j)-(1*p_EF(g,d,i,j)*pi/6))/(0.027*(R+1)*shield(g,d,i,j)),0))^0.5;
%                 Ca_EF(g,d,i,j) = 0.65/(1+labda_EF(g,d,i,j)^-1)^3;
%                 Ca_vR(g,d,i,j) = 0.015*(D50(d)/a_vR(g,i,j))*S0(g,d,i,j)^1.5/Dst(g,d)^0.3;
%                 Ca_dL(g,d,i,j) = 0.000474*(ust(g,i,j)/ws(g,d))^1.77*Fr(g,i,j)^1.18;
                
                %Integral
                fun = @(z) Ca_E50(g,d,i,j) * ( ((h(g,i,j)-z)./z) * (a_E50(d)/(h(g,i,j)-a_E50(d))) ).^Rouse(g,d,i,j);
                qs_E50(g,d,i,j) = integral(fun,a_E50(d),h(g,i,j)); %Einstein 1950
%                 fun = @(z) Ca_EF(g,d,i,j) * ( ((h(g,i,j)-z)./z) * (a_EF(d)/(h(g,i,j)-a_EF(d))) ).^Rouse(g,d,i,j);
%                 qs_EF(g,d,i,j) = integral(fun,a_EF(d),h(g,i,j)); %Engelund and Fredsoe 1976
%                 fun = @(z) Ca_dL(g,d,i,j) * ( ((h(g,i,j)-z)./z) * (a_dL(g,i,j)/(h(g,i,j)-a_dL(g,i,j))) ).^Rouse(g,d,i,j);
%                 qs_dL(g,d,i,j) = integral(fun,a_dL(g,i,j),h(g,i,j)); %de Leeuw 2020
%                 fun = @(z) Ca_vR(g,d,i,j) * ( ((h(g,i,j)-z)./z) * (a_vR(g,i,j)/(h(g,i,j)-a_vR(g,i,j))) ).^Rouse(g,d,i,j);
%                 qs_vR(g,d,i,j) = integral(fun,a_vR(g,i,j),h(g,i,j)); %van Rijn 1984
            end
        end
    end
end
toc
clear g d i Kappa Qb_E50 Qb_FB Qb_EF Qb_vR Qb_EH p_EF labda_EF S0 Fr
clear a_E50 a_EF a_dL a_vR Rouse Ca_E50 Ca_EF Ca_vR Ca_dL fun

%% multiply sediment classes with bed availability from lognormal sediment curve
for a = 1
    for i = 1:LQ
        for g = 1:Lg
            for j = 1:LS
                qb_E50_frac(g,:,i,j) = qb_E50(g,:,i,j).*bins_y_perc/100;
                %qb_FB_frac(g,:,i,j) = qb_FB(g,:,i,j).*bins_y_perc/100;
                %qb_EF_frac(g,:,i,j) = qb_EF(g,:,i,j).*bins_y_perc/100;
                %qb_vR_frac(g,:,i,j) = qb_vR(g,:,i,j).*bins_y_perc/100;
                qs_E50_frac(g,:,i,j) = qs_E50(g,:,i,j).*bins_y_perc/100;
                %qs_dL_frac(g,:,i,j) = qs_dL(g,:,i,j).*bins_y_perc/100;
                %qs_EF_frac(g,:,i,j) = qs_EF(g,:,i,j).*bins_y_perc/100;
                %qs_vR_frac(g,:,i,j) = qs_vR(g,:,i,j).*bins_y_perc/100;
                %qt_EH_frac(g,:,i,j) = qb_EH(g,:,i,j).*bins_y_perc/100;
            end
        end
    end
    
%     qt_vR_frac = qb_vR_frac + qs_vR_frac;
    qt_E50_frac = qb_E50_frac + qs_E50_frac;
%     qt_EF_frac = qb_EF_frac + qs_EF_frac;
%     qt_dL_frac = qb_FB_frac + qs_dL_frac;
%     
%     qt_vR_sum = sum(qt_vR_frac,2);
    qt_E50_sum = sum(qt_E50_frac,2);
%     qt_EF_sum = sum(qt_EF_frac,2);
%     qt_dL_sum = sum(qt_dL_frac,2);
%     qt_EH_sum = sum(qt_EH_frac,2);
    
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
    CMAPOBJ = clrmap('read','Discharge.clrmap');
    %clr2 = clrmap(CMAPOBJ,Lh);
    clr3 = clrmap(CMAPOBJ,LQ);
    CMAPOBJ = clrmap('read','GrainSize.clrmap');
    clr4 = clrmap(CMAPOBJ,LD);
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
    hold on
    plot([1.351 1.352],[10^-1 10^4],':','color',[0 0.348 0.15],'linewidth',pl.line)
    plot([3.75 3.75],[10^-1 10^4],':','color',clr1(1,:),'linewidth',pl.line)
    plot([8.87 8.87],[10^-1 10^4],':','color',[0.298 0.2 0],'linewidth',pl.line)
    plot([9.81 9.81],[10^-1 10^4],':','color',clr1(end,:),'linewidth',pl.line)
    plot(grav,qt_E50_sum(:,:)*W,'-','color',[0 0 0],'linewidth',pl.line)
%     plot(grav,qt_EF_sum(:,:,Qselect,Sselect),'--','color',[0 0 0],'linewidth',pl.line)
%     plot(grav,qt_dL_sum(:,:,Qselect,Sselect),'-.','color',[0 0 0],'linewidth',pl.line)
%     plot(grav,qt_vR_sum(:,:,Qselect,Sselect),':','color',[0 0 0],'linewidth',pl.line)
%     plot(grav,qt_EH_sum(:,:,Qselect,Sselect),'-','color',[0 0 0 0.2],'linewidth',pl.line)
    set(gca,'Yscale','log');
    xlim([grav(1) grav(end)]);
    ylim([10^1 10^4]);
    l(1) = xlabel('Gravity acceleration (g) [m/s^2]');
    l(2) = ylabel('Total transport [m^3/s]');
    t(1) = text(xg,10^(1+(4-1)*y),' a');
    t(2) = text(1.351+0.25,500,'Titan','color',[0 0.348 0.15],'HorizontalAlignment', 'right');
    t(3) = text(3.75+0.25,280,'Mars','color',[clr1(1,:)],'HorizontalAlignment', 'right');
    t(4) = text(8.87+0.25,155,'Venus','color',[0.298 0.2 0],'HorizontalAlignment', 'right');
    t(5) = text(9.81+0.25,140,'Earth','color',[clr1(end,:)],'HorizontalAlignment', 'right');
    set(t(2:5),'Rotation',90);
    t(6) = text(grav(end),1.5*10^1,'W = 200 m + S = 0.001 m/m + Q = 2000 m^3/s   ','HorizontalAlignment', 'right');

    s1(2) = axes('Position',[pl.lmarge+pl.wt+pl.inth 1-pl.tmarge-pl.ht pl.wt pl.ht]);
    hold on
    g = gselect(1);
    plot(D50,qt_E50_frac(g,:)/qt_E50_sum(g,:)*100,'-','color',clr1(1,:),'linewidth',pl.line);
    g = gselect(2);
    plot(D50,qt_E50_frac(g,:)/qt_E50_sum(g,:)*100,'-','color',clr1(end,:),'linewidth',pl.line);
%     g = gselect(1);
%     plot(D50,qt_EF_frac(g,:,Qselect,Sselect)/qt_EF_sum(g,:,Qselect,Sselect)*100,'--','color',clr1(1,:),'linewidth',pl.line);
%     plot(D50,qt_dL_frac(g,:,Qselect,Sselect)/qt_dL_sum(g,:,Qselect,Sselect)*100,'-.','color',clr1(1,:),'linewidth',pl.line);
%     plot(D50,qt_vR_frac(g,:,Qselect,Sselect)/qt_vR_sum(g,:,Qselect,Sselect)*100,':','color',clr1(1,:),'linewidth',pl.line);
%     plot(D50,qt_EH_frac(g,:,Qselect,Sselect)/qt_EH_sum(g,:,Qselect,Sselect)*100,'-','color',[clr1(1,:) 0.2],'linewidth',pl.line);
%     g = gselect(2);
%     plot(D50,qt_EF_frac(g,:,Qselect,Sselect)/qt_EF_sum(g,:,Qselect,Sselect)*100,'--','color',clr1(end,:),'linewidth',pl.line);
%     plot(D50,qt_dL_frac(g,:,Qselect,Sselect)/qt_dL_sum(g,:,Qselect,Sselect)*100,'-.','color',clr1(end,:),'linewidth',pl.line);
%     plot(D50,qt_vR_frac(g,:,Qselect,Sselect)/qt_vR_sum(g,:,Qselect,Sselect)*100,':','color',clr1(end,:),'linewidth',pl.line);
%     plot(D50,qt_EH_frac(g,:,Qselect,Sselect)/qt_EH_sum(g,:,Qselect,Sselect)*100,'-','color',[clr1(end,:) 0.2],'linewidth',pl.line);
    %set(gca,'Yscale','log');
    set(gca,'Xscale','log');
    xlim([D50(1) D50(end)]);
    ylim([0 50]);
    l(3) = xlabel('Sediment class');
    l(4) = ylabel('% of total');
    t(7) = text(xd,50*y,' b');
    set(gca,'XTick',D50,'XtickLabel',bins_txt);
    
%     s1(3) = axes('Position',[pl.lmarge pl.bmarge pl.wt pl.ht]);
%     hold on
%     plot(Q,squeeze(qt_E50_sum(gselect(1),:,:,Sselect)),'-','color',clr1(1,:),'linewidth',pl.line)
% %     plot(Q,squeeze(qt_EF_sum(gselect(1),:,:,Sselect)),'--','color',clr1(1,:),'linewidth',pl.line)
% %     plot(Q,squeeze(qt_dL_sum(gselect(1),:,:,Sselect)),'-.','color',clr1(1,:),'linewidth',pl.line)
% %     plot(Q,squeeze(qt_vR_sum(gselect(1),:,:,Sselect)),':','color',clr1(1,:),'linewidth',pl.line)
% %     plot(Q,squeeze(qt_EH_sum(gselect(1),:,:,Sselect)),'-','color',[clr1(1,:) 0.2],'linewidth',pl.line)
%     plot(Q,squeeze(qt_E50_sum(gselect(2),:,:,Sselect)),'-','color',clr1(end,:),'linewidth',pl.line)
% %     plot(Q,squeeze(qt_EF_sum(gselect(2),:,:,Sselect)),'--','color',clr1(end,:),'linewidth',pl.line)
% %     plot(Q,squeeze(qt_dL_sum(gselect(2),:,:,Sselect)),'-.','color',clr1(end,:),'linewidth',pl.line)
% %     plot(Q,squeeze(qt_vR_sum(gselect(2),:,:,Sselect)),':','color',clr1(end,:),'linewidth',pl.line)
% %     plot(Q,squeeze(qt_EH_sum(gselect(2),:,:,Sselect)),'-','color',[clr1(end,:) 0.2],'linewidth',pl.line)
%     set(gca,'Yscale','log');
%     %set(gca,'Xscale','log');
%     xlim([Q(1) Q(end)]);
%     ylim([10^-1 10^2]);
%     l(5) = xlabel('Discharge [m^3/s]');
%     l(6) = ylabel('Total transport [m^3/ms]');
%     t(3) = text(xq,10^(-1+(2--1)*y),' c');
%     t(9) = text(Q(end),1.5*10^-1,'S = 0.001 m/m   ','HorizontalAlignment', 'right');
%     
%     s1(4) = axes('Position',[pl.lmarge+pl.wt+pl.inth pl.bmarge pl.wt pl.ht]);
%     hold on
%     plot(S,squeeze(qt_E50_sum(gselect(1),:,Qselect,:)),'-','color',clr1(1,:),'linewidth',pl.line)
% %     plot(S,squeeze(qt_EF_sum(gselect(1),:,Qselect,:)),'--','color',clr1(1,:),'linewidth',pl.line)
% %     plot(S,squeeze(qt_dL_sum(gselect(1),:,Qselect,:)),'-.','color',clr1(1,:),'linewidth',pl.line)
% %     plot(S,squeeze(qt_vR_sum(gselect(1),:,Qselect,:)),':','color',clr1(1,:),'linewidth',pl.line)
% %     plot(S,squeeze(qt_EH_sum(gselect(1),:,Qselect,:)),'-','color',[clr1(1,:) 0.2],'linewidth',pl.line)
%     plot(S,squeeze(qt_E50_sum(gselect(2),:,Qselect,:)),'-','color',clr1(end,:),'linewidth',pl.line)
% %     plot(S,squeeze(qt_EF_sum(gselect(2),:,Qselect,:)),'--','color',clr1(end,:),'linewidth',pl.line)
% %     plot(S,squeeze(qt_dL_sum(gselect(2),:,Qselect,:)),'-.','color',clr1(end,:),'linewidth',pl.line)
% %     plot(S,squeeze(qt_vR_sum(gselect(2),:,Qselect,:)),':','color',clr1(end,:),'linewidth',pl.line)
% %     plot(S,squeeze(qt_EH_sum(gselect(2),:,Qselect,:)),'-','color',[clr1(end,:) 0.2],'linewidth',pl.line)
%     set(gca,'Yscale','log');
%     set(gca,'Xscale','log');
%     xlim([S(1) S(end)]);
%     ylim([10^-1 10^2]);
%     l(7) = xlabel('Slope [m/m]');
%     l(8) = ylabel('Total transport [m^3/ms]');
%     t(4) = text(xq,10^(-1+(2--1)*y),' d');
%     t(10) = text(S(end),1.5*10^-1,'Q = 2000 m^3/s   ','HorizontalAlignment', 'right');

    %leg1 = legend(s1(1),'Einstein 1950','Engelund and Fredsoe 1976','de Leeuw 2020','van Rijn 1984','Engelund and Hansen 1967','location','northeast');
    leg2 = legend(s1(2),'Mars','Earth','location','northeast');
    %leg2 = legend(s1(2),'Einstein 1950','Engelund and Fredsoe 1976','de Leeuw 2020','van Rijn 1984','Engelund and Hansen 1967','location','northeast');
    %leg3 = legend(s1(3),'Einstein 1950','Engelund and Fredsoe 1976','de Leeuw 2020','van Rijn 1984','Engelund and Hansen 1967','location','northeast');
    %leg4 = legend(s1(4),'Einstein 1950','Engelund and Fredsoe 1976','de Leeuw 2020','van Rijn 1984','Engelund and Hansen 1967','location','northeast');
    
    xtickangle(s1(2),90);
    set(s1,'box','on','Layer','top', ...
        'XMinorTick','on','YMinorTick','on', ...
        'FontSize',pl.fsz,'LineWidth',pl.line_ax,'TickDir','in','YDir','normal', ...
        'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0], ...
        'MinorGridColor',[0 0 0],'GridColor',[0.5 0.5 0.5]);
    
    set(t,'FontSize',pl.fsz2);
    set(l,'FontSize',pl.fsz2);
    
    set(gcf,'renderer','painters');
    print(f1,'-dpng',[output '/FlowOnMars6_SedimentMixtureFlux_QS'],'-r400');
    print(f1,'-dpdf',[output '/FlowOnMars6_SedimentMixtureFlux_QS'],'-r400');
    
end
