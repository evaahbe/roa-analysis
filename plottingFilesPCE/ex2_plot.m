function ex2_plot(filenames)

    load(filenames.filenamestoplot)
    contourmatrix = contourmatrix + xiep0;
    x = sdpvar(sys.xdim,1);
    xi = sdpvar((sys.p+1)*sys.xdim,1);

    nr_mu = length(sys.mu{:});
    mu = sdpvar(nr_mu,1);

    % PCE preprocessing
    vartot = [mu;x];

    varitot = [];
    for j = 1:nr_mu
        varitot = [varitot;sys.mu{j}.mu_coefs];
    end
    for j = 1:sys.xdim
        varitot = [varitot;xi((sys.p+1)*(j-1)+1:(sys.p+1)*j)'];
    end

    fxstoch = dynamics_iannelli(x,mu);
    fxi = clean(castPCEdynamics(vartot,varitot,sys,fxstoch)',numsetsRE.clean_thresh);

    T = 5;
    tinterval = linspace(0,T,500);

    Nx0 = 3; %number initial points
    Nm0 = 8; %number of values in mu range
    [xvec_mat,xvec_mat_spec] = WC_distribution(xi,Nx0,Nm0,tinterval,sys,contourmatrix,fxi);

    colors = getCustomColors();
    mycolors = colors.mycolors;


    figure(1)
    hold on   
    
    plot(xiep0(1), xiep0(2), 'ko','LineWidth',1.5)

    plot(contourmatrix(1,:),contourmatrix(2,:),'-r','LineWidth',2)

    lw = 1.7;

    plot(10,10,'--k','LineWidth',2) % quick legend workaround, to be ignored
    for j = 1:length(xvec_mat)-3
        for jj = 1:length(xvec_mat{j})
            hold on
            plot(xvec_mat{j}{jj}(:,2),xvec_mat{j}{jj}(:,3),'LineWidth',lw,'Color',mycolors(6,:))
        end
        plot(xvec_mat{j}{jj}(1,2),xvec_mat{j}{jj}(1,3),'x','LineWidth',3,'Color',mycolors(6,:))
    end
    plot(xvec_mat{j+1}{1}(1,2),xvec_mat{j+1}{1}(1,3),'X','LineWidth',3,'Color',mycolors(6,:))
    plot(xvec_mat{j+1}{1}(:,2),xvec_mat{j+1}{1}(:,3),'LineWidth',lw,'Color',mycolors(6,:))
    plot(xvec_mat{j+2}{1}(1,2),xvec_mat{j+2}{1}(1,3),'X','LineWidth',3,'Color',mycolors(6,:))
    plot(xvec_mat{j+2}{1}(:,2),xvec_mat{j+2}{1}(:,3),'LineWidth',lw,'Color',mycolors(6,:))
    plot(xvec_mat{j+3}{1}(1,2),xvec_mat{j+3}{1}(1,3),'x','LineWidth',3,'Color',mycolors(6,:))
    plot(xvec_mat{j+3}{1}(:,2),xvec_mat{j+3}{1}(:,3),'LineWidth',lw,'Color',mycolors(6,:))


    for j = 1:length(xvec_mat)-3
        plot(xvec_mat_spec{j}(:,2),xvec_mat_spec{j}(:,3),'--k','LineWidth',2.2)
    end
    lw = 1.2;
    set(gca,...
    'LineWidth', 1.5,...
    'Units','normalized',...                  
    'FontUnits','points',...
    'FontSize',20)
    xlabel('$\bar{x}_{1_0}$','FontSize',35,'interpreter','latex');
    ylabel('$\bar{x}_{2_0}$','FontSize',35,'interpreter','latex');  
    xlim([-1. 2]);
    ylim([-3 5]);
    legend(  '$\mathcal{R}_0$','$\bar{x}_0(t)$', '$x(t,c \in [0.9,1.1])$', 'LineWidth',lw,'FontSize',20,'NumColumns',1,'interpreter','latex')


end



function [xvec_mat,xvec_mat_spec] = WC_distribution(xi,Nx0,Nm0,tinterval,sys,contourmatrix,fxi)
        

    inioff = 1/Nx0;

    paravec = linspace(0,1,length(contourmatrix(1,:)));
    evalvec = linspace(0,1-inioff,Nx0);

    mean_spline = spline(paravec,contourmatrix);

    for i = 1:length(evalvec)
        roa_mean(:,i) = fnval(mean_spline,evalvec(i));
    end

    %add unstable initial points for illustration
    roa_mean = [roa_mean,[-0.7434;-1.247]];
    roa_mean = [roa_mean,[-0.793;-1.15]];
    roa_mean = [roa_mean,[-0.676;-1.35]];


    for j = 1:length(roa_mean)

        xini = kron(roa_mean(:,j),[1;zeros(sys.p,1)]);

        %compute mean trajectory from PCE
        [fxi_sym_t,xst] = comp_fxi_sym_t(sys,xi,fxi);
        fxi_sym_fun = odeFunction(fxi_sym_t,xst);
        [tvec_PCE, xvec_PCE] = ode45(fxi_sym_fun, tinterval,xini);

        xvec0_PCE = [xvec_PCE(:,1),xvec_PCE(:,sys.p+2)];
        xvec_mat_spec{j} = [tvec_PCE,xvec0_PCE];

        %compute Monte Carlo trajectories
        mu_array = linspace(sys.mu{1}.mu_A,sys.mu{1}.mu_B,Nm0);            
        for h = 1:length(mu_array)
            [tvec_minus, xvec_minus] = ode45(@(t,xx) dynamics_iannelli_INT(t,xx,mu_array(h)),tinterval,[xini(1);xini(sys.p+2)]);
            xvec_mat{j}{h} = [tvec_minus,xvec_minus];
        end

    end
                  
end




function [fxi_sym_t,xst] = comp_fxi_sym_t(sys,xi,fxi)
    syms xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10 xs11 xs12 xs13 xs14 xs15 xs16 xs17 xs18 xs19 xs20 xs21 xs22 xs23 xs24 xs25 xs26 xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t) xst15(t) xst16(t) xst17(t) xst18(t) xst19(t) xst20(t) xst21(t) xst22(t) xst23(t) xst24(t) xst25(t) xst26(t)
    if sys.p == 0
        xs = [xs1;xs2];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2],[xst1(t) xst2(t)]);
        xst = [xst1(t) xst2(t)];
    elseif sys.p == 1
        xs = [xs1;xs2;xs3;xs4];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4],[xst1(t) xst2(t) xst3(t) xst4(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t)];
    elseif sys.p == 2
        xs = [xs1;xs2;xs3;xs4;xs5;xs6];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t)];
    elseif sys.p == 3
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7; xs8];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t)];
    elseif sys.p == 4
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7; xs8; xs9; xs10];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t)];
    end
end
