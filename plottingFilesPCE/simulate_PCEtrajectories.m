clear all
clear yalmip

sys          = struct();
fns          = struct();
clean_thresh = 1e-6;

sys.p = 3;
[sys,fns] = initializeIannelli(sys,fns); %choose the initialization file of the desired system 
xini = kron([1;1],[1;zeros(sys.p,1)]);
integration_time = 10;



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

fxstoch = fns.system_dyn_poly(x,mu);
fxi = clean(castPCEdynamics(vartot,varitot,sys,fxstoch)',clean_thresh);

[fxi_sym_fun] = comp_fxi_sym_t(sys,xi,fxi);

tinterval = [0 integration_time];
[tvec, xvec] = ode45(fxi_sym_fun, tinterval,xini);


colors = getCustomColors();

figure
hold on
title('PCE trajectories','FontSize',20)
set(gca,...
    'LineWidth', 1.5,...
    'Units','normalized',...                  
    'FontUnits','points',...
    'FontSize',20)

for i = 1:(sys.xdim*(sys.p+1))
    
    xdim_nr = ceil(i/(sys.p+1));
    p_nr = mod((i-1),(sys.p+1)); 
    disp(p_nr)
    plot(tvec,xvec(:,i),'LineWidth',1.5,'Color',colors.mycolors(i,:))
    
    legendentries{i} =  strcat('$\bar{x}_{',num2str(xdim_nr),'_',num2str(p_nr),'}$');
    
end
legend(legendentries,'LineWidth',1.2,'FontSize',25,'NumColumns',2,'interpreter','latex')
xlabel('$t$','FontSize',25,'interpreter','latex');
ylabel('$\bar{x}(t)$','FontSize',25,'interpreter','latex');
hold off


    
function [fxi_sym_fun] = comp_fxi_sym_t(sys,xi,fxi)

    syms xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10 xs11 xs12 xs13 xs14 xs15 xs16 xs17 xs18 xs19 xs20 xs21 xs22 xs23 xs24 xs25 xs26 xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t) xst15(t) xst16(t) xst17(t) xst18(t) xst19(t) xst20(t) xst21(t) xst22(t) xst23(t) xst24(t) xst25(t) xst26(t)
    if ((sys.p+1)*sys.xdim) == 1
        xs = [xs1];
        fxi_sym = makesdpvarsym(fxi,xi,xs); 
        fxi_sym_t = subs(fxi_sym,[xs1],[xst1(t)]);
        xst = [xst1(t)];
    elseif ((sys.p+1)*sys.xdim) == 2
        xs = [xs1;xs2];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2],[xst1(t) xst2(t)]);
        xst = [xst1(t) xst2(t)];
    elseif ((sys.p+1)*sys.xdim) == 3
        xs = [xs1;xs2;xs3];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3],[xst1(t) xst2(t) xst3(t)]);
        xst = [xst1(t) xst2(t) xst3(t)];
    elseif ((sys.p+1)*sys.xdim) == 4
        xs = [xs1;xs2;xs3;xs4];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4],[xst1(t) xst2(t) xst3(t) xst4(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t)];
    elseif ((sys.p+1)*sys.xdim) == 5
        xs = [xs1;xs2;xs3;xs4;xs5];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t)];
    elseif ((sys.p+1)*sys.xdim) == 6
        xs = [xs1;xs2;xs3;xs4;xs5;xs6];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t)];
    elseif ((sys.p+1)*sys.xdim) == 7
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t)];
    elseif ((sys.p+1)*sys.xdim) == 8
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7;xs8];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t)];
    elseif ((sys.p+1)*sys.xdim) == 9
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7;xs8;xs9];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t)];
    elseif ((sys.p+1)*sys.xdim) == 10
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7;xs8;xs9;xs10];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t)];
    elseif ((sys.p+1)*sys.xdim) == 11
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7; xs8; xs9; xs10; xs11];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10 xs11],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t)];
    elseif ((sys.p+1)*sys.xdim) == 12
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7; xs8; xs9; xs10; xs11; xs12];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10 xs11 xs12],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t)];
    elseif ((sys.p+1)*sys.xdim) == 13
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7; xs8; xs9; xs10; xs11; xs12; xs13];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10 xs11 xs12 xs13],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t)];
    elseif ((sys.p+1)*sys.xdim) == 14
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7; xs8; xs9; xs10; xs11; xs12; xs13; xs14];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10 xs11 xs12 xs13 xs14],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t)];
    elseif ((sys.p+1)*sys.xdim) == 15
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7; xs8; xs9; xs10; xs11; xs12; xs13; xs14;xs15];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10 xs11 xs12 xs13 xs14 xs15],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t) xst15(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t) xst15(t)];
    elseif ((sys.p+1)*sys.xdim) == 16
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7; xs8; xs9; xs10; xs11; xs12; xs13; xs14;xs15;xs16];
        fxi_sym = comp_dynamics(x,mu,sys,xim,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10 xs11 xs12 xs13 xs14 xs15 xs16],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t) xst15(t) xst16(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t) xst15(t) xst16(t) ];
    elseif ((sys.p+1)*sys.xdim) == 17
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7; xs8; xs9; xs10; xs11; xs12; xs13; xs14;xs15;xs16;xs17];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10 xs11 xs12 xs13 xs14 xs15 xs16 xs17],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t) xst15(t) xst16(t) xst17(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t) xst15(t) xst16(t) xst17(t)];
    elseif ((sys.p+1)*sys.xdim) == 18
        xs = [xs1;xs2;xs3;xs4;xs5;xs6;xs7; xs8; xs9; xs10; xs11; xs12; xs13; xs14;xs15;xs16;xs17;xs18];
        fxi_sym = makesdpvarsym(fxi,xi,xs);
        fxi_sym_t = subs(fxi_sym,[xs1 xs2 xs3 xs4 xs5 xs6 xs7 xs8 xs9 xs10 xs11 xs12 xs13 xs14 xs15 xs16 xs17 xs18],[xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t) xst15(t) xst16(t) xst17(t) xst18(t)]);
        xst = [xst1(t) xst2(t) xst3(t) xst4(t) xst5(t) xst6(t) xst7(t) xst8(t) xst9(t) xst10(t) xst11(t) xst12(t) xst13(t) xst14(t) xst15(t) xst16(t) xst17(t) xst18(t)];
    else
        error('System and PCE dimensions not supported')
    end
    
    fxi_sym_fun = odeFunction(fxi_sym_t,xst);


end



