clear all
clear yalmip

%Initialize structs
sys         = struct();
fns         = struct();
numsets     = struct();
numsetsRE   = struct();
filenames   = struct();

cd ..
addfolders = addpath(genpath('examples_automatica'));
addfolders = addpath(genpath('PCEprojection'));
addfolders = addpath(genpath('plottingFilesPCE'));
addfolders = addpath(genpath('resultsFiles'));
addfolders = addpath(genpath('systemFiles'));
addfolders = addpath(genpath('verificationROAPCE'));




%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%%%%%%% USER INPUT HERE %%%%%%%%%%

only_plot = 1; %set to 1 to produce plots from previously computed results files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% for ROA PCE computation

verimethod      = 'EEdegdV';   
system          = 'VDPEP';
sys.p           = 3; %PCE truncation order (p in {0,1,2,3..})

numsets.initVscale = 1e-3;   % this needs to be tuned once in the very beginning such that feasible for multiplier search! 

numsets.clean_thresh    = 1e-6;
numsets.iteration_max   = 50;
numsets.sdpsetting1     = sdpsettings('solver','mosek','verbose',0); 
numsets.sdpsetting2     = sdpsettings('solver','mosek','verbose',0); 
numsets.convCrit        = 1e-2;
numsets.gfac            = 1e-4; %

numsetsRE.clean_thresh    = 1e-6;
numsetsRE.iteration_max   = 15;
numsetsRE.sdpsetting1     = sdpsettings('solver','mosek','verbose',0); 
numsetsRE.sdpsetting2     = sdpsettings('solver','mosek','verbose',0); 
numsetsRE.convCrit        = 1e-2;


% compute the ROA results first for degree-2 V and then for degree-4 V
degVarray = [2, 4];

for jj = 1:length(degVarray)
    

    % for PCE ROA computation   
    numsets.degs.V_dU  = degVarray(jj);
    numsets.degs.V_dL  = 2;
    numsets.degs.s1_dU = 4;
    numsets.degs.s1_dL = 2;
    numsets.degs.s2_dU = 2;
    numsets.degs.s2_dL = 0;    


    % for stoch ROA computation
    numsetsRE.degs.V_dU  = numsets.degs.V_dU;
    numsetsRE.degs.V_dL  = numsets.degs.V_dL;
    numsetsRE.degs.s1_dU = 0;
    numsetsRE.degs.s1_dL = 0;
    numsetsRE.degs.s2_dU = 2;
    numsetsRE.degs.s2_dL = 0; 
    numsetsRE.degs.hi_dU = 2;
    numsetsRE.degs.hi_dL = 0; 

    numsetsRE.initQ0scale = 1e-3;


    %--------------------------------------------------------------------------
    %--------------------------------------------------------------------------


    filenames.initializationFile = strcat('initialize',system);

    filenames.initializationFile = str2func(filenames.initializationFile);
    [sys,fns] = filenames.initializationFile(sys,fns);


    filenames.resultsdirectory        = 'resultsFiles/';
    filenames.intermedresultsROAPCE   = strcat('intermediate_results_PCE_',system,'_',verimethod,'_p',num2str(sys.p),'.mat');
    filenames.finalresultsROAPCE      = strcat( filenames.resultsdirectory,'final_results_PCE_',system,'_',verimethod,'_V',num2str(numsets.degs.V_dU),'_p',num2str(sys.p),'.mat');

    % Compute PCE ROA
    veri_filename = str2func(strcat('veriROA_PCE_',verimethod));
    if only_plot == 0
        veri_filename(sys,fns,filenames,numsets);
        delete(filenames.intermedresultsROAPCE)
    end

    filenames.filenamestoplot{jj}      = strcat('final_results_',system,'_V',num2str(numsets.degs.V_dU)); %or enter each file manually
    
    
    % Compute stochastic ROA for a fixed variance on initial condition    
    veri_filename = str2func('recover_ROAstoch');
    
    if numsets.degs.V_dU == 2
        varfixarray = 0;
    else
        varfixarray = [0, 0.01, 0.03, 0.05, 0.08, 0.1];
    end

    for i = 1:length(varfixarray)
        
        sys.varfix = varfixarray(i)*eye(sys.xdim);
        
        filenames.finalresultsROAstoch = strcat(filenames.resultsdirectory,'final_results_',system,'_V',num2str(numsets.degs.V_dU),'_var',num2str(sys.varfix(1,1)),'.mat');
        
        if only_plot == 0
            veri_filename(sys,filenames,numsetsRE);
        end

    end
end


ex1_plot(filenames)










