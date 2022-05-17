        load('S00002.mat');
        [pln,cst] = matRad_setPlanParameters(ct,cst,'protons_HIT_APM','constRBE',0,pln);
        stf = matRad_generateStf(ct,cst,pln);

%         pln.propHeterogeneity.calcBothModes = true;
        %% ANALYTICAL
        %%% Analytical constRBE
        dij = matRad_calcParticleDose(ct,stf,pln,cst);
        resultGUI = matRad_fluenceOptimization(dij,cst,pln);      
        
        spotRemoval_cfg = MatRad_spotRemovalDij(dij,resultGUI.w);  
        resultGUI2 = spotRemoval_cfg.reoptimize(cst,pln);