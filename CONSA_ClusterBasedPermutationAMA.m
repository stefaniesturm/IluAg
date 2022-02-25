function CONSA_ClusterBasedPermutationAMA

%%%%%% Main script to launch  analyses
% created by Marta for SGMem, modified by Nad for SGMem2 and later modified
% by Trisia for CONSA
%% 0) DO THIS EVERY TIME YOU START DOING STH WITH THIS MAIN SCRIPT
%% Paths and codes
clear
ft_defaults
eeglab
close all

Exp_Code = '';
hd = '';

Analysis_folder = '/DATA2/BLB_EXP_202101_CONSA/Analysis/MEEG/';
scripts_folder = '/DATA2/BLB_EXP_202101_CONSA/BLB_BackUp_files/Protocols/Analysis_Scripts';
elecs_file = '/DATA2/BLB_EXP_202101_CONSA/BLB_BackUp_files/Protocols/Configuration/CONSA.asc';

addpath(genpath('/DATA2/BLB_EXP_202101_CONSA/BLB_BackUp_files/Protocols/Analysis_Scripts/plugins'));


%% 1) Export to fieldtrip for ERP analyses

session = 'LtH';

subarray = [1 2 3 4 6 7 8 9 10 12 14 15 16 17 18 19 21 22 23 24]
% participants 5, 11, 13, 20 excluded

indir = [Analysis_folder 'ClusterBased/'];

%% Don't rerun everything. If data is prepared skip the following if-statement

preparation=1 %  1 --> Data is prepared (28.05.2020); 0 --> Data is not prepared
if preparation ==0
    cfg=[];
    for subj = 1:length(subarray)
        for cond = {'highMA' 'midMA' 'lowMA' 'audhighmeanAMA' 'audmidmeanAMA' 'audlowmeanAMA' 'motorM'}%{'selfMA' 'highMA' 'midMA' 'lowMA' 'meanallMA' 'motorM'}% 'selfMA-motorM' 'highMA-motorM' 'midMA-motorM' 'lowMA-motorM'} %{'NoSound' 'eA' 'eMA' 'eM' 'tAretr' 'tMAretr' 'tAenc_F_2T' 'tAenc_R_2T' 'tMAenc_F_2T' 'tMAenc_R_2T' 'tAretr1T' 'tMAretr1T' 'Different_MA_Rem' 'Different_A_Rem' 'Same_MA_Rem' 'Same_MA_Rem' 'tAretr_R_2T' 'tAretr_F_2T' 'tMAretr_R_2T' 'tMAretr_F_2T'}
            EEG=pop_loadset([Analysis_folder num2str(subarray(subj),'%0.2d') '_' session '_' cond{:} '.set']);
            EEG= pop_selectevent(EEG, 'latency','-0.01<=0.01','deleteevents','on','deleteepochs','on','invertepochs','off');%delete duplicates if any
            %EEG=pop_rmbase(EEG,[-100 0]); %not sure why
            data=eeglab2fieldtrip(EEG,'preprocessing','none');
            % no need to include the eye channels
          
            erpdata=ft_timelockanalysis(cfg,data);
            save([indir num2str(subarray(subj),'%0.2d') '_' session '_' cond{:} '_ERP.mat'],'erpdata');
        end
    end
end

%% 2)Prepare neighbours
if preparation == 0
    %%% Do this once to get the layout and the neighbors and save them
    load([indir '10_HtL_highMA_ERP.mat']);
    cfg=[];
    cfg.elec =erpdata.elec;
    cfg.rotate=90;
    %    cfg.headshape =  ft_read_headshape([Analysis_folder num2str(subarray(subj),'%0.2d') '_' cond{:} '.set'])
    
    cfg.channel ={'all', '-HEOGR', '-HEOGL', '-VEOGU','-VEOGO','-VEOG','-HEOG','-HL1','-CB2','-CB1'}%,'-COMNT','-SCALE'}
    cfg.skipscale='yes';
    cfg.skipcomnt='yes'
    layout=ft_prepare_layout(cfg,erpdata);
    % Make sure that layout is correct
    
    cfg.layout=layout;
    
    ft_layoutplot(cfg)%,erpdata)
    
    % Prepare neighbours
    cfg=[];
    cfg.method='triangulation';
    
    cfg.layout=layout;
    
    neighbours=ft_prepare_neighbours(cfg);
    cfg.neighbours = neighbours;
    ft_neighbourplot(cfg)
    % Sanity check: confirm that the electrodes in neighbours follow the
    % same order as specified in the protocol doc (DONE: 10.06)
    cd(indir)
    save('layout.mat','layout');
    save('neighbours.mat', 'neighbours');
    cd(scripts_folder)
end

    %% 3) Perform Grand Averages of ERPs per condition
    % Do not keep individuals
    
    cfg=[]
    for cond = {'highMA' 'midMA' 'lowMA' 'audhighmeanAMA' 'audmidmeanAMA' 'audlowmeanAMA' 'highMA_motorM' 'midMA_motorM' 'lowMA_motorM'}%{'selfMA' 'highMA' 'midMA' 'lowMA' 'meanallMA' 'selfMA_motorM' 'highMA_motorM' 'midMA_motorM' 'lowMA_motorM'}%'eA' 'eMA' 'eM' 'eMA_eM' 'tAretr' 'tMAretr' 'tAenc_F_2T', 'tAenc_R_2T', ...
        %'tAretr_F_2T','tAretr_R_2T','tMAenc_F_2T_eM', 'tMAenc_R_2T_eM', 'tMAretr_F_2T', ...
        %'tMAretr_R_2T', 'NoSound' ,'tAretr1T','tMAretr1T'}%,   'Same_MA_Rem', 'Different_MA_Rem','Same_A_Rem', 'Different_A_Rem'}
        i=1;
        for subj = 1:length(subarray)
            
            load([indir num2str(subarray(subj),'%0.2d') '_' session '_' cond{:} '_ERP.mat']);
            eval(['data_' num2str(i) ' = erpdata;']);
            i = i+1;
        end
        cfg.keepindividual = 'no';
        % In this following line you should write as many data_X as
        % subjects that you are including in the GA
        eval(['GA_' session '_' cond{:} '_ERP = ft_timelockgrandaverage(cfg,data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11, data_12, data_13, data_14, data_15, data_16,data_17,data_18, data_19, data_20)']);
        eval(['save(''' indir  'GA_' session '_' cond{:} '_ERP'',''GA_' session '_' cond{:} '_ERP'')']);
        eval(['clear(''GA_' session '_' cond{:} '_ERP'')']);
        
    end
    
    % for iDiff = 1:length(condArray1)
    %
    %    load([indir 'GA_' condArray1{iDiff} '_ERP']);
    %    load([indir 'GA_' condArray2{iDiff} '_ERP']);
    %
    % cfg = [];
    % cfg.operation = 'subtract';
    % cfg.parameter = 'avg';
    % eval(['GA_' condArray1{iDiff} '_' condArray2{iDiff} '_ERP = ft_math(cfg, GA_' condArray1{iDiff} '_ERP, GA_' condArray2{iDiff} '_ERP)']);
    % save([indir 'GA_' condArray1{iDiff} '_' condArray2{iDiff} '_ERP.mat'],['GA_' condArray1{iDiff} '_' condArray2{iDiff} '_ERP']);
    %
    %
    % %this is the same as avobe, but manually, and variance should be recomputed
    % %( the variance of two random variables subtracted from each other is equal
    % %to the variances added, minus twice the covariance between them)
    % %     eval(['ERPdiff = GA_' condArray1{iDiff} '_ERP.avg - GA_' condArray2{iDiff} '_ERP.avg;']);
    % %     eval(['GA_' condArray1{iDiff} '_ERP.avg = ERPdiff;']);
    % %     eval(['ERPvardiff = GA_' condArray1{iDiff} '_ERP.var - GA_' condArray2{iDiff} '_ERP.var;']);
    % %     eval(['GA_' condArray1{iDiff} '_ERP.var = ERPvardiff;']);
    % %     eval(['GA_' condArray1{iDiff} '_' condArray2{iDiff} '_ERP = GA_' condArray1{iDiff} '_ERP;']);
    % %     save([indir 'GA_' condArray1{iDiff} '_' condArray2{iDiff} '_ERP'],['GA_' condArray1{iDiff} '_' condArray2{iDiff} '_ERP']);
    %
    % end
    
    %% Correct Motor
    %indir = 'C:\Users\Marta\Documents\BLB_EXP_201706_SGMem\Analysis\';
    %subarray = [2 3 4 5 6 7 8 9 10 11];
    
    % conditions to subtract
%     condArray1 = {'selfMA', 'highMA' ,'midMA', 'lowMA'}; %   % 1 minus 2
%     condArray2 = {'motorM', 'motorM', 'motorM' , 'motorM'};

    condArray1 = {'highMA' ,'midMA', 'lowMA'}; %   % 1 minus 2
    condArray2 = {'motorM', 'motorM' , 'motorM'};
    
    for subj = 1:length(subarray)
        for iDiff = 1:length(condArray1)
            
            load([indir num2str(subarray(subj),'%0.2d') '_' session '_' condArray1{iDiff} '_ERP']);
            cond1 = erpdata;
            load([indir num2str(subarray(subj),'%0.2d') '_' session '_' condArray2{iDiff} '_ERP']);
            cond2 = erpdata;
            
            cfg = [];
            cfg.operation = 'subtract';
            cfg.parameter = 'avg';
            
            erpdata = ft_math(cfg, cond1, cond2)
            save([indir num2str(subarray(subj),'%0.2d') '_' session '_' condArray1{iDiff} '_' condArray2{iDiff} '_ERP.mat'],'erpdata');
            
        end
    end
    
 %% 4) GRAND AVERAGES OF THE CORRECT MOTOR
    cd(indir)
   
      
    cfg=[]
    for cond =   {'highMA_motorM' 'midMA_motorM' 'lowMA_motorM'} %{'meanallMA_highMA_motorM' 'meanallMA_highMA_motorM' 'meanallMA_midMA_motorM' 'meanallMA_lowMA_motorM'}
        % 'Same_A_Rem_Same_MA_Rem', 'Different_A_Rem_Different_MA_Rem'}; %   % 1 minus 2
   
        i=1
        for subj = 1:length(subarray)
            
            load([indir num2str(subarray(subj),'%0.2d') '_' session '_' cond{:} '_ERP.mat']);
            eval(['data_' num2str(i) ' = erpdata;']);
            i = i+1;
        end
        cfg.keepindividual = 'no';
        % In this following line you should write as many data_X as
        % subjects that you are including in the GA
        eval(['GA_' session '_' cond{:} '_ERP = ft_timelockgrandaverage(cfg,data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11, data_12, data_13, data_14, data_15, data_16,data_17,data_18, data_19, data_20)']);
        eval(['save(''' indir  'GA_' session '_' cond{:} '_ERP'',''GA_' session '_' cond{:} '_ERP'')']);
        eval(['clear(''GA_' session '_' cond{:} '_ERP'')']);
        
    end

    %% Raw effect (for viz later): eA-[eMA-eM]  Difference waves
    % conditions to subtract
    
%     condArray1 = {'meanallMA' 'meanallMA' 'meanallMA' 'meanallMA'}  % by Mem for int   
%     % ,   'Same_A_Rem', 'Different_A_Rem'}; %   % 1 minus 2
%     condArray2 = {'selfMA_motorM' 'highMA_motorM' 'midMA_motorM' 'lowMA_motorM'}; % By Mem for int
%     %'Same_MA_Rem', 'Different_MA_Rem'};

    condArray1 = {'audhighmeanAMA' 'audmidmeanAMA' 'audlowmeanAMA'}  % by Mem for int   
    % ,   'Same_A_Rem', 'Different_A_Rem'}; %   % 1 minus 2
    condArray2 = {'highMA_motorM' 'midMA_motorM' 'lowMA_motorM'}; % By Mem for int
    %'Same_MA_Rem', 'Different_MA_Rem'};
    
    for subj = 1:length(subarray)
        for iDiff = 1:length(condArray1)
            
            load([indir num2str(subarray(subj),'%0.2d') '_' session '_' condArray1{iDiff} '_ERP']);
            cond1 = erpdata;
            load([indir num2str(subarray(subj),'%0.2d') '_' session '_' condArray2{iDiff} '_ERP']);
            cond2 = erpdata;
            
            cfg = [];
            cfg.operation = 'subtract';
            cfg.parameter = 'avg';
            
            raw_effect = ft_math(cfg, cond1, cond2)
            save([indir num2str(subarray(subj),'%0.2d') '_' session '_' condArray1{iDiff} '_' condArray2{iDiff} '_ERP.mat'],'raw_effect');
            
        end
    end
    

    %% 4) GRAND AVERAGES OF THE DIFFERENCE WAVES
    cd(indir)
   
      
    cfg=[]
    for cond =   {'audhighmeanAMA_highMA_motorM' 'audmidmeanAMA_midMA_motorM' 'audlowmeanAMA_lowMA_motorM'} %{'meanallMA_highMA_motorM' 'meanallMA_highMA_motorM' 'meanallMA_midMA_motorM' 'meanallMA_lowMA_motorM'}
        % 'Same_A_Rem_Same_MA_Rem', 'Different_A_Rem_Different_MA_Rem'}; %   % 1 minus 2
   
        i=1
        for subj = 1:length(subarray)
            
            load([indir num2str(subarray(subj),'%0.2d') '_' session '_' cond{:} '_ERP.mat']);
            eval(['data_' num2str(i) ' = raw_effect;']);
            i = i+1;
        end
        cfg.keepindividual = 'no';
        % In this following line you should write as many data_X as
        % subjects that you are including in the GA
        eval(['GA_' session '_' cond{:} '_ERP = ft_timelockgrandaverage(cfg,data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11, data_12, data_13, data_14, data_15, data_16,data_17,data_18, data_19, data_20)']);
        eval(['save(''' indir  'GA_' session '_' cond{:} '_ERP'',''GA_' session '_' cond{:} '_ERP'')']);
        eval(['clear(''GA_' session '_' cond{:} '_ERP'')']);
        
    end
  
    %% 4) ERP PLOTS
cd(indir)
session = 'LtH';
cond1 = 'audhighmeanAMA';%'meanallMA';
cond2 = 'highMA_motorM';
load layout
eval(['load ' 'GA_' session '_' cond1 '_ERP'])
eval(['load ' 'GA_' session '_' cond2 '_ERP'])
cfg=[];
cfg.layout=layout;
cfg.xlim=[-0.1 0.5];
cfg.ylim=[-4 4];
figure;
eval(['ft_multiplotER(cfg,GA_' session '_' cond1 '_ERP,' 'GA_' session '_' cond2 '_ERP)'])

%% STATISTICS
    % --------------------------------
    % PREPARE DATA CELLS 
    % --------------------------------
%  T-TESTS

% Conditions
cd(indir)
    %Normal
    HtL_High_MA = cell(length(subarray),1);
    HtL_Mid_MA = cell(length(subarray),1);
    HtL_Low_MA = cell(length(subarray),1);
    
    HtL_High_A = cell(length(subarray),1);
    HtL_Mid_A = cell(length(subarray),1);
    HtL_Low_A = cell(length(subarray),1);
    
    LtH_High_MA = cell(length(subarray),1);
    LtH_Mid_MA = cell(length(subarray),1);
    LtH_Low_MA = cell(length(subarray),1);
    
    LtH_High_A = cell(length(subarray),1);
    LtH_Mid_A = cell(length(subarray),1);
    LtH_Low_A = cell(length(subarray),1);
    
    %Difference Waves
    HtL_Diff_High = cell(length(subarray),1);
    HtL_Diff_Mid = cell(length(subarray),1);
    HtL_Diff_Low = cell(length(subarray),1);
    
    LtH_Diff_High = cell(length(subarray),1);
    LtH_Diff_Mid = cell(length(subarray),1);
    LtH_Diff_Low = cell(length(subarray),1);
    
% CONDITIONS
for cond = {'HtL_High_MA' 'HtL_Mid_MA' 'HtL_Low_MA' 'LtH_High_MA' 'LtH_Mid_MA' 'LtH_Low_MA' 'HtL_High_A' 'HtL_Mid_A' 'HtL_Low_A' 'LtH_High_A' 'LtH_Mid_A' 'LtH_Low_A'}
    switch cond{:}
        case 'HtL_High_MA'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_highMA_motorM_ERP.mat']);
                HtL_High_MA{i}=erpdata;
                i = i+1
            end
        case 'HtL_Mid_MA'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_midMA_motorM_ERP.mat']);
                HtL_Mid_MA{i}=erpdata;
                i = i+1
            end
        case 'HtL_Low_MA'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_lowMA_motorM_ERP.mat']);
                HtL_Low_MA{i}=erpdata;
                i = i+1
            end
        case 'LtH_High_MA'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_highMA_motorM_ERP.mat']);
                LtH_High_MA{i}=erpdata;
                i = i+1
            end
        case 'LtH_Mid_MA'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_midMA_motorM_ERP.mat']);
                LtH_Mid_MA{i}=erpdata;
                i = i+1
            end
        case 'LtH_Low_MA'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_lowMA_motorM_ERP.mat']);
                LtH_Low_MA{i}=erpdata;
                i = i+1
            end
            
        case 'HtL_High_A'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_audhighmeanAMA_ERP.mat']);
                HtL_High_A{i}=erpdata;
                i = i+1
            end
        case 'HtL_Mid_A'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_audmidmeanAMA_ERP.mat']);
                HtL_Mid_A{i}=erpdata;
                i = i+1
            end
        case 'HtL_Low_A'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_audlowmeanAMA_ERP.mat']);
                HtL_Low_A{i}=erpdata;
                i = i+1
            end
        case 'LtH_High_A'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_audhighmeanAMA_ERP.mat']);
                LtH_High_A{i}=erpdata;
                i = i+1
            end
        case 'LtH_Mid_A'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_audmidmeanAMA_ERP.mat']);
                LtH_Mid_A{i}=erpdata;
                i = i+1
            end
        case 'LtH_Low_A'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_audlowmeanAMA_ERP.mat']);
                LtH_Low_A{i}=erpdata;
                i = i+1
            end
    end
     
        cd(indir)
        eval(['save ' 'HtL_High_MA HtL_High_MA']);
        eval(['save ' 'HtL_Mid_MA HtL_Mid_MA']);
        eval(['save ' 'HtL_Low_MA HtL_Low_MA']);

        eval(['save ' 'LtH_High_MA LtH_High_MA']);
        eval(['save ' 'LtH_Mid_MA LtH_Mid_MA']);
        eval(['save ' 'LtH_Low_MA LtH_Low_MA']);

        eval(['save ' 'HtL_High_A HtL_High_A']);
        eval(['save ' 'HtL_Mid_A HtL_Mid_A']);
        eval(['save ' 'HtL_Low_A HtL_Low_A']);

        eval(['save ' 'LtH_High_A LtH_High_A']);
        eval(['save ' 'LtH_Mid_A LtH_Mid_A']);
        eval(['save ' 'LtH_Low_A LtH_Low_A']);


end
 
% DIFFERENCE WAVES
for cond = {'HtL_Diff_High' 'HtL_Diff_Mid' 'HtL_Diff_Low' 'LtH_Diff_High' 'LtH_Diff_Mid' 'LtH_Diff_Low'}
    switch cond{:}
        case 'HtL_Diff_High'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_audhighmeanAMA_highMA_motorM_ERP.mat']);
                HtL_Diff_High{i}=raw_effect;
                i = i+1
            end
        case 'HtL_Diff_Mid'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_audmidmeanAMA_midMA_motorM_ERP.mat']);
                HtL_Diff_Mid{i}=raw_effect;
                i = i+1
            end
        case 'HtL_Diff_Low'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_audlowmeanAMA_lowMA_motorM_ERP.mat']);
                HtL_Diff_Low{i}=raw_effect;
                i = i+1
            end
        case 'LtH_Diff_High'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_audhighmeanAMA_highMA_motorM_ERP.mat']);
                LtH_Diff_High{i}=raw_effect;
                i = i+1
            end
        case 'LtH_Diff_Mid'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_audmidmeanAMA_midMA_motorM_ERP.mat']);
                LtH_Diff_Mid{i}=raw_effect;
                i = i+1
            end
        case 'LtH_Diff_Low'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_audlowmeanAMA_lowMA_motorM_ERP.mat']);
                LtH_Diff_Low{i}=raw_effect;
                i = i+1
            end
    end

        cd(indir)
        eval(['save ' 'HtL_Diff_High HtL_Diff_High']);
        eval(['save ' 'HtL_Diff_Mid HtL_Diff_Mid']);
        eval(['save ' 'HtL_Diff_Low HtL_Diff_Low']);

        
        eval(['save ' 'LtH_Diff_High LtH_Diff_High']);
        eval(['save ' 'LtH_Diff_Mid LtH_Diff_Mid']);
        eval(['save ' 'LtH_Diff_Low LtH_Diff_Low']);

end

%% % We perform the stats CONDITION
cd(indir)

load layout
load neighbours

cfg = [];
cfg.parameter        = 'avg'; 
cfg.latency          = [0 0.5];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; 
cfg.clusterstatistic = 'maxsum';
cfg.minnbchan        = 2; 
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05; 
cfg.numrandomization = 10000; % Pump up
cfg.correcttail = 'prob';

cfg.layout = layout; 
cfg.neighbours = neighbours; 

% nº of levels of the other factor (there is no other factor so 1) * length(subarray)
subj = length(subarray);
% nº of variables (2 [Self or High or Mid or Low]) * subj
design = zeros(2,2*subj);
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i)=i;
end
design(2,1:subj)=1; design(2,subj+1:2*subj)=2;
cfg.design = design;
cfg.ivar=2; cfg.uvar=1;

% HERE YOU WRITE THE NAME YOU WANT OF THE STATS DATA STRUCTURE
% HtL_High_MAvsA
[ERP_HtL_High_MAvsA] = ft_timelockstatistics(cfg, HtL_High_MA{:}, HtL_High_A{:}); 
save ([indir '/ERP_HtL_High_MAvsA.mat'], 'ERP_HtL_High_MAvsA')
% HtL_Mid_MAvsA
[ERP_HtL_Mid_MAvsA] = ft_timelockstatistics(cfg, HtL_Mid_MA{:}, HtL_Mid_A{:}); 
save ([indir '/ERP_HtL_Mid_MAvsA.mat'], 'ERP_HtL_Mid_MAvsA')
% HtL_Low_MAvsA
[ERP_HtL_Low_MAvsA] = ft_timelockstatistics(cfg, HtL_Low_MA{:}, HtL_Low_A{:}); 
save ([indir '/ERP_HtL_Low_MAvsA.mat'], 'ERP_HtL_Low_MAvsA')

% LtH_High_MAvsA
[ERP_LtH_High_MAvsA] = ft_timelockstatistics(cfg, LtH_High_MA{:}, LtH_High_A{:}); 
save ([indir '/ERP_LtH_High_MAvsA.mat'], 'ERP_LtH_High_MAvsA')
% LtH_Mid_MAvsA
[ERP_LtH_Mid_MAvsA] = ft_timelockstatistics(cfg, LtH_Mid_MA{:}, LtH_Mid_A{:}); 
save ([indir '/ERP_LtH_Mid_MAvsA.mat'], 'ERP_LtH_Mid_MAvsA')
% LtH_Low_MAvsA
[ERP_LtH_Low_MAvsA] = ft_timelockstatistics(cfg, LtH_Low_MA{:}, LtH_Low_A{:}); 
save ([indir '/ERP_LtH_Low_MAvsA.mat'], 'ERP_LtH_Low_MAvsA')



% HtL_Diff_HighvsMid
[ERP_HtL_Diff_HighvsMid] = ft_timelockstatistics(cfg, HtL_Diff_High{:}, HtL_Diff_Mid{:}); 
save ([indir '/ERP_HtL_Diff_HighvsMid.mat'], 'ERP_HtL_Diff_HighvsMid')
% HtL_Diff_HighvsLow
[ERP_HtL_Diff_HighvsLow] = ft_timelockstatistics(cfg, HtL_Diff_High{:}, HtL_Diff_Low{:}); 
save ([indir '/ERP_HtL_Diff_HighvsLow.mat'], 'ERP_HtL_Diff_HighvsLow')
% HtL_Diff_MidvsLow
[ERP_HtL_Diff_MidvsLow] = ft_timelockstatistics(cfg, HtL_Diff_Mid{:}, HtL_Diff_Low{:}); 
save ([indir '/ERP_HtL_Diff_MidvsLow.mat'], 'ERP_HtL_Diff_MidvsLow')

% LtH_Diff_HighvsMid
[ERP_LtH_Diff_HighvsMid] = ft_timelockstatistics(cfg, LtH_Diff_High{:}, LtH_Diff_Mid{:}); 
save ([indir '/ERP_LtH_Diff_HighvsMid.mat'], 'ERP_LtH_Diff_HighvsMid')
% LtH_Diff_HighvsLow
[ERP_LtH_Diff_HighvsLow] = ft_timelockstatistics(cfg, LtH_Diff_High{:}, LtH_Diff_Low{:}); 
save ([indir '/ERP_LtH_Diff_HighvsLow.mat'], 'ERP_LtH_Diff_HighvsLow')
% LtH_Diff_MidvsLow
[ERP_LtH_Diff_MidvsLow] = ft_timelockstatistics(cfg, LtH_Diff_Mid{:}, LtH_Diff_Low{:}); 
save ([indir '/ERP_LtH_Diff_MidvsLow.mat'], 'ERP_LtH_Diff_MidvsLow')
%% PLOTS T-TESTS
cd(indir)
load layout

load ERP_HtL_High_MAvsA
load ERP_HtL_Mid_MAvsA
load ERP_HtL_Low_MAvsA

load ERP_LtH_High_MAvsA
load ERP_LtH_Mid_MAvsA
load ERP_LtH_Low_MAvsA



load ERP_HtL_Diff_HighvsMid
load ERP_HtL_Diff_HighvsLow
load ERP_HtL_Diff_MidvsLow


load ERP_LtH_Diff_HighvsMid
load ERP_LtH_Diff_HighvsLow
load ERP_LtH_Diff_MidvsLow

cfg=[];
cfg.layout=layout;
cfg.highlightsizeseries = [15 15 15 15 15];
cfg.highlightsymbolseries=['.','.','.','.','.'];
cfg.zlim='maxabs';
cfg.alpha=0.05;


% HtL_High_MA_A
ft_clusterplot(cfg,ERP_HtL_High_MAvsA);
% HtL_Mid_MA_A
ft_clusterplot(cfg,ERP_HtL_Mid_MAvsA);
% HtL_Low_MA_A
ft_clusterplot(cfg,ERP_HtL_Low_MAvsA);


% LtH_High_MA_A
ft_clusterplot(cfg,ERP_LtH_High_MAvsA);
% LtH_Mid_MA_A
ft_clusterplot(cfg,ERP_LtH_Mid_MAvsA);
% LtH_Low_MA_A
ft_clusterplot(cfg,ERP_LtH_Low_MAvsA);




% HtL_Diff_High_Mid
ft_clusterplot(cfg,ERP_HtL_Diff_HighvsMid);
% HtL_Diff_High_Low
ft_clusterplot(cfg,ERP_HtL_Diff_HighvsLow);
% HtL_Diff_Mid_Low
ft_clusterplot(cfg,ERP_HtL_Diff_MidvsLow);


% LtH_Diff_High_Mid
ft_clusterplot(cfg,ERP_LtH_Diff_HighvsMid);
% LtH_Diff_High_Low
ft_clusterplot(cfg,ERP_LtH_Diff_HighvsLow);
% LtH_Diff_Mid_Low
ft_clusterplot(cfg,ERP_LtH_Diff_MidvsLow);

 end