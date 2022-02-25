function CONSA_ClusterBasedPermutation

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

session = 'HtL';

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
    
    %% Raw effect (for viz later): HtL-LtH  Difference waves (for doing the ANOVA [if you substract both you can do ANOVA 1x4 instead ANOVA 2x4])
    % conditions to subtract
    
    condArray1 = {'HtL_meanallMA_selfMA_motorM' 'HtL_meanallMA_highMA_motorM' 'HtL_meanallMA_midMA_motorM' 'HtL_meanallMA_lowMA_motorM'}  % by Mem for int   
    % ,   'Same_A_Rem', 'Different_A_Rem'}; %   % 1 minus 2
    condArray2 = {'LtH_meanallMA_selfMA_motorM' 'LtH_meanallMA_highMA_motorM' 'LtH_meanallMA_midMA_motorM' 'LtH_meanallMA_lowMA_motorM'}; % By Mem for int
    %'Same_MA_Rem', 'Different_MA_Rem'};
    
    for subj = 1:length(subarray)
        for iDiff = 1:length(condArray1)
            
            load([indir num2str(subarray(subj),'%0.2d') '_' condArray1{iDiff} '_ERP']);
            cond1 = raw_effect;
            load([indir num2str(subarray(subj),'%0.2d') '_' condArray2{iDiff} '_ERP']);
            cond2 = raw_effect;
            
            cfg = [];
            cfg.operation = 'subtract';
            cfg.parameter = 'avg';
            
            raw_effect = ft_math(cfg, cond1, cond2)
            save([indir num2str(subarray(subj),'%0.2d') '_' condArray1{iDiff} '_' condArray2{iDiff} '_ERP.mat'],'raw_effect');
            
        end
    end 
    
     %% GRAND AVERAGES OF THE RAW EFFECT
    cd(indir)
   
      
    cfg=[]
    for cond =   {'HtL_meanallMA_selfMA_motorM_LtH_meanallMA_selfMA_motorM' 'HtL_meanallMA_highMA_motorM_LtH_meanallMA_highMA_motorM' 'HtL_meanallMA_midMA_motorM_LtH_meanallMA_midMA_motorM' 'HtL_meanallMA_lowMA_motorM_LtH_meanallMA_lowMA_motorM'}
        % 'Same_A_Rem_Same_MA_Rem', 'Different_A_Rem_Different_MA_Rem'}; %   % 1 minus 2
   
        i=1
        for subj = 1:length(subarray)
            
            load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
            eval(['data_' num2str(i) ' = raw_effect;']);
            i = i+1;
        end
        cfg.keepindividual = 'no';
        % In this following line you should write as many data_X as
        % subjects that you are including in the GA
        eval(['GA_' cond{:} '_ERP = ft_timelockgrandaverage(cfg,data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11, data_12, data_13, data_14, data_15, data_16,data_17,data_18, data_19, data_20)']);
        eval(['save(''' indir  'GA_' cond{:} '_ERP'',''GA_' cond{:} '_ERP'')']);
        eval(['clear(''GA_' cond{:} '_ERP'')']);
        
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
% %%%%Difference waves
    
    % Main effect of condition
%     Self = cell(length(subarray),1);
    High = cell(length(subarray),1);
    Mid = cell(length(subarray),1);
    Low = cell(length(subarray),1);
   
    % Main effect of session
    HtL = cell(length(subarray),1);
    LtH= cell(length(subarray),1);
  
    % Interaction
%     Diff_Self = cell(length(subarray),1);
    Diff_High = cell(length(subarray),1);
    Diff_Mid = cell(length(subarray),1);
    Diff_Low = cell(length(subarray),1);
    
%% MAIN EFFECT OF SESSION (for ANOVA)
 for cond = {'HtL' 'LtH'}
        switch cond{:}
            case 'HtL'
                i=1;
                for subj = 1:length(subarray)
                    load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_audhighmeanAMA_highMA_motorM_ERP.mat']);
                    HtL{i}=raw_effect;
                    i = i+1;
                    load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_audmidmeanAMA_midMA_motorM_ERP.mat']);
                    HtL{i}=raw_effect;
                    i = i+1
                    load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_audlowmeanAMA_lowMA_motorM_ERP.mat']);
                    HtL{i}=raw_effect;
                    i = i+1;
                end
            case 'LtH'
                i=1;
                for subj = 1:length(subarray)
                    load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_audhighmeanAMA_highMA_motorM_ERP.mat']);
                    HtL{i}=raw_effect;
                    i = i+1;
                    load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_audmidmeanAMA_midMA_motorM_ERP.mat']);
                    HtL{i}=raw_effect;
                    i = i+1
                    load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_audlowmeanAMA_lowMA_motorM_ERP.mat']);
                    HtL{i}=raw_effect;
                    i = i+1;
                end
        end
     
  
        cd(indir)
        eval(['save ' 'HtL HtL']);
        eval(['save ' 'LtH LtH']);
 end
   
%  for cond = {'HtL' 'LtH'}
%         switch cond{:}
%             case 'HtL'
%                 i=1;
%                 for subj = 1:length(subarray)
%                     load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_meanallMA_selfMA_motorM_ERP.mat']);
%                     HtL{i}=raw_effect;
%                     i = i+1
%                     load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_meanallMA_highMA_motorM_ERP.mat']);
%                     HtL{i}=raw_effect;
%                     i = i+1;
%                     load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_meanallMA_midMA_motorM_ERP.mat']);
%                     HtL{i}=raw_effect;
%                     i = i+1
%                     load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_meanallMA_lowMA_motorM_ERP.mat']);
%                     HtL{i}=raw_effect;
%                     i = i+1;
%                 end
%             case 'LtH'
%                 i=1;
%                 for subj = 1:length(subarray)
%                     load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_meanallMA_selfMA_motorM_ERP.mat']);
%                     LtH{i}=raw_effect;
%                     i = i+1
%                     load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_meanallMA_highMA_motorM_ERP.mat']);
%                     LtH{i}=raw_effect;
%                     i = i+1;
%                     load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_meanallMA_midMA_motorM_ERP.mat']);
%                     LtH{i}=raw_effect;
%                     i = i+1
%                     load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_meanallMA_lowMA_motorM_ERP.mat']);
%                     LtH{i}=raw_effect;
%                     i = i+1;
%                 end
%         end
%      
%   
%         cd(indir)
%         eval(['save ' 'HtL HtL']);
%         eval(['save ' 'LtH LtH']);
%  end
%    
     %% Grand average for the main effects of SESSION (for ANOVA)
     cd(indir)
     
     cfg=[]
    for cond =   {'HtL' 'LtH'}
        
        eval(['load ' cond{:}])
        
        eval(['tmp =' cond{:}])
        i=1
        
        for subj = 1:length(subarray)*4 % times four because they are main effecs here 
            % and 4 conditions of the other variable inside each condition
            tmp2 = tmp{subj};
            eval(['data_' num2str(i) ' = tmp2;']);
            i = i+1;
        end
        cfg.keepindividual = 'no';
        % In this following line you should write as many data_X as
        % subjects that you are including in the GA
        eval(['GA_' cond{:} '_ERP = ft_timelockgrandaverage(cfg,data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11, data_12, data_13, data_14, data_15, data_16,data_17,data_18, data_19, data_20)']);
        eval(['save(''' indir  'GA_' cond{:} '_ERP'',''GA_' cond{:} '_ERP'')']);
        eval(['clear(''GA_' cond{:} '_ERP'')']);
        
    end

%% % We perform the stats of Main Effects of SESSION (ANOVA)
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

% nº of levels of the other factor (4 [Self High Mid Low]) * length(subarray)
subj = 4*length(subarray);
% nº of factors of the independent variable (2 [Htl & LtH]) * subj
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
[ERP_HtLvsLtH] = ft_timelockstatistics(cfg, HtL{:}, LtH{:}); 
save ([indir '/ERP_HtLvsLtH.mat'], 'ERP_HtLvsLtH')

%% STATS PLOTS SESSION (ANOVA)
cd(indir)
load layout
load ERP_HtLvsLtH
cfg=[];
cfg.layout=layout;
cfg.highlightsizeseries = [15 15 15 15 15];
cfg.highlightsymbolseries=['.','.','.','.','.'];
cfg.zlim='maxabs';
cfg.alpha=0.05;
ft_clusterplot(cfg,ERP_HtLvsLtH);

%% MAIN EFFECT OF CONDITION (for ANOVA)
 for cond = {'meanallMA_selfMA_motorM' 'meanallMA_highMA_motorM' 'meanallMA_midMA_motorM' 'meanallMA_lowMA_motorM'}
        switch cond{:}
            case 'meanallMA_selfMA_motorM'
                i=1;
                for subj = 1:length(subarray)
                    load([indir num2str(subarray(subj),'%0.2d') '_HtL_' cond{:} '_ERP.mat']);
                    Self{i}=raw_effect;
                    i = i+1
                    load([indir num2str(subarray(subj),'%0.2d') '_LtH_' cond{:} '_ERP.mat']);
                    Self{i}=raw_effect;
                    i = i+1;
                end
            case 'meanallMA_highMA_motorM'
                i=1;
                for subj = 1:length(subarray)
                    load([indir num2str(subarray(subj),'%0.2d') '_HtL_' cond{:} '_ERP.mat']);
                    High{i}=raw_effect;
                    i = i+1
                    load([indir num2str(subarray(subj),'%0.2d') '_LtH_' cond{:} '_ERP.mat']);
                    High{i}=raw_effect; 
                    i = i+1;
                end
            case 'meanallMA_midMA_motorM'
                i=1;
                for subj = 1:length(subarray)
                    load([indir num2str(subarray(subj),'%0.2d') '_HtL_' cond{:} '_ERP.mat']);
                    Mid{i}=raw_effect;
                    i = i+1
                    load([indir num2str(subarray(subj),'%0.2d') '_LtH_' cond{:} '_ERP.mat']);
                    Mid{i}=raw_effect;   
                    i = i+1;
                end
            case 'meanallMA_lowMA_motorM'
                i=1;
                for subj = 1:length(subarray)
                    load([indir num2str(subarray(subj),'%0.2d') '_HtL_' cond{:} '_ERP.mat']);
                    Low{i}=raw_effect;
                    i = i+1
                    load([indir num2str(subarray(subj),'%0.2d') '_LtH_' cond{:} '_ERP.mat']);
                    Low{i}=raw_effect;     
                    i = i+1;
                end
        end
     
  
        cd(indir)
        eval(['save ' 'Self Self']);
        eval(['save ' 'High High']);
        eval(['save ' 'Mid Mid']);
        eval(['save ' 'Low Low']);
 end
   
     %% Grand average for the main effects of CONDITION (for ANOVA)
     cd(indir)
     
     cfg=[]
    for cond =   {'Self' 'High' 'Mid' 'Low'}
        
        eval(['load ' cond{:}])
        
        eval(['tmp =' cond{:}])
        i=1
        
        for subj = 1:length(subarray)*2 % times two because they are main effecs here
            tmp2 = tmp{subj};
            eval(['data_' num2str(i) ' = tmp2;']);
            i = i+1;
        end
        cfg.keepindividual = 'no';
        % In this following line you should write as many data_X as
        % subjects that you are including in the GA
        eval(['GA_' cond{:} '_ERP = ft_timelockgrandaverage(cfg,data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11, data_12, data_13, data_14, data_15, data_16,data_17,data_18, data_19, data_20)']);
        eval(['save(''' indir  'GA_' cond{:} '_ERP'',''GA_' cond{:} '_ERP'')']);
        eval(['clear(''GA_' cond{:} '_ERP'')']);
        
    end

%% % We perform the stats of Main Effects of CONDITION (ANOVA)
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

load Self
load High
load Mid
load Low
% Cluster for One way anova 1T sequences
cfg.statistic='depsamplesFmultivariate';

% nº of levels of the other factor (2 [HtL & LtH]) * length(subarray)
subj = 2*length(subarray);
% nº of factors of the independent variable (4 [Self High Mid Low]) * subj
design = zeros(2,4*subj);
% One per factor
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i)=i;
end
for i = 1:subj
    design(1,subj*2+i)=i;
end
for i = 1:subj
    design(1,subj*3+i)=i;
end

design(2,1:subj)=1; design(2,subj+1:2*subj)=2; design(2,2*subj+1:3*subj)=3; design(2,3*subj+1:4*subj)=4;
cfg.design = design;
cfg.ivar=2; cfg.uvar=1;

cfg.clusterthreshold='nonparametric';
ERP_oneway_condition=ft_timelockstatistics(cfg, Self{:}, High{:}, Mid{:}, Low{:});
save ([indir 'ERP_oneway_condition.mat'], 'ERP_oneway_condition')

% %%    
%     session = 'HtL';
%     
%     for cond = {'selfMA' 'highMA' 'midMA' 'lowMA'}
%         switch cond{:}
%             case 'selfMA'
%                 i=1;
%                 for subj = 1:length(subarray)
%                     load([indir num2str(subarray(subj),'%0.2d') '_' session '_' cond{:} '_ERP.mat']);
%                     eval(['Self_' session '{i}=erpdata']);
%                     i = i+1
%                     
%                 end
%             case 'highMA'
%                 i=1;
%                 for subj = 1:length(subarray)
%                     load([indir num2str(subarray(subj),'%0.2d') '_' session '_' cond{:} '_ERP.mat']);
%                     eval(['High_' session '{i}=erpdata']);
%                     i = i+1
%                     
%                 end
%             case 'midMA'
%                 i=1;
%                 for subj = 1:length(subarray)
%                     load([indir num2str(subarray(subj),'%0.2d') '_' session '_' cond{:} '_ERP.mat']);
%                     eval(['Mid_' session '{i}=erpdata']);
%                     i = i+1
%                     
%                 end
%             case 'lowMA'
%                 i=1;
%                 for subj = 1:length(subarray)
%                     load([indir num2str(subarray(subj),'%0.2d') '_' session '_' cond{:} '_ERP.mat']);
%                     eval(['Low_' session '{i}=erpdata']);
%                     i = i+1
%                     
%                 end
% %             % Sound ME
% %             case 'remembered_enc'
% %                 i = 1;
% %                 for subj = 1:length(subarray)
% %                     load([indir num2str(subarray(subj),'%0.2d') '_' session '_tMAenc_R_2T_eM_ERP.mat']);
% %                     remembered_enc{i}=erpdata;
% %                     i = i+1
% %                     load([indir num2str(subarray(subj),'%0.2d') '_' session '_tAenc_R_2T_ERP.mat']);
% %                     remembered_enc{i}=erpdata;
% %                     
% %                     i = i+1;
% %                 end
% %             case 'forgotten_enc'
% %                  i = 1;
% %                 for subj = 1:length(subarray)
% %                     load([indir num2str(subarray(subj),'%0.2d') '_' session '_tMAenc_F_2T_eM_ERP.mat']);
% %                     forgotten_enc{i}=erpdata;
% %                     i = i+1
% %                     load([indir num2str(subarray(subj),'%0.2d') '_' session '_tAenc_F_2T_ERP.mat']);
% %                     forgotten_enc{i}=erpdata;
% %                     
% %                     i = i+1;
% %                 end
%         end
%      
%   
%         cd(indir)
%         eval(['save ' session '_self' ' ' 'Self_' session]);
%         eval(['save ' session '_high' ' ' 'High_' session]);
%         eval(['save ' session '_mid' ' ' 'Mid_' session]);
%         eval(['save ' session '_low' ' ' 'Low_' session]);
%     end

%     %% Grand average for the main effects 
%      cd(indir)
%    
%       load all_dataSetsForClusters
%     cfg=[]
%     for cond =   {'tAenc' 'tMAenc_eM' 'tMAretr_ME_sound' 'tAretr_ME_sound' 'remembered_enc' 'forgotten_enc' 'remembered_retr' 'forgotten_retr' }
%         
%         eval(['tmp =' cond{:}])
%         i=1
%         
%         for subj = 1:length(subarray)*2 % times two because they are main effecs here
%             tmp2 = tmp{subj};
%             eval(['data_' num2str(i) ' = tmp2;']);
%             i = i+1;
%         end
%         cfg.keepindividual = 'no';
%         % In this following line you should write as many data_X as
%         % subjects that you are including in the GA
%         eval(['GA_' session '_' cond{:} '_ERP = ft_timelockgrandaverage(cfg,data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11, data_12, data_13, data_14, data_15, data_16,data_17,data_18, data_19, data_20)']);
%         eval(['save(''' indir  'GA_' session '_' cond{:} '_ERP'',''GA_' session '_' cond{:} '_ERP'')']);
%         eval(['clear(''GA_' session '_' cond{:} '_ERP'')']);
%         
%     end

%% STATS PLOTS CONDITION
cd(indir)
load layout
load ERP_oneway_condition
cfg=[];
cfg.layout=layout;
cfg.highlightsizeseries = [15 15 15 15 15];
cfg.highlightsymbolseries=['.','.','.','.','.'];
cfg.zlim='maxabs';
cfg.alpha=0.05;
ft_clusterplot(cfg,ERP_oneway_condition);

%%  POST-HOC (ANOVA)
% Conditions
cd(indir)
    HtL_Self = cell(length(subarray),1);
    HtL_High = cell(length(subarray),1);
    HtL_Mid = cell(length(subarray),1);
    HtL_Low = cell(length(subarray),1);
    
    LtH_Self = cell(length(subarray),1);
    LtH_High = cell(length(subarray),1);
    LtH_Mid = cell(length(subarray),1);
    LtH_Low = cell(length(subarray),1);
    
% CONDITIONS
for cond = {'HtL_Self' 'HtL_High' 'HtL_Mid' 'HtL_Low' 'LtH_Self' 'LtH_High' 'LtH_Mid' 'LtH_Low'}
    switch cond{:}
        case 'HtL_Self'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_meanallMA_selfMA_motorM_ERP.mat']);
                HtL_Self{i}=raw_effect;
                i = i+1
            end
        case 'HtL_High'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_meanallMA_highMA_motorM_ERP.mat']);
                HtL_High{i}=raw_effect;
                i = i+1
            end
        case 'HtL_Mid'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_meanallMA_midMA_motorM_ERP.mat']);
                HtL_Mid{i}=raw_effect;
                i = i+1
            end
        case 'HtL_Low'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'HtL_meanallMA_lowMA_motorM_ERP.mat']);
                HtL_Low{i}=raw_effect;
                i = i+1
            end
        case 'LtH_Self'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_meanallMA_selfMA_motorM_ERP.mat']);
                LtH_Self{i}=raw_effect;
                i = i+1
            end
        case 'LtH_High'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_meanallMA_highMA_motorM_ERP.mat']);
                LtH_High{i}=raw_effect;
                i = i+1
            end
        case 'LtH_Mid'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_meanallMA_midMA_motorM_ERP.mat']);
                LtH_Mid{i}=raw_effect;
                i = i+1
            end
        case 'LtH_Low'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' 'LtH_meanallMA_lowMA_motorM_ERP.mat']);
                LtH_Low{i}=raw_effect;
                i = i+1
            end
        end
     
  
%         cd(indir)
%         eval(['save ' 'HtL_Self HtL_Self']);
%         eval(['save ' 'HtL_High HtL_High']);
%         eval(['save ' 'HtL_Mid HtL_Mid']);
%         eval(['save ' 'HtL_Low HtL_Low']);
%         
%         eval(['save ' 'LtH_Self LtH_Self']);
%         eval(['save ' 'LtH_High LtH_High']);
%         eval(['save ' 'LtH_Mid LtH_Mid']);
%         eval(['save ' 'LtH_Low LtH_Low']);
 end

%% % We perform the stats of Main Effects of SESSION (ANOVA)
    %----------------------------------
    % Post hocs condition
    %----------------------------------
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
% HtL_Self_High
[ERP_HtL_SelfvsHigh] = ft_timelockstatistics(cfg, HtL_Self{:}, HtL_High{:}); 
save ([indir '/ERP_HtL_SelfvsHigh.mat'], 'ERP_HtL_SelfvsHigh')
% HtL_Self_Mid
[ERP_HtL_SelfvsMid] = ft_timelockstatistics(cfg, HtL_Self{:}, HtL_Mid{:}); 
save ([indir '/ERP_HtL_SelfvsMid.mat'], 'ERP_HtL_SelfvsMid')
% HtL_Self_Low
[ERP_HtL_SelfvsLow] = ft_timelockstatistics(cfg, HtL_Self{:}, HtL_Low{:}); 
save ([indir '/ERP_HtL_SelfvsLow.mat'], 'ERP_HtL_SelfvsLow')
% HtL_High_Mid
[ERP_HtL_HighvsMid] = ft_timelockstatistics(cfg, HtL_High{:}, HtL_Mid{:}); 
save ([indir '/ERP_HtL_HighvsMid.mat'], 'ERP_HtL_HighvsMid')
% HtL_High_low
[ERP_HtL_HighvsLow] = ft_timelockstatistics(cfg, HtL_High{:}, HtL_Low{:}); 
save ([indir '/ERP_HtL_HighvsLow.mat'], 'ERP_HtL_HighvsLow')
% HtL_Mid_Low
[ERP_HtL_MidvsLow] = ft_timelockstatistics(cfg, HtL_Mid{:}, HtL_Low{:}); 
save ([indir '/ERP_HtL_MidvsLow.mat'], 'ERP_HtL_MidvsLow')

% LtH_Self_High
[ERP_LtH_SelfvsHigh] = ft_timelockstatistics(cfg, LtH_Self{:}, LtH_High{:}); 
save ([indir '/ERP_LtH_SelfvsHigh.mat'], 'ERP_LtH_SelfvsHigh')
% LtH_Self_Mid
[ERP_LtH_SelfvsMid] = ft_timelockstatistics(cfg, LtH_Self{:}, LtH_Mid{:}); 
save ([indir '/ERP_LtH_SelfvsMid.mat'], 'ERP_LtH_SelfvsMid')
% LtH_Self_Low
[ERP_LtH_SelfvsLow] = ft_timelockstatistics(cfg, LtH_Self{:}, LtH_Low{:}); 
save ([indir '/ERP_LtH_SelfvsLow.mat'], 'ERP_LtH_SelfvsLow')
% LtH_High_Mid
[ERP_LtH_HighvsMid] = ft_timelockstatistics(cfg, LtH_High{:}, LtH_Mid{:}); 
save ([indir '/ERP_LtH_HighvsMid.mat'], 'ERP_LtH_HighvsMid')
% LtH_High_low
[ERP_LtH_HighvsLow] = ft_timelockstatistics(cfg, LtH_High{:}, LtH_Low{:}); 
save ([indir '/ERP_LtH_HighvsLow.mat'], 'ERP_LtH_HighvsLow')
% LtH_Mid_Low
[ERP_LtH_MidvsLow] = ft_timelockstatistics(cfg, LtH_Mid{:}, LtH_Low{:}); 
save ([indir '/ERP_LtH_MidvsLow.mat'], 'ERP_LtH_MidvsLow')
%% PLOTS POST-HOC (ANOVA)
cd(indir)
load layout
load ERP_HtL_SelfvsHigh
load ERP_HtL_SelfvsMid
load ERP_HtL_SelfvsLow
load ERP_HtL_HighvsMid
load ERP_HtL_HighvsLow
load ERP_HtL_MidvsLow

load ERP_LtH_SelfvsHigh
load ERP_LtH_SelfvsMid
load ERP_LtH_SelfvsLow
load ERP_LtH_HighvsMid
load ERP_LtH_HighvsLow
load ERP_LtH_MidvsLow

cfg=[];
cfg.layout=layout;
cfg.highlightsizeseries = [15 15 15 15 15];
cfg.highlightsymbolseries=['.','.','.','.','.'];
cfg.zlim='maxabs';
cfg.alpha=0.05;

% HtL_Self_High
ft_clusterplot(cfg,ERP_HtL_SelfvsHigh);
% HtL_Self_Mid
ft_clusterplot(cfg,ERP_HtL_SelfvsMid);
% HtL_Self_Low
ft_clusterplot(cfg,ERP_HtL_SelfvsLow);
% HtL_High_Mid
ft_clusterplot(cfg,ERP_HtL_HighvsMid);
% HtL_High_Low
ft_clusterplot(cfg,ERP_HtL_HighvsLow);
% HtL_Mid_Low
ft_clusterplot(cfg,ERP_HtL_MidvsLow);

% LtH_Self_High
ft_clusterplot(cfg,ERP_LtH_SelfvsHigh);
% LtH_Self_Mid
ft_clusterplot(cfg,ERP_LtH_SelfvsMid);
% LtH_Self_Low
ft_clusterplot(cfg,ERP_LtH_SelfvsLow);
% LtH_High_Mid
ft_clusterplot(cfg,ERP_LtH_HighvsMid);
% LtH_High_Low
ft_clusterplot(cfg,ERP_LtH_HighvsLow);
% LtH_Mid_Low
ft_clusterplot(cfg,ERP_LtH_MidvsLow);

%% MAIN EFFECT OF INTERACTION (ANOVA)

for cond = {'meanallMA_selfMA_motorM' 'meanallMA_highMA_motorM' 'meanallMA_midMA_motorM' 'meanallMA_lowMA_motorM'}
        switch cond{:}
            case 'meanallMA_selfMA_motorM'
                i=1;
                for subj = 1:length(subarray)
                    load([indir num2str(subarray(subj),'%0.2d') '_HtL_' cond{:} '_' 'LtH_' cond{:} '_ERP.mat']);
                    Diff_Self{i}=raw_effect;
                    i = i+1
                end
            case 'meanallMA_highMA_motorM'
                i=1;
                for subj = 1:length(subarray)
                    load([indir num2str(subarray(subj),'%0.2d') '_HtL_' cond{:} '_' 'LtH_' cond{:} '_ERP.mat']);
                    Diff_High{i}=raw_effect;
                    i = i+1
                end
            case 'meanallMA_midMA_motorM'
                i=1;
                for subj = 1:length(subarray)
                    load([indir num2str(subarray(subj),'%0.2d') '_HtL_' cond{:} '_' 'LtH_' cond{:} '_ERP.mat']);
                    Diff_Mid{i}=raw_effect;
                    i = i+1
                end
            case 'meanallMA_lowMA_motorM'
                i=1;
                for subj = 1:length(subarray)
                    load([indir num2str(subarray(subj),'%0.2d') '_HtL_' cond{:} '_' 'LtH_' cond{:} '_ERP.mat']);
                    Diff_Low{i}=raw_effect;
                    i = i+1
                end
        end
     
  
        cd(indir)
        eval(['save ' 'Diff_Self Diff_Self']);
        eval(['save ' 'Diff_High Diff_High']);
        eval(['save ' 'Diff_Mid Diff_Mid']);
        eval(['save ' 'Diff_Low Diff_Low']);
 end
   
     %% Grand average for the main effects of INTERACTION (ANOVA)
     cd(indir)
     
     cfg=[]
    for cond =   {'Diff_Self' 'Diff_High' 'Diff_Mid' 'Diff_Low'}
        
        eval(['load ' cond{:}])
        
        eval(['tmp =' cond{:}])
        i=1
        
        for subj = 1:length(subarray) % NOT times two because they are NOT main effecs here
            tmp2 = tmp{subj};
            eval(['data_' num2str(i) ' = tmp2;']);
            i = i+1;
        end
        cfg.keepindividual = 'no';
        % In this following line you should write as many data_X as
        % subjects that you are including in the GA
        eval(['GA_' cond{:} '_ERP = ft_timelockgrandaverage(cfg,data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11, data_12, data_13, data_14, data_15, data_16,data_17,data_18, data_19, data_20)']);
        eval(['save(''' indir  'GA_' cond{:} '_ERP'',''GA_' cond{:} '_ERP'')']);
        eval(['clear(''GA_' cond{:} '_ERP'')']);
        
    end

%% % We perform the stats of Main Effects of INTERACTION (ANOVA)
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

load Diff_Self
load Diff_High
load Diff_Mid
load Diff_Low
% Cluster for One way anova 1T sequences
cfg.statistic='depsamplesFmultivariate';

% nº of levels of the other factor (1 [the difference]) * length(subarray)
subj = length(subarray);
% nº of factors of the independent variable (4 [Self High Mid Low]) * subj
design = zeros(2,4*subj);
% One per factor
for i = 1:subj
    design(1,i) = i;
end
for i = 1:subj
    design(1,subj+i)=i;
end
for i = 1:subj
    design(1,subj*2+i)=i;
end
for i = 1:subj
    design(1,subj*3+i)=i;
end

design(2,1:subj)=1; design(2,subj+1:2*subj)=2; design(2,2*subj+1:3*subj)=3; design(2,3*subj+1:4*subj)=4;
cfg.design = design;
cfg.ivar=2; cfg.uvar=1;

cfg.clusterthreshold='nonparametric';
ERP_oneway_interaction=ft_timelockstatistics(cfg, Diff_Self{:}, Diff_High{:}, Diff_Mid{:}, Diff_Low{:});
save ([indir 'ERP_oneway_interaction.mat'], 'ERP_oneway_interaction')

%% STATS PLOTS INTERACTION (ANOVA)
cd(indir)
load layout
load ERP_oneway_interaction
cfg=[];
cfg.layout=layout;
cfg.highlightsizeseries = [15 15 15 15 15];
cfg.highlightsymbolseries=['.','.','.','.','.'];
cfg.zlim='maxabs';
cfg.alpha=0.05;
ft_clusterplot(cfg,ERP_oneway_interaction);


    %% Difference waves
% cd(indir); 
% % Load the DW for the maincomps (already have them
% load GA_eA_eMA_eM_ERP; load GA_tAretr_tMAretr_ERP;
% 
% % For the interactions, load the GAvs of interest
% load GA_HtL_selfMA_LtH_selfMA_ERP
% load GA_HtL_highMA_LtH_highMA_ERP
% load GA_HtL_midMA_LtH_midMA_ERP
% load GA_LtH_lowMA_LtH_lowMA_ERP
% % ANd calculate DW of the DWs:
% % i.e., (tAencR-(tMAencR-eM)) - (tAencF-(tMAencF-eM))
% % (tAretrR-tMAretrR) - (tAretrF-tMAretrF)
% % (Same_A_Rem-Same_MA_Rem) - (Different_A_Rem-Different_MA_Rem)
% cfg = [];
% cfg.operation = 'subtract';
% cfg.parameter = 'avg';
% GA_DW_EncodingR_F = ft_math(cfg,GA_tAenc_R_2T_tMAenc_R_2T_eM_ERP, GA_tAenc_F_2T_tMAenc_F_2T_eM_ERP);
% GA_DW_RetrievalR_F = ft_math(cfg, GA_tAretr_R_2T_tMAretr_R_2T_ERP, GA_tAretr_F_2T_tMAretr_F_2T_ERP);
% GA_DW_Encoding_ME_sound = ft_math(cfg,GA_tAenc_R_2T_tAenc_F_2T_ERP, GA_tMAenc_R_2T_eM_tMAenc_F_2T_eM_ERP);
% 
% save ([indir '\GA_DW_EncodingR_F.mat'], 'GA_DW_EncodingR_F');
% save([indir '\GA_DW_RetrievalR_F.mat'], 'GA_DW_RetrievalR_F');
% 
% save([indir '\GA_DW_Encoding_ME_sound.mat'], 'GA_DW_Encoding_ME_sound');

%% Prepare the configuration file for cluster (same for all comps)
cd(indir)
load('all_dataSetsForClusters.mat');
load layout
load neighbours

    cd(scripts_folder)
    cfg = [];
    cfg.parameter        = 'avg';
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'depsamplesT';
    cfg.correctm         = 'cluster';
    cfg.clusteralpha     = 0.05;
    cfg.clusterstatistic = 'maxsum';
   
    cfg.minnbchan        = 2;  % should be set to 2 at least (Jordi's comment at meeting)
    cfg.tail             = 0; %-1, 1 or 0 one-sided or two-sided test
    cfg.clustertail      = 0;
    %cfg.alpha            = 0.025; %VERY IMPORTANT: If you have two-sided
    %test, set cfg.alpha to 0.025 if correcttail='alpha'
    cfg.numrandomization = 10000; % Pump up
    cfg.correcttail = 'prob';
    
    cfg.layout = layout;
    cfg.neighbours = neighbours;
    
    %cfg.design = [ ones(1,length(subarray)) ones(1,length(subarray))*2; 1:length(subarray) 1:length(subarray) ];
    %cfg.ivar     = 1; % Row in which the independent variable (group number) is set
    %cfg.uvar     = 23; % Row in which the unit of observation variable (subject number) is set
    
    %ncol = number of subjects * 2
    %nrow = 2
    subj = length(subarray);
    design = zeros(2,2*subj);
    %row1 = nº subjs * 2
    %*1
    for i = 1:subj
        design(1,i) = i;
    end
    %*2
    for i = 1:subj
        design(1,subj+i)=i;
    end
    %row2 = 1 * nº subjs; 2 * nº subjs
    design(2,1:subj)=1; design(2,subj+1:2*subj)=2;
    cfg.design = design;
    %ivar = independent variable; uvar = unit variable;
    cfg.ivar=2; cfg.uvar=1;
    %----------------------------------
    % Cluster for encoding: eA vs. eMA
    %----------------------------------
    [stats_eAvseMA_eM] = ft_timelockstatistics(cfg, C1{:}, C2{:});
    save ([indir '1_stats_eAvseMA_eM.mat'], 'stats_eAvseMA_eM')
    %----------------------------------
    % Cluster for retrieval: tA vs. tMA for 2T
    %----------------------------------
    [stats_tAvstMA_2T] = ft_timelockstatistics(cfg, C3{:}, C4{:});
    save ([indir '2a_stats_tAvstMA2T.mat'], 'stats_tAvstMA_2T')
 %----------------------------------
    % Cluster for retrieval: tA vs. tMA for 1T
    %----------------------------------
    [stats_tAvstMA_1T] = ft_timelockstatistics(cfg, tAretr1T{:}, tMAretr1T{:});
    save ([indir '2b_stats_tAvstMA.mat'], 'stats_tAvstMA_1T')
 
     %----------------------------------
    % Cluster for Interaction in enc and then retr: (AR-AF)-(MAR-MAF) as in FT
    % tutorial
    % a. T-test on the difference waves
    %----------------------------------
    [stats_encAvsencMA] = ft_timelockstatistics(cfg, eA_RvsF{:}, eMA_RvsF{:});
    save ([indir '3a_stats_interaction_encoding.mat'], 'stats_encAvsencMA')
       
    [stats_retrAvsretrMA_RvsF] = ft_timelockstatistics(cfg, tA_RvsF_2T{:}, tMA_RvsF_2T{:});
    save ([indir '4a_stats_interaction_retrieval.mat'], 'stats_retrAvsretrMA_RvsF')
    
    %'''''' MAIN EFFECTS FOR ANOVAS '''''
    
    % for main effects in "ANOVA " we need to change the designM
    subj_ME = length(subarray)*2;
    design_ME = zeros(2,2*subj_ME);
    for i = 1:subj_ME,  design_ME(1,i) = i;end
    
    for i = 1:subj_ME, design_ME(1,subj_ME+i)=i;end
    design_ME(2,1:subj_ME)=1; design_ME(2,subj_ME+1:2*subj_ME)=2;
    cfg.design = design_ME;
    
    % Main effects of sound and memory at encoding
    [stats_ME_sound_enc] =ft_timelockstatistics(cfg, tAenc{:}, tMAenc_eM{:})
    [stats_ME_memory_enc] = ft_timelockstatistics(cfg,remembered_enc{:}, forgotten_enc{:})
    save ([indir '3b_ME_sound_enc.mat'], 'stats_ME_sound_enc')
    save ([indir '3c_ME_memory_enc.mat'], 'stats_ME_memory_enc')

        % Main effects of sound and memory at retrieval
    [stats_ME_sound_retr] =ft_timelockstatistics(cfg, tAretr_ME_sound{:}, tMAretr_ME_sound{:})
    [stats_ME_memory_retr] = ft_timelockstatistics(cfg,remembered_retr{:}, forgotten_retr{:})
        save ([indir '4b_ME_sound_retr.mat'], 'stats_ME_sound_retr')
    save ([indir '4c_ME_memory_retr.mat'], 'stats_ME_memory_retr')
 
     %----------------------------------
    % Post hocs encoding and retrieval
    % - T-test for possible differences between Remembered MA vs Forgotten MA
    % - T-test for possible differences between Remembered A vs Forgotten A
    %----------------------------------
    subj = length(subarray);
    design = zeros(2,2*subj);
    for i = 1:subj
        design(1,i) = i;
    end
    for i = 1:subj
        design(1,subj+i)=i;
    end
    design(2,1:subj)=1; design(2,subj+1:2*subj)=2;
    cfg.design = design;
   
    [stats_encMA_RvsF] = ft_timelockstatistics(cfg, EncodingMA_R{:}, EncodingMA_F{:});    
    [stats_encA_RvsF] = ft_timelockstatistics(cfg, EncodingA_R{:}, EncodingA_F{:});

    [stats_retrMA_RvsF] = ft_timelockstatistics(cfg, RetrievalMA_R{:}, RetrievalMA_F{:});
 
    [stats_retrA_RvsF] = ft_timelockstatistics(cfg, RetrievalA_R{:}, RetrievalA_F{:});
    save ([indir '4b_stats_retrA_RvsF.mat'], 'stats_retrA_RvsF')
    %----------------------------------
    % Cluster for retrieval by order:
    % T-test on the difference waves
    %----------------------------------
   % [stats_SameDWvsDiffDW] = ft_timelockstatistics(cfg, Same{:}, Different{:});
   % save ([indir '5_stats_SameDWvsDiffDW.mat'], 'stats_SameDWvsDiffDW')
    % -------
%     % Cluster for One way anova 1T sequences
%     cfg.statistic='depsamplesFmultivariate';
%     design = zeros(2,3*length(subarray));nsub = length(subarray);
%     for i =1:nsub, design(1,i)=i; end
%     for i = 1:nsub, design(1,nsub+i)=i; end
%     for i =1:nsub, design(1,(nsub*2)+i)=i;end
%     design(2,1:nsub)=1;design(2,nsub+1:2*nsub)=2;
%     design(2,(2*nsub)+1:3*nsub)=3;
%     cfg.design = design;
%    % cfg.alpha=0.05;
%     %
%     cfg.clusterthreshold='nonparametric'; 
%     stats_oneway_1T=ft_timelockstatistics(cfg, tMAretr1T{:}, tAretr1T{:}, nosounds{:});
%     save ([indir '6_stats_oneway.mat'], 'stats_oneway_1T')

    
    %% PLOT 
% if not loaded, load them
cd(indir)
load('1_stats_eAvseMA_eM')
load('2a_stats_tAvstMA_2T')
load('2b_stats_tAvstMA_1T')
load('3a_stats_interaction_encoding')
load('3b_ME_sound_enc')
load('3c_ME_memory_enc')
load('4a_stats_interaction_retrieval')
load('4b_ME_sound_retr')
load('4c_ME_memory_retr')

%% Prepare the configuration file for cluster (same for all comps)
cd(indir)
load('all_dataSetsForClusters.mat');
load layout
load neighbours

% Load all saved analyses/stats/ERPs
cd(indir)
files = dir('*.mat');
tables = {};
for i = 1:length(files)
    load(files(i).name);
end
%t = significant clusters as shown by ft_clusterplot
t{1} = [.116 .234]; t{2}=[0.052  0.388];t{3}=[0.566  0.698];  
plotting2_avgWave(GA_eMA_eM_ERP, GA_eA_ERP, [],[],GA_eA_eMA_eM_ERP,stats_eAvseMA_eM, [], {'eMA','eA'},'Encoding',0.5,1,1, 0,t)
t=[];
plotting2_avgWave(GA_tMAretr_ERP, GA_tAretr_ERP, [],[],GA_tAretr_tMAretr_ERP,stats_tAvstMA_2T, [], {'tMA','tA'},'Retrieval 2T Sequences',0.5,1,1, 0,t)
t=[];
plotting2_avgWave(GA_tMAretr1T_ERP, GA_tAretr1T_ERP, [],[],GA_tAretr1T_tMAretr1T_ERP,stats_tAvstMA_1T, [], {'tMA','tA'},'Retrieval 1T Sequences',0.5,1,1, 0,t)

clear t; t{1} = [0.132 0.22]; t{2}=[0.096 0.31];
load('GA_tMAenc_R_2T_eM_ERP'); load('GA_tAenc_R_2T_ERP'); load('GA_tMAenc_F_2T_eM_ERP');  load('GA_tAenc_F_2T_ERP');
load('GA_DW_Encoding_ME_sound'); load('3b_ME_sound_enc.mat');
plotting2_avgWave(GA_tMAenc_R_2T_eM_ERP, GA_tMAenc_F_2T_eM_ERP,GA_tAenc_R_2T_ERP,GA_tAenc_F_2T_ERP,...
    GA_DW_Encoding_ME_sound,stats_ME_sound_enc, [], {'MA_R','MA_F','A_R','A_F'},'Encoding Interaction',0.5,1,1, 1,t)

clear t; t{1} = [0.426 0.698];
load('GA_tMAretr_R_2T_ERP'); load('GA_tMAretr_F_2T_ERP'); load('GA_tAretr_R_2T_ERP');  load('GA_tAretr_F_2T_ERP');
load('GA_DW_RetrievalR_F'); load('4c_ME_memory_retr.mat');

plotting2_avgWave(GA_tMAretr_R_2T_ERP, GA_tMAretr_F_2T_ERP,GA_tAretr_R_2T_ERP,GA_tAretr_F_2T_ERP,...
    GA_DW_RetrievalR_F,stats_ME_memory_retr, [], {'MA_R','MA_F','A_R','A_F'},'Retrieval Interaction',0.5,1,1, 1,t)
t=[];
plotting2_avgWave(GA_tMAretr1T_ERP, GA_tAretr1T_ERP, GA_NoSound_ERP,[],[],stats_oneway_1T, [], {'tMA','tA','New'},'Retrieval 1T',0.5,1,1, 2,t)


%% More plots
load layout;
load GA_HtL_ERP
load GA_LtH_ERP
load ERP_HtLvsLtH

cfg=[]
cfg.layout = layout;
cfg.parameter = 'stat';
ft_clusterplot(cfg, ERP_HtLvsLtH)
%ft_clusterplot(cfg, stats_tAvstMA_2T) % Not significant
%ft_clusterplot(cfg, stats_tAvstMA_1T) % Not significant

% ---Interactions
% ft_clusterplot(cfg, stats_encAvsencMA) %not significant
% ft_clusterplot(cfg, stats_retrAvsretrMA_RvsF) %not significant

%---Main effects at encoding
ft_clusterplot(cfg, stats_ME_sound_enc) %Significant
% ft_clusterplot(cfg, stats_ME_memory_enc) % not significant

%---Main effects at retrieval
%ft_clusterplot(cfg, stats_ME_sound_retr) % not significant
ft_clusterplot(cfg, stats_ME_memory_retr) %significant

% Plot ERP averages across all electrodes for the significant cluster-based
% comparisons first
% encoding main

% --- for plotting with variance: boundedline 
% eA = squeeze(mean(GA_eA_ERP.avg,1)) ; boundseA = std(eA, [], 2);
% eMA = squeeze(mean(GA_eMA_eM_ERP.avg,1)) ; boundseMA = std(eMA, [], 2);
% figure();
% boundedline(GA_eA_ERP.time, eA, boundseA, 'alpha', 'b','transparency', 0.1); hold on;
% boundedline(GA_eA_ERP.time, eMA, boundseMA, 'alpha', 'r','transparency', 0.1); hold on;
% xlim([GA_eA_ERP.time(1)  GA_eA_ERP.time(end)]); hold on;
% set(gca,'Ydir','reverse')
% ----------------------
cfg = []
cfg.channel = 'all';
cfg.linewidth = 1;
cfg.fontsize = 10
figureMainComps = figure;

subplot(222)
cfg.title = 'Main effect of Session'
ft_singleplotER(cfg, GA_HtL_ERP, GA_LtH_ERP);
set(gca,'Ydir','reverse')
legend({'HtL', 'LtH'});
xlabel('Time'); 
ylabel('mV');

subplot(121)
cfg.title = 'Encoding'
ft_singleplotER(cfg, GA_eA_ERP, GA_eMA_eM_ERP)
set(gca,'Ydir','reverse')
legend({'Passive', 'Self-generated'}, 'Location', 'NorthEast')
xlabel('Time'); 
ylabel('mV');
subplot(122)
cfg.title = 'Retrieval'
ft_singleplotER(cfg, GA_tAretr_ERP, GA_tMAretr_ERP)
set(gca,'Ydir','reverse')
legend({'Encoded as Passive', 'Encoded as Self-generated'}, 'Location', 'NorthEast')
xlabel('Time'); 
ylabel('mV');

figureMainEffects_ENC_RETR = figure;
subplot(221)
cfg.title = 'Main effect of sound at encoding'
ft_singleplotER(cfg, GA_tAenc_ERP, GA_tMAenc_eM_ERP)
set(gca,'Ydir','reverse')
legend({'test as A', 'test as MA'});
xlabel('Time'); 
ylabel('mV');

subplot(222)
cfg.title = 'Main effect of Memory at encoding'
ft_singleplotER(cfg, GA_remembered_enc_ERP, GA_forgotten_enc_ERP);
set(gca,'Ydir','reverse')
legend({'Later Remembered', 'Later Forgotten'});
xlabel('Time'); 
ylabel('mV');

subplot(223)
cfg.title = 'Main effect of sound at retrieval'
ft_singleplotER(cfg, GA_tAretr_ME_sound_ERP, GA_tMAretr_ME_sound_ERP)
set(gca,'Ydir','reverse')
legend({'Encoded as Passive', 'Encoded as SG'});
xlabel('Time'); 
ylabel('mV');


subplot(224)
cfg.title = 'Main effect of memory at retrieval'
ft_singleplotER(cfg, GA_remembered_retr_ERP, GA_forgotten_retr_ERP)
set(gca,'Ydir','reverse')
legend({'Remembered', 'Forgotten'});
xlabel('Time'); 
ylabel('mV'); 
    
    plot_clusters(stats_eAvseMA_eM, GA_eA_eMA_eM_ERP)
      plot_clusters(stats_encAvsencMA, GA_)

%% OLD CODE FOR STATS PLOTS - Encoding first
% cd(indir)
% load layout
% cfg=[];
% cfg.layout=layout;
% 
% % cfg.highlightsymbolseries=['*','*','*','*','*'];
% cfg.zlim='maxabs';
% cfg.parameter='stat'
% %cfg.maskparameter = 'mask';
% cfg.highlight = 'on';
% cfg.highlightcolorpos = [0 0 0];
% cfg.hightlightcolorneg = [1 1 1];
% ft_clusterplot(cfg,stats_eAvseMA_eM);
% % ''' with cfg.alpha not specified
% % There are 3 clusters smaller than alpha (0.05)
% % Positive cluster: 1, pvalue: 0.011599 (x), t = 0.116 to 0.234
% % Negative cluster: 1, pvalue: 0.0009999 (*), t = 0.052 to 0.388
% % Negative cluster: 2, pvalue: 0.037396 (x), t = 0.566 to 0.698
% 
% % +++  with cfg alpha set to 0.025
% % There are 3 clusters smaller than alpha (0.05)
% % Positive cluster: 1, pvalue: 0.013799 (x), t = 0.116 to 0.234
% % Negative cluster: 1, pvalue: 0.00059994 (*), t = 0.052 to 0.388
% % Negative cluster: 2, pvalue: 0.036596 (x), t = 0.566 to 0.698
% ft_clusterplot(cfg, stats_tAvstMA);
% %% From these stat plots we know now what time windows we want to plot. in this case, from 650ms to 1000ms in steps of 50ms.
% %AGExCOND_T2BvsT1B_YvsO_2sec_18062019.maskedstat=AGExCOND_T2BvsT1B_YvsO_2sec_18062019.stat .* AGExCOND_T2BvsT1B_YvsO_2sec_18062019.mask;
% % figure;
% % cfg=[];
% % cfg.layout=layout;
% % %cfg.xlim=[.200 .400];
% % cfg.colorbar='WestOutside';
% % cfg.parameter='stat';
% % cfg.highlightcolorpos = [1 1 1];
% % cfg.highlightcolorneg = [0.5 0.5 0.5];
% % 
% % cfg.zlim= 'maxabs';%'absmax';
% % cfg.latency =[.05 .150];
% % %cfg.maskparameter = 'mask';
% % cfg.highlight ='on';
% % cfg.opacitymap='rampup';
% % cfg.gridscale=300;
% % cfg.colormap = 'jet';
% % cfg.shading = 'interp';
% % cfg.style='fill';
% % cfg.interplimits='electrodes';
% % cfg.interpolatenan='no';
% % ft_topoplotER(cfg,ERP_eAvseMA_eM_P2)
% 
% %% PLOTTING
% cd(indir); 
% % Load the DW for the maincomps (already have them
% load GA_eA_eMA_eM_ERP; load GA_tAretr_tMAretr_ERP;
% 
% % For the interactions, load the GAvs of interest
% load GA_tAenc_R_2T_tMAenc_R_2T_eM_ERP
% load GA_tAenc_F_2T_tMAenc_F_2T_eM_ERP
% load GA_tAretr_R_2T_tMAretr_R_2T_ERP
% load GA_tAretr_F_2T_tMAretr_F_2T_ERP
% load GA_Same_A_Rem_Same_MA_Rem_ERP
% load GA_Different_A_Rem_Different_MA_Rem_ERP
% % ANd calculate DW of the DWs:
% % i.e., (tAencR-(tMAencR-eM)) - (tAencF-(tMAencF-eM))
% % (tAretrR-tMAretrR) - (tAretrF-tMAretrF)
% % (Same_A_Rem-Same_MA_Rem) - (Different_A_Rem-Different_MA_Rem)
% cfg = [];
% cfg.operation = 'subtract';
% cfg.parameter = 'avg';
% GA_DW_EncodingR_F = ft_math(cfg,GA_tAenc_R_2T_tMAenc_R_2T_eM_ERP, GA_tAenc_F_2T_tMAenc_F_2T_eM_ERP);
% GA_DW_RetrievalR_F = ft_math(cfg, GA_tAretr_R_2T_tMAretr_R_2T_ERP, GA_tAretr_F_2T_tMAretr_F_2T_ERP);
% GA_DW_Same_Different = ft_math(cfg, GA_Same_A_Rem_Same_MA_Rem_ERP, GA_Different_A_Rem_Different_MA_Rem_ERP);
% save ([indir 'GA_DW_EncodingR_F.mat'], 'GA_DW_EncodingR_F');
% save([indir 'GA_DW_RetrievalR_F.mat'], 'GA_DW_RetrievalR_F');
% save ([indir 'GA_DW_Same_Different.mat'], 'GA_DW_Same_Different');
% 
% % ---- USE SUBFUNCTION FOR PLOTTING OTHERWISE IT'S A MESS

% %% Create plots for retrieval
% % Cluster plot first
% cd(indir)
% load layout
% cfg=[];
% cfg.layout=layout;
% cfg.latency = [0.05 0.3]
% cfg.highlightcolorpos = [1 1 1];
% cfg.highlightcolorneg = [0.5 0.5 0.5];
% % cfg.highlightsymbolseries=['*','*','*','*','*'];
% cfg.zlim='maxabs';
% cfg.parameter='stat'
% cfg.maskparameter = 'mask';
% cfg.alpha=0.05;
% %cfg.colorbar ='yes';
% ft_clusterplot(cfg,ERP_tAvstMA);
% % then topoplot
% cd(indir); load GA_tAretr_ERP; load GA_tMAretr_ERP;
% 
% cfg = [];
% cfg.operation = 'subtract';
% cfg.parameter = 'avg';
% GA_tAretr_vs_tMAretr = ft_math(cfg,GA_tAretr_ERP, GA_tMAretr_ERP)
% 
% 
% 
% 
% 
% %----
% figure;
% % define params for plotting 
% timestep = 0.05 %in secs
% sampling_rate = 500; %
% sample_count = length(ERP_tAvstMA.time);
% % j and M must have the same temporal length
% j = [0:timestep:ERP_tAvstMA.time(end)]; %temporal endpoints (in secs) of the ERP avg computed in each subplot
% m = [1:timestep*sampling_rate:sample_count]; %temporal endpoints in EEG samples
% % get relevant (significant) values for positive
% pos_cluster_pvals = [ERP_tAvstMA.posclusters(:).prob];
% pos_signif_clust = find(pos_cluster_pvals < ERP_tAvstMA.cfg.alpha);
% pos = ismember(ERP_tAvstMA.posclusterslabelmat, pos_signif_clust);
% 
% % get relevant (significant) values for negative (if any)
% if ~isempty(ERP_tAvstMA.negclusters)
%     neg_cluster_pvals = [ERP_tAvstMA.negclusters(:).prob];
%     neg_signif_clust = find(neg_cluster_pvals < ERP_tAvstMA.cfg.alpha);
%     neg = ismember(ERP_tAvstMA.negclusterslabelmat, neg_signif_clust);
% 
% end
% % First ensure the channels to have the same order in the avg and in the
% % stat output. This might not be the case, because ft_math might shuffle
% % the order
% [i1,i2] = match_str(GA_tAretr_vs_tMAretr.label, ERP_tAvstMA.label);
%     
% % plot pliz for positive/negative clusters
% for k = 1:length(m)-1
%     subplot(2,5,k);
%     cfg=[];
%     cfg.xlim = [j(k) j(k+1)];
%    % cfg.zlim = [-5e-14 5e-14];
%     pos_int = zeros(numel(GA_tAretr_vs_tMAretr.label),1);
%     pos_int(i1) = all(pos(i2, m(k):m(k+1)),2);
% %     neg_int = zeros(numel(GA_eA_vs_eMA_eM.label),1);
% %     neg_int(i1) = all(neg(i2, m(k):m(k+1)),2);
%     cfg.highlight = 'on';
%     %cfg.highlightchannel = find(pos_int);
%     cfg.comment = 'xlim';
%     cfg.maskparameter = 'mask';
%     cfg.commentpos = 'title';
%     %cfg.interactive ='no';
%     cfg.layout = layout;
%     ft_topoplotER(cfg,GA_tAretr_vs_tMAretr);
% end
% %% try other plots
% cfg = [];
% cfg.layout= 'ordered';
% layout2=ft_prepare_layout(cfg,GA_eA_ERP)
% cfg=[];
% cfg.channel = 'all' % Cz is 28
% ft_singleplotER (cfg, GA_eA_ERP, GA_eMA_eM_ERP)
% 
% chan = 28;
% time = [0 400];
% %Scale the vertical ax
% 
% figure;
% for iSub = 1:numel(subarray)
%     subplot(4,6,iSub)
%     plot(C1{iSub}.time, C1{iSub}.avg(chan,:), 'b'); %eA
%     hold on;
%     plot(C2{iSub}.time, C2{iSub}.avg(chan,:), 'r') %eMA
%     title(strcat('Sub: ', num2str(subarray(iSub))));
%     hold on;
%   
%     xlim([C1{iSub}.time(1) C1{iSub}.time(end)])
% end
% %legend({'eA','eMA'}, 'Location', 'SouthEast');  
% 
% % find data points for the effect of interest in the GA
% chan = 28;
% time = [0.100 0.180]; %N1 timewindow
% timesel=find(GA_eA_ERP.time >= time(1) & GA_eA_ERP.time <= time(2));
% 
% for iSub = 1:numel(subarray)
%     values_eA(iSub) = mean(C1{iSub}.avg(chan,timesel));
%     values_eMA(iSub) = mean(C2{iSub}.avg(chan,timesel));
%     sub{iSub} = ['s. ' num2str(subarray(iSub))];
% end
% 
% M = [values_eA', values_eMA'];
% figure;
% plot(M','o-'); xticks([1 2]); xticklabels({'eA', 'eMA'});xlabel(['Condition']);
% legend(sub, 'location', 'eastoutside')
 end