function step_2_2_GA_cluster
% We use script clusterStatistics_nadia.
% Creates data for cluster based and performs the following stats:
% a) Permutation paired t-test for differences between eA vs. eMA
% b) Permutation paired t-test for differences between tAretr vs. tMAretr
% c) Somehow we need to do the anova for Memory by Generation

% For interactions effects see the following:
% https://mailman.science.ru.nl/pipermail/fieldtrip/2011-January/003447.html
% https://mailman.science.ru.nl/pipermail/fieldtrip/2010-December/003338.html

%% Paths
clear; clc; close all;
analysis_path = '/DATA2/BLB_EXP_201902_SGMem2/BLB_BackUp_files/Protocols/Analysis_Scripts/Pupillometry';
data_path = '/DATA2/BLB_EXP_201902_SGMem2/Analysis/EyeTracking/first-level';
raw_mat = '/DATA2/BLB_EXP_201902_SGMem2/Analysis/EyeTracking';
raw_asc=raw_mat;
% add plugins needed
addpath([analysis_path '/Tools-master_AU/eye']);
addpath([analysis_path '/Tools_Nadia/boundedline-pkg-master/boundedline-pkg-master/boundedline']);
addpath([analysis_path '/Tools_Nadia/boundedline-pkg-master/boundedline-pkg-master/catuneven']);
addpath([analysis_path '/Tools_Nadia/boundedline-pkg-master/boundedline-pkg-master/Inpaint_nans']);
addpath([analysis_path '/Tools_Nadia/boundedline-pkg-master/boundedline-pkg-master/singlepatch']);

addpath([analysis_path '/Tools-master_AU/plotting']);
addpath([analysis_path '/Tools-master_AU/plotting/cbrewer'])
analysis_log = [analysis_path '/analysis_log.txt'];

all_subs_folder = ['/DATA2/BLB_EXP_201902_SGMem2/Analysis/EyeTracking/all_subjects/'];
ft_defaults;
% File
nBlocks=24;
subjects = [6 7 8 9 10 11 12 13 14 16 18 19 20 21 22 23 24 25 26 27 28 29 30];
tic;
fileID = fopen(analysis_log,'wt');
fprintf(fileID,'Data path: %s*\n\n',data_path);
fprintf(fileID, '*********************\n\n');
fprintf(fileID,'Creating Grand averages per condition...');

preparation = 0;

if preparation == 0 % if the data is not prepared
    % Create empty cell arrays for the conditions of interest
    for cond = {'eA_2T' 'eMA_2T' 'eM_2T' 'eA' 'eMA' 'eM' 'tAretr' 'tMAretr' 'tAretr1T' 'tMAretr1T' 'nosounds' 'tAenc_R_2' 'tAenc_F_2' 'tMAenc_R_2' 'tMAenc_F_2' 'tAretr_R_2' 'tAretr_F_2' 'tMAretr_R_2' 'tMAretr_F_2' 'response' 'retention' 'trialstart' 'blockstart'}
        eval( [cond{:}  ' = cell(length(subjects),1)'])
    end
    %% Prepare data for cluster based stats
    % if already done don't continue
    cd(all_subs_folder)
    
    
    for cond = {'eA_2T' 'eMA_2T' 'eM_2T' 'eA' 'eMA' 'eM' 'tAretr' 'tMAretr' 'tAretr1T' 'tMAretr1T' 'nosounds' 'tAenc_R_2' 'tAenc_F_2' 'tMAenc_R_2' 'tMAenc_F_2' 'tAretr_R_2' 'tAretr_F_2' 'tMAretr_R_2' 'tMAretr_F_2' 'response' 'retention' 'trialstart' 'blockstart'}
        
        % Don't continue if already done
        cd(all_subs_folder)
        i=1;
        for subj = 1:length(subjects)
            evoked_folder = [data_path '/evoked_' num2str(subjects(subj))];
            
            load([evoked_folder '/' num2str(subjects(subj),'%0.2d') '_AVG_' cond{:} '_forClust']);
           
            eval([[['all_' cond{:}] '{' num2str(i) '}'] '=avg']')
            
            i = i+1;
        end
        eval(['save(''' all_subs_folder  'GA_' cond{:} '_all'', ''all_' cond{:} ''')']);
        
    end
    
    %% Get the eA-eMA differnece
      condArray1 = {'eA_2T'}
    condArray2 = {'eMA_2T'}
    for iDiff = 1:length(condArray1)
        for subj = 1:length(subjects)
            
            eval(['cond1 = ' [['all_' condArray1{iDiff}] '{' num2str(subj) '}']]');
            eval(['cond2 = ' [['all_' condArray2{iDiff}] '{' num2str(subj) '}']]');
        
            cfg = [];
            cfg.operation = 'subtract';
            cfg.parameter = 'avg';
            eval([['all_' condArray1{iDiff} '_' condArray2{iDiff}] '{' num2str(subj) '}' '= ft_math(cfg, cond1, cond2)']');
        end
        eval(['save(''' all_subs_folder  'GA_' condArray1{iDiff} '_' condArray2{iDiff} '_all'', ''all_' condArray1{iDiff} '_' condArray2{iDiff} ''')']);
    end
    
    % Now for ME at encoding and at retrieval
    for main_effect = {'tAenc_2ME', 'tMAenc_2ME', 'tAretr_2ME', 'tMAretr_2ME'  'enc_R_2ME', 'enc_F_2ME', 'retr_R_2ME', 'retr_F_2ME'}
        switch main_effect{:}
            case 'tAenc_2ME'
                eval(['all_' main_effect{:} ' = [all_tAenc_R_2, all_tAenc_F_2]']);
            case 'tMAenc_2ME'
                eval(['all_' main_effect{:} ' = [all_tMAenc_R_2, all_tMAenc_F_2]']);
            case 'tAretr_2ME'
                eval(['all_' main_effect{:} ' = [all_tAretr_R_2, all_tAretr_F_2]']);
            case 'tMAretr_2ME'
                eval(['all_' main_effect{:} ' =  [all_tMAretr_R_2, all_tMAretr_F_2]']);
            case  'enc_R_2ME'
                eval(['all_' main_effect{:} ' = [all_tAenc_R_2, all_tMAenc_R_2]']);
            case 'enc_F_2ME'
                eval(['all_' main_effect{:} ' = [all_tAenc_F_2, all_tMAenc_F_2]']);
            case 'retr_R_2ME'
                eval(['all_' main_effect{:} ' = [all_tAretr_R_2, all_tMAretr_R_2]']);
            case 'retr_F_2ME'
                eval(['all_' main_effect{:} ' = [all_tAretr_F_2, all_tMAretr_F_2]']);
        end
        eval(['save(''' all_subs_folder  'GA_' main_effect{:} '_all'', ''all_' main_effect{:} ''')']);
        
    end
    %% Test for interaction effects as in fieldtrip tutorial:
    % We have a factorial design:
    % The four cells in this design are denoted by (1,1), (1,2), (2,1) and (2,2)
    % First number indicates sound type (1 = A, 2= MA)
    % Second number indicates memory (1 = Remembered, 2= Forgotten)
    % So, for example (1,1) is for A-Remembered
    
    % I have to create to structures for the following comparisons:
    % GAdiff11_12 to test A-Remembered vs. A-Forgotten
    % GAdiff21_22 to test MA-Remembered vs. MA-Forgotten
    
    % Comparing GAdiff11_12 and GAdiff21_22 would be testing an interaction effect.
    
    condArray1 = { 'tAenc_R_2',   'tMAenc_R_2', 'tAretr_R_2', 'tMAretr_R_2'}
    condArray2 = { 'tAenc_F_2',   'tMAenc_F_2', 'tAretr_F_2', 'tMAretr_F_2'}
    for iDiff = 1:length(condArray1)
        for subj = 1:length(subjects)
            
            eval(['cond1 = ' [['all_' condArray1{iDiff}] '{' num2str(subj) '}']]');
            eval(['cond2 = ' [['all_' condArray2{iDiff}] '{' num2str(subj) '}']]');
        
            cfg = [];
            cfg.operation = 'subtract';
            cfg.parameter = 'avg';
            eval([['all_' condArray1{iDiff} '_' condArray2{iDiff}] '{' num2str(subj) '}' '= ft_math(cfg, cond1, cond2)']');
        end
        eval(['save(''' all_subs_folder  'GA_' condArray1{iDiff} '_' condArray2{iDiff} '_all'', ''all_' condArray1{iDiff} '_' condArray2{iDiff} ''')']);
    end
else
    %% If the data is already prepared just load it
    cd(all_subs_folder)
    files = dir('*.mat');
    tables = {};
    for i = 1:length(files)
        load(files(i).name);
    end
    %% GRAND AVERAGE ACROSS SUBJECTS
    cfg =[]
    cd(all_subs_folder)
    
    
    conds = {'eA_2T' 'eMA_2T' 'eM_2T' 'eA' 'eMA' 'eM' 'tAretr' 'tMAretr' 'tAretr1T' 'tMAretr1T' 'nosounds' 'tAenc_R_2' 'tAenc_F_2' 'tMAenc_R_2' 'tMAenc_F_2' 'tAretr_R_2' 'tAretr_F_2' 'tMAretr_R_2' 'tMAretr_F_2' 'tAenc_2ME', 'tMAenc_2ME', 'tAretr_2ME', 'tMAretr_2ME'  'enc_R_2ME', 'enc_F_2ME', 'retr_R_2ME', 'retr_F_2ME'}
    tAenc_R_2_tAenc_F_2' 'tMAenc_R_2_tMAenc_F_2' 'tAretr_R_2_tAretr_F_2' 'tMAretr_R_2_tMAretr_F_2' 'response' 'retention' 'trialstart' 'blockstart'}
    for cond = conds
        % Just load if already done
        %     if exist(sprintf(['GA_' cond{:} '_EPR.mat']), 'file'),
        %         load(['GA_' cond{:} '_EPR']);
        %         continue;
        %     else        
        disp('Rerunning...Grand average')
        i=1;
        for subj = 1:length(subjects)
            %evoked_folder = [data_path '\' num2str(subjects(subj)) '\evoked_' num2str(subjects(subj))];
            load([all_subs_folder 'GA_' cond{:} '_all'])
            eval([' tmp  = ' [[['all_' cond{:}] '{' num2str(subj) '}']]]')
            eval(['data_' num2str(i) ' = tmp']);
            i = i+1;
        end
        cfg.keepindividual = 'no';
        cfg.parameter = 'avg';
        % In this following line you should write as many data_X as
        % subjects that you are including in the GA
        eval(['GA_' cond{:} '_EPR = ft_timelockgrandaverage(cfg,data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11, data_12, data_13, data_14, data_15, data_16,data_17,data_18, data_19, data_20, data_21, data_22, data_23)']);
        eval(['save(''' all_subs_folder  'GA_' cond{:} '_EPR'',''GA_' cond{:} '_EPR'')']);       
        eval(['clear(''GA_' cond{:} '_EPR'')']);
    end
    
    % Grand average for difference wave
    cond= {'eA_2T_eMA_2T' }
    load([all_subs_folder 'GA_' cond{:} '_all'])
    
    disp('Rerunning...Grand average')
    i=1;
    for subj = 1:length(subjects)
        %evoked_folder = [data_path '\' num2str(subjects(subj)) '\evoked_' num2str(subjects(subj))];
        eval([' tmp  = ' [[['all_' cond{:}] '{' num2str(subj) '}']]]')
        eval(['data_' num2str(i) ' = tmp']);
        i = i+1;
    end
    cfg.keepindividual = 'no';
    cfg.parameter = 'avg';
    % In this following line you should write as many data_X as
    % subjects that you are including in the GA
    eval(['GA_' cond{:} '_EPR = ft_timelockgrandaverage(cfg,data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11, data_12, data_13, data_14, data_15, data_16,data_17,data_18, data_19, data_20, data_21, data_22, data_23)']);
    eval(['save(''' all_subs_folder  'GA_' cond{:} '_EPR'',''GA_' cond{:} '_EPR'')']);
    eval(['clear(''GA_' cond{:} '_EPR'')']);
    
end


%% Non parametric t-based permutation
nsubj = 23;
% --------------------------------------------------------------------------
% MAIN COMPARISONS : ENCODING
[stat_eA_vs_eMA_2T] = clusterStatistics_nadia(all_eA_2T, all_eMA_2T, nsubj)
[stat_eA_vs_eMA_all] = clusterStatistics_nadia(all_eA, all_eMA, nsubj);
[stat_tA_vs_tMA_2T] = clusterStatistics_nadia(all_tAretr, all_tMAretr, nsubj);
[stat_tA_vs_tMA_1T] = clusterStatistics_nadia(all_tAretr1T, all_tMAretr1T, nsubj);
save([all_subs_folder 'stat_eA_vs_eMA_2T'],'stat_eA_vs_eMA_2T');

save([all_subs_folder 'stat_eA_vs_eMA_all'],'stat_eA_vs_eMA_all');
save([all_subs_folder 'stat_tA_vs_tMA_2T'],'stat_tA_vs_tMA_2T');
save([all_subs_folder 'stat_tA_vs_tMA_1T'],'stat_tA_vs_tMA_1T');

% INTERACTIONS
[stat_interaction_encoding] = clusterStatistics_nadia(all_tAenc_R_2_tAenc_F_2, all_tMAenc_R_2_tMAenc_F_2, nsubj);
[stat_interaction_retrieval] = clusterStatistics_nadia(all_tAretr_R_2_tAretr_F_2, all_tMAretr_R_2_tMAretr_F_2, nsubj);
save([all_subs_folder 'stat_interaction_encoding'],'stat_interaction_encoding');
save([all_subs_folder 'stat_interaction_retrieval'],'stat_interaction_retrieval');

% MAIN EFFECTS AT ENCODING
nsubj = length(subjects)*2
stat_ME_memory_at_enc = clusterStatistics_nadia(all_enc_R_2ME, all_enc_F_2ME, nsubj);
stat_ME_sound_at_enc = clusterStatistics_nadia(all_tAenc_2ME, all_tMAenc_2ME, nsubj);
% MAIN EFFECTS AT retrieval
stat_ME_memory_at_retr = clusterStatistics_nadia(all_retr_R_2ME, all_retr_F_2ME, nsubj);
stat_ME_sound_at_retr = clusterStatistics_nadia(all_tAretr_2ME, all_tMAretr_2ME, nsubj);

save([all_subs_folder 'stat_ME_memory_at_enc'],'stat_ME_memory_at_enc');
save([all_subs_folder 'stat_ME_sound_at_enc'],'stat_ME_sound_at_enc');
save([all_subs_folder 'stat_ME_memory_at_retr'],'stat_ME_memory_at_retr');
save([all_subs_folder 'stat_ME_sound_at_retr'],'stat_ME_sound_at_retr');

%% ===========================================
%% One way anova for 1T sequences
%% ===========================================

cfg                  = [];
cfg.method           = 'montecarlo'; % permutation test
cfg.statistic        = 'depsamplesFmultivariate'; %
% do cluster correction
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; % threshold to consider smth a cluster. If cfg.clusteralpha is equal to 0.05, the t-values are
% thresholded at the 95-th quantile for a one-sided t-test, and at the 2.5-th and the 97.5-th quantiles for a two-sided t-test.
% cfgstats.clusterstatistic = 'maxsize'; % weighted cluster mass needs cfg.wcm_weight...
% cfgstats.minnbchan        = 1; % average over chans
cfg.tail             = 0; % two-tailed!
cfg.clustertail      = 0; % two-tailed!
cfg.alpha            = 0.025;
cfg.clusterthreshold = 'nonparametric'; %added by Nad
%cfg.correcttail = 'prob'; %added by Nad
cfg.numrandomization = 10000; % make sure this is large enough
cfg.randomseed       = 1; % make the stats reproducible!

% use only our preselected sensors for the time being
cfg.channel          = 'EyePupil';
nsubj=23;
design = zeros(2,3*nsubj);
for i = 1:nsubj,  design(1,i) = i;        end
for i = 1:nsubj,  design(1,nsubj+i) = i;  end
for i = 1:nsubj,  design(1,(nsubj*2)+i) = i;  end

design(2,1:nsubj)         = 1;
design(2,nsubj+1:2*nsubj) = 2;
design(2,(2*nsubj)+1:3*nsubj) = 3;

cfg.design = design;             % design matrix
cfg.uvar     = 1;               % Row in which the unit of observation variable (subject number) is set
cfg.ivar  = 2;                   % Row in which the independent variable (group number) is set
% specifies with which sensors other sensors can form clusters

stat_oneway_1T = ft_timelockstatistics(cfg, all_tAretr1T{:}, all_tMAretr1T{:}, all_nosounds{:})
save([all_subs_folder 'stat_oneway_1T'],'stat_oneway_1T');

%% ===========================================
cd(all_subs_folder)
files = dir('*EPR.mat');
tables = {};
for i = 1:length(files)
    load(files(i).name);
end
cd(analysis_path)
enc= figure;
legnames = {'Motor-auditory', 'Auditory-only', 'Difference'}
titl = 'Encoding'
% Find where was the significant effect
time = linspace(min(GA_eMA_2T_EPR.time), max(GA_eMA_2T_EPR.time), numel(GA_eMA_2T_EPR.time));
significant_effect = time(find(stat_eA_vs_eMA_2T.mask==1));
significant_time_window = [significant_effect(1) significant_effect(end)]
%subplot(121)
step_2_2_a_plotsEvoked(GA_eMA_2T_EPR, GA_eA_2T_EPR, GA_eA_2T_eMA_2T_EPR, 0,stat_eA_vs_eMA_2T,0,legnames,titl,0.1,1,1,0)
cd(all_subs_folder)
saveas(enc,'encoding_2T.fig');

% enc1= figure;
% legnames = {'Motor-auditory', 'Auditory-only'}
% titl = 'Encoding'
% %subplot(121)
% step_2_2_a_plotsEvoked(GA_eMA_EPR, GA_eA_EPR, 0, 0,stat_eA_vs_eMA_all,0,legnames,titl,0.1,1,1,0)
% cd(all_subs_folder)
% saveas(enc1,'encoding_all.fig');

% retr2T=figure;
% titl2 = 'Retrieval 2T'
% legnames = {'Encoded as Self-generated', 'Encoded as Externally-generated'}
% step_2_2_a_plotsEvoked(GA_tMAretr_EPR , GA_tAretr_EPR, 0, 0 ,stat_tA_vs_tMA_2T, 0,legnames,titl2,0.1,1,1,0)
% saveas(retr2T,'retr2T.fig');
% 
% retr1T=figure;
% titl2a = 'Retrieval 1T'
% legnames = {'Encoded as Self-generated', 'Encoded as Externally-generated'}
% step_2_2_a_plotsEvoked(GA_tMAretr1T_EPR , GA_tAretr1T_EPR, 0, 0 ,stat_tA_vs_tMA_1T, 0,legnames,titl2a,0.1,1,1,0)
% saveas(retr1T,'retr1T.fig');

cd(analysis_path)
% Main effect encoding & interaction
figure2 = figure;
% % titl2 = {'Sound type at encoding as a function of later memory performance'}
% % legnames = {'EG_F_o_r_g_o_t_t_e_n','EG_R_e_m_e_m_b_e_r_e_d','SG_F_o_r_g_o_t_t_e_n','SG_R_e_m_e_m_b_e_r_e_d'}
% % subplot(121)
% % step_2_2_a_plotsEvoked(  GA_tAenc_F_2_EPR ,GA_tAenc_R_2_EPR, GA_tMAenc_F_2_EPR, GA_tMAenc_R_2_EPR ,stat_ME_sound_at_enc, stat_ME_memory_at_enc,legnames,titl2,0.5,1,1,1)
% % subplot(122)
% Find where the significant main effect of memory is
time = linspace(min(GA_tAretr_F_2_EPR.time), max(GA_tAretr_F_2_EPR.time), numel(GA_tAretr_F_2_EPR.time));
significant_effect_retr = time(find(stat_ME_memory_at_retr.mask==1));
significant_time_window_retr = [significant_effect_retr(1) significant_effect_retr(end)]

titl2 = {'Sound type at retrieval as a function of memory performance'}

legnames = {'Encoded as A and remembered','Encoded as A and forgotten','Encoded as MA and remembered','Encoded as MA and forgotten'}
step_2_2_a_plotsEvoked( GA_tAretr_R_2_EPR, GA_tAretr_F_2_EPR , GA_tMAretr_R_2_EPR, GA_tMAretr_F_2_EPR ,stat_ME_sound_at_retr, stat_ME_memory_at_retr,legnames,titl2,0,1,1,1)
cd(all_subs_folder)
saveas(figure2,'figure2.fig');


%%------
figure3 = figure;

titl2 = {'Sound type at retrieval for 1T sequences'}
legnames = {'Encoded as self-generated','Encoded as externally-generated', 'New sounds'}
step_2_2_a_plotsEvoked( GA_tMAretr1T_EPR , GA_tAretr1T_EPR, GA_nosounds_EPR ,0, stat_oneway_1T, 0,legnames,titl2,0.2,1,1,2)
saveas(figure3,'figure3.fig');

%% Extract peak for each participant and condition in the window 0 to 1s post stimulus. Take latency as well
matrix_peaks={};
conds = {'eA_2T_eMA_2T', 'eA_2T' 'eMA_2T' 'eM_2T'}
    %'eA' 'eMA' 'eM' 'tAretr' 'tMAretr' 'tAretr1T' 'tMAretr1T' 'nosounds' ...
   % 'tAenc_R_2' 'tAenc_F_2' 'tMAenc_R_2' 'tMAenc_F_2' ...
    %'tAretr_R_2' 'tAretr_F_2' 'tMAretr_R_2' 'tMAretr_F_2'}
for cond = conds
    load([all_subs_folder 'GA_' cond{:} '_all'])
    for sub = 1:length(subjects)
       
        eval([' tmp  = ' [[['all_' cond{:}] '{' num2str(sub) '}']]]')
        dummy = find(tmp.time>=-0.18 & tmp.time <=1.23); % take samples after stim onset up to 1 s poststimulus
        maxpupil = max(tmp.avg(dummy));
        peak_latency = find(tmp.avg==maxpupil);
        eval( ['matrix_peaks' '{' num2str(sub) '}.subject' '= subjects(sub)']');
        eval( [['matrix_peaks' '{' num2str(sub) '}.' cond{:} 'peak'] '= maxpupil ']');
        eval( [['matrix_peaks' '{' num2str(sub) '}.' cond{:} 'latency'] '=  peak_latency']');
    end
end
cd(all_subs_folder)
matrix_peaks1 = cell2mat(matrix_peaks)
writetable(struct2table(matrix_peaks1), 'Peaks082021.xlsx')


end