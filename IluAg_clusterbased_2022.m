%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Cluster-based permutation tests for IluAg (Stefanie, 2022) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set paths
clear
ft_defaults
eeglab
close all

Exp_Code = '/BLB_EXP_201705_IluAg';
hd = '/DATA3';
RawData_folder = [hd Exp_Code '/BLB_BackUp_files/Raw_Data/'];
Analysis_folder = [hd Exp_Code '/Analysis'];
scripts_folder = [hd Exp_Code '/BLB_BackUp_files/Protocols/Analysis_Scripts/'];
elecs_file = [hd Exp_Code '/BLB_BackUp_files/Protocols/Configuration/IluAg.asc'];
anal_logfile = [Analysis_folder '/analysis_log.txt'];
Eeprobe_folder = [ Analysis_folder '/Eeprobe/experimento_IA/'];

addpath(genpath('/DATA3/BLB_EXP_201705_IluAg/BLB_BackUp_files/Protocols/Analysis_Scripts/plugins'));

subarray = [4 6 7 8 9 10 11 12 13 14 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31]


%% Adapt the data for Fieldtrip

indir = [Analysis_folder '/ClusterBased/'];

cfg=[];
for subj = 1:length(subarray)
    for cond = {'CONc','INCc','YESc','NOc','POSCONc','POSINCc','NEGCONc','NEGINCc','NEGc','POSc','TWonec','TWtwoc','TWthreec','TWfourc'}  
        EEG=pop_loadset([Analysis_folder num2str(subarray(subj),'%0.2d') '_' cond{:} '.set']);
        erpdata=eeglab2fieldtrip(EEG,'timelockanalysis','none');
        save([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat'],'erpdata');        
    end
end

%% 2)Prepare neighbours
%%% Do this once to get the layout and the neighbors and save them

load([indir '04_CONc_ERP.mat']);
cfg=[];
cfg.rotate = 90;
layout=ft_prepare_layout(cfg,erpdata); 
cfg=[];
cfg.method='triangulation';
cfg.layout=layout;
neighbours=ft_prepare_neighbours(cfg,erpdata);
cd(indir)
save('layout.mat','layout');
save('neighbours.mat', 'neighbours');
cd(scripts_folder)
        
load([indir 'GA_YESc_ERP.mat']);
cfg=[];
%cfg.rotate = 90;
layout2=ft_prepare_layout(cfg,erpdata); 
cfg=[];
cfg.method='triangulation';
cfg.layout=layout2;
neighbours2=ft_prepare_neighbours(cfg,erpdata);
cd(indir)
save('layout2.mat','layout2');
save('neighbours2.mat', 'neighbours2');
%% 3) Perform Grand Averages of ERPs

indir = 'D:\IluAg\ClusterBased\';

cfg=[]
    for cond = {'CONc','INCc','YESc','NOc','POSCONc','POSINCc','NEGCONc','NEGINCc','NEGc','POSc','TWonec','TWtwoc','TWthreec','TWfourc'}  
        i=1;
                for subj = 1:length(subarray)
                    
                    load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                    eval(['data_' num2str(i) ' = erpdata;']);
                    i = i+1;
                end
                cfg.keepindividual = 'no';
                % In this following line you should write as many data_X as
                % subjects that you are including in the GA
                eval(['GA_' cond{:} '_ERP = ft_timelockgrandaverage(cfg,data_1,data_2,data_3,data_4,data_5,data_6,data_7,data_8,data_9,data_10,data_11, data_12, data_13, data_14, data_15, data_16,data_17,data_18, data_19, data_20, data_21, data_22, data_23, data_24, data_25, data_26)']);
                eval(['save(''' indir  'GA_' cond{:} '_ERP'',''GA_' cond{:} '_ERP'')']);
                eval(['clear(''GA_' cond{:} '_ERP'')']);
                
    end

%% Load the preprocessed data

subarray = [4 6 7 8 9 10 11 12 13 14 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31]
indir = 'D:\IluAg\ClusterBased\';

CON=cell(length(subarray),1);  
INC=cell(length(subarray),1);
YES=cell(length(subarray),1);  
NO=cell(length(subarray),1);
POS=cell(length(subarray),1);
NEG=cell(length(subarray),1);
POSCON=cell(length(subarray),1);
POSINC=cell(length(subarray),1);
NEGCON=cell(length(subarray),1);
NEGINC=cell(length(subarray),1);
TW1=cell(length(subarray),1);
TW2=cell(length(subarray),1);
TW3=cell(length(subarray),1);
TW4=cell(length(subarray),1);

for cond = {'CONc','INCc','YESc','NOc','POSCONc','POSINCc','NEGCONc','NEGINCc','NEGc','POSc','TWonec','TWtwoc','TWthreec','TWfourc'}
    switch cond{:}
        case 'CONc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                CON{i}=erpdata;
                i = i+1;
            end
        case 'INCc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                INC{i}=erpdata;
                i = i+1;
            end
        case 'YESc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                YES{i}=erpdata;
                i = i+1;
            end
        case 'NOc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                NO{i}=erpdata;
                i = i+1;
            end
        case 'POSc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                POS{i}=erpdata;
                i = i+1;
            end
        case 'NEGc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                NEG{i}=erpdata;
                i = i+1;
            end
        case 'POSCONc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                POSCON{i}=erpdata;
                i = i+1;
            end
        case 'POSINCc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                POSINC{i}=erpdata;
                i = i+1;
            end
        case 'NEGCONc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                NEGCON{i}=erpdata;
                i = i+1;
            end
        case 'NEGINCc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                NEGINC{i}=erpdata;
                i = i+1;
            end
        case 'TWonec'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                TW1{i}=erpdata;
                i = i+1;
            end
        case 'TWtwoc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                TW2{i}=erpdata;
                i = i+1;
            end
        case 'TWthreec'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                TW3{i}=erpdata;
                i = i+1;
            end
        case 'TWfourc'
            i=1;
            for subj = 1:length(subarray)
                load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
                TW4{i}=erpdata;
                i = i+1;
            end
    end
end

% 
% for cond = {'YESc' 'NOc'}
%     switch cond{:}
%         case 'YESc'
%             i=1;
%             for subj = 1:length(subarray)
%                 load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
%                 YES{i}=erpdata;
%                 i = i+1;
%             end
%         case 'NOc'
%             i=1;
%             for subj = 1:length(subarray)
%                 load([indir num2str(subarray(subj),'%0.2d') '_' cond{:} '_ERP.mat']);
%                 NO{i}=erpdata;
%                 i = i+1;
%             end
%     end
% end

%% Configuration for statistics

cd(indir)
load layout
load neighbours
% cd(scripts_folder)
cfg = [];
cfg.parameter        = 'avg'; 
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; 
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'parametric'; %added by Nad
cfg.minnbchan        = 1;  %Marta had 1
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025; 
cfg.numrandomization = 1000; % Pump up
cfg.correcttail = 'prob';
cfg.layout = layout; 
cfg.neighbours = neighbours; 

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
cfg.ivar=2; cfg.uvar=1;

% Windows
N1 = [0.08 0.12];
N2 = [0.18 0.22];
P3 = [0.3 0.4];
P3_short = [0.34 0.37];

% N1 
cfg.latency = N1;
yes_no_N1 = ft_timelockstatistics(cfg, YES{:}, NO{:});
pos_neg_N1 = ft_timelockstatistics(cfg, POS{:}, NEG{:});
con_inc_N1 = ft_timelockstatistics(cfg, CON{:}, INC{:}); % nada
poscon_posinc_N1 = ft_timelockstatistics(cfg, POSCON{:}, POSINC{:});
negcon_neginc_N1 = ft_timelockstatistics(cfg, NEGCON{:}, NEGINC{:});
poscon_negcon_N1 = ft_timelockstatistics(cfg, POSCON{:}, NEGCON{:});
posinc_neginc_N1 = ft_timelockstatistics(cfg, POSINC{:}, NEGINC{:});
TW1_TW2_N1 = ft_timelockstatistics(cfg, TW1{:}, TW2{:});
TW2_TW3_N1 = ft_timelockstatistics(cfg, TW2{:}, TW3{:});
TW3_TW4_N1 = ft_timelockstatistics(cfg, TW3{:}, TW4{:});
 
% N2 
cfg.latency = N2;
yes_no_N2 = ft_timelockstatistics(cfg, YES{:}, NO{:});
pos_neg_N2 = ft_timelockstatistics(cfg, POS{:}, NEG{:});
con_inc_N2 = ft_timelockstatistics(cfg, CON{:}, INC{:}); 
poscon_posinc_N2 = ft_timelockstatistics(cfg, POSCON{:}, POSINC{:});
negcon_neginc_N2 = ft_timelockstatistics(cfg, NEGCON{:}, NEGINC{:});
poscon_negcon_N2 = ft_timelockstatistics(cfg, POSCON{:}, NEGCON{:});
posinc_neginc_N2 = ft_timelockstatistics(cfg, POSINC{:}, NEGINC{:});
TW1_TW2_N2 = ft_timelockstatistics(cfg, TW1{:}, TW2{:});
TW2_TW3_N2 = ft_timelockstatistics(cfg, TW2{:}, TW3{:});
TW3_TW4_N2 = ft_timelockstatistics(cfg, TW3{:}, TW4{:});

% P3 
cfg.latency = P3;
yes_no_P3 = ft_timelockstatistics(cfg, YES{:}, NO{:});
pos_neg_P3 = ft_timelockstatistics(cfg, POS{:}, NEG{:});
poscon_posinc_P3 = ft_timelockstatistics(cfg, POSCON{:}, POSINC{:});
poscon_negcon_P3 = ft_timelockstatistics(cfg, POSCON{:}, NEGCON{:});
posinc_neginc_P3 = ft_timelockstatistics(cfg, POSINC{:}, NEGINC{:});
TW2_TW3_P3 = ft_timelockstatistics(cfg, TW2{:}, TW3{:});
TW3_TW4_P3 = ft_timelockstatistics(cfg, TW3{:}, TW4{:});

% P3_short
cfg.latency = P3_short;
con_inc_P3 = ft_timelockstatistics(cfg, CON{:}, INC{:}); 
negcon_neginc_P3 = ft_timelockstatistics(cfg, NEGCON{:}, NEGINC{:});
TW1_TW2_P3 = ft_timelockstatistics(cfg, TW1{:}, TW2{:});


cfg=[];
cfg.layout=layout;
cfg.highlightcolorpos = [1 1 1];
cfg.highlightcolorneg = [0 0 0];
% cfg.highlightsymbolseries=['*','*','*','*','*'];
cfg.zlim='maxabs';
cfg.parameter='stat'
cfg.maskparameter = 'mask';
cfg.alpha=0.050;
ft_clusterplot(cfg,comparisons{28});
% ft_clusterplot(cfg,yes_no_N2);
% ft_clusterplot(cfg,yes_no_P3);
% ft_clusterplot(cfg, pos_neg_N2);
% ft_clusterplot(cfg, poscon_posinc_N1); % if we raise the alpha for plotting to 0.052 we can see an almost significant, left lateralised N1 attenuation
% ft_clusterplot(cfg, negcon_neginc_P3);

%% Plot a single scalp plot with significant electrodes highlighted 
% cfg = [];
% cfg.operation = 'subtract';
% cfg.parameter = 'avg';
% GA_YES_NO = ft_math(cfg,GA_YESc_ERP, GA_NOc_ERP)
% 
% figure;
% % define parameters for plotting
% timestep      = 0.05; %(in seconds)
% sampling_rate = 500;
% sample_count  = length(yes_no_N1.time);
% j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
% m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples

%% One scalp plot with an average of activation and significant electrodes highlighted

% YES-NO-N1

% Specify which comparison to plot and which cluster you want
% comparison = yes_no_N1;

comparisons = {yes_no_N1,  pos_neg_N1, con_inc_N1, poscon_posinc_N1, negcon_neginc_N1, poscon_negcon_N1, posinc_neginc_N1, TW1_TW2_N1, TW2_TW3_N1, TW3_TW4_N1, yes_no_N2, pos_neg_N2, con_inc_N2, poscon_posinc_N2, negcon_neginc_N2, poscon_negcon_N2, posinc_neginc_N2, TW1_TW2_N2, TW2_TW3_N2, TW3_TW4_N2, yes_no_P3, pos_neg_P3, con_inc_P3, poscon_posinc_P3, negcon_neginc_P3, poscon_negcon_P3, posinc_neginc_P3, TW1_TW2_P3, TW2_TW3_P3, TW3_TW4_P3}
names = {'yes_no_N1', ' pos_neg_N1', 'con_inc_N1', 'poscon_posinc_N1', 'negcon_neginc_N1', 'poscon_negcon_N1', 'posinc_neginc_N1', 'TW1_TW2_N1', 'TW2_TW3_N1', 'TW3_TW4_N1', 'yes_no_N2', 'pos_neg_N2', 'con_inc_N2', 'poscon_posinc_N2', 'negcon_neginc_N2', 'poscon_negcon_N2', 'posinc_neginc_N2', 'TW1_TW2_N2', 'TW2_TW3_N2', 'TW3_TW4_N2', 'yes_no_P3', 'pos_neg_P3', 'con_inc_P3', 'poscon_posinc_P3', 'negcon_neginc_P3', 'poscon_negcon_P3', 'posinc_neginc_P3', 'TW1_TW2_P3', 'TW2_TW3_P3', 'TW3_TW4_P3'}
windows = {N1, N1, N1, N1, N1, N1, N1, N1, N1, N1, N2, N2, N2, N2, N2, N2, N2, N2, N2, N2, P3, P3, P3_short, P3, P3_short, P3, P3, P3_short, P3, P3}

%for i = 1:length(comparisons)
for i = 23
    % check if there are any clusters at all
    if isfield(comparisons{i},'posclusters')
        % get relevant values
        pos_cluster_pvals = [comparisons{i}.posclusters(:).prob];
        pos_clust = find(pos_cluster_pvals < 0.05);
        pos       = ismember(comparisons{i}.posclusterslabelmat, pos_clust);
        neg_cluster_pvals = [comparisons{i}.negclusters(:).prob];
        neg_clust = find(neg_cluster_pvals < 0.05);
        neg       = ismember(comparisons{i}.negclusterslabelmat, neg_clust);

        if pos_clust == 1
            highlight = 1;
        elseif neg_clust == 1
            highlight = 0;
        else
            highlight = [];
            continue
        end

        if highlight == 0
            color = [0 0 0];
            highlightstat = neg;
        else
            color = [1 1 1];
            highlightstat = pos;
        end

        cfg                     = [];
        cfg.parameter           = 'stat';
        cfg.layout              = layout;
        cfg.xlim                = windows{i};
        cfg.maskparameter       = 'mask';
        cfg.zlim                = 'maxabs';
        cfg.comment             = 'no';
        %cfg.colorbar            = 'yes';
        cfg.highlight           = 'on';
        cfg.highlightchannel    = {find(highlightstat)};
%         cfg.highlightcolor      = color;
        cfg.highlightsize       = 12;
        f = figure;
        f.Position = [500 500 270 200];
        ft_topoplotER(cfg, comparisons{i});
        %colorbar('FontSize',20);
        colorbar;
        saveas(f, [sprintf(names{i}) '.jpg'])
    else
        continue
    end
end 

% YES-NO-N2

% Specify which comparison to plot and which cluster you want
comparison = yes_no_N1;
window = N2;

% get relevant values
pos_cluster_pvals = [comparison.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.05);
pos       = ismember(comparison.posclusterslabelmat, pos_clust);
neg_cluster_pvals = [comparison.negclusters(:).prob];
neg_clust = find(neg_cluster_pvals < 0.05);
neg       = ismember(comparison.negclusterslabelmat, neg_clust);

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
% [i1,i2] = match_str(GA_YES_NO.label, yes_no_N1.label);

% Specify here what to plot
highlight = pos;

if highlight == neg
    color = [0 0 0];
else 
    color = [1 1 1];
end 

cfg                     = [];
cfg.parameter           = 'stat';
cfg.layout              = layout;
cfg.xlim                = window;
cfg.maskparameter       = 'mask';
cfg.zlim                = 'maxabs';
cfg.comment             = 'no';
%cfg.colorbar            = 'yes';
cfg.highlight           = 'on';
cfg.highlightchannel    = find(highlight);
cfg.highlightcolor      = color;
cfg.highlightsize       = 12;
f = figure;
f.Position = [500 500 270 200];
ft_topoplotER(cfg, comparison);
%colorbar('FontSize',20);
colorbar;
