
% ILUAG: CLUSTER-BASED PERMUTATION pos_negS

% Script for cluster-based permutation pos_negs in IluAg
% Created by Stefanie, December 2021
% Adapted from Marta

%% Paths and codes
clear
ft_defaults
eeglab
close all

Exp_Code = '';
hd = '';
RawData_folder = '/DATA3/BLB_EXP_201705_IluAg/BLB_BackUp_files/Raw_Data/';
Analysis_folder = '/DATA3/BLB_EXP_201705_IluAg/Analysis/';
scripts_folder = '/DATA3/BLB_EXP_201705_IluAg/BLB_BackUp_files/Protocols/Analysis_Scripts/';
elecs_file = '/DATA3/BLB_EXP_201705_IluAg/BLB_BackUp_files/Protocols/Configuration/IluAg.asc';
anal_logfile = [Analysis_folder '/analysis_log.txt'];
Eeprobe_folder = '/DATA3/BLB_EXP_201705_IluAg/Analysis/Eeprobe/experimento_IA/';

addpath(genpath('/DATA3/BLB_EXP_201705_IluAg/BLB_BackUp_files/Protocols/Analysis_Scripts/plugins'));

%% 1) Export to fieldtrip for ERP analyses

subarray = [4 6 7 8 9 10 11 12 13 14 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31]
% participants 5, 15 excluded

indir = [Analysis_folder 'ClusterBased/'];

cfg=[];
for subj = 1:length(subarray)
    for cond = {'CONc','INCc','YESc','NOc','POSCONc','POSINCc','NEGCONc','NEGINCc','NEGc','POSc','TWonec','TWtwoc','TWthreec','TWfourc'}  
        EEG=pop_loadset([Analysis_folder num2str(subarray(subj),'%0.2d') '_' cond{:} '.set']);
        % EEG= pop_selectevent(EEG, 'latency','-0.01<=0.01','deleteevents','on','deleteepochs','on','invertepochs','off');%delete duplicates if any
        % EEG=pop_rmbase(EEG,[-100 0]); %not sure why
        erpdata=eeglab2fieldtrip(EEG,'timelockanalysis','none');
        %erpdata=ft_timelockanalysis(cfg,data);
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
        
% Correct Motor 
% indir = 'C:\Users\Marta\Documents\BLB_EXP_201706_SGMem\Analysis\';
% subarray = [2 3 4 5 6 7 8 9 10 11];
% 
% conditions to subtract 
% condArray1 = {'eMA'}; %   % 1 minus 2
% condArray2 = {'eM'};
% 
% for subj = 1:length(subarray)
%     for iDiff = 1:length(condArray1)
%         
%         load([indir num2str(subarray(subj),'%0.2d') '_' condArray1{iDiff} '_ERP']);
%         cond1 = erpdata;
%         load([indir num2str(subarray(subj),'%0.2d') '_' condArray2{iDiff} '_ERP']);
%         cond2 = erpdata;
%         
%         cfg = [];
%         cfg.operation = 'subtract';
%         cfg.parameter = 'avg';
%         
%         erpdata = ft_math(cfg, cond1, cond2)
%         save([indir num2str(subarray(subj),'%0.2d') '_' condArray1{iDiff} '_' condArray2{iDiff} '_ERP.mat'],'erpdata');
%         
%     end
% end

% %% Raw effect (for viz later): eA-[eMA-eM]
% % conditions to subtract 
% condArray1 = {'eA'}; %   % 1 minus 2
% condArray2 = {'eMA_eM'};
% 
% for subj = 1:length(subarray)
%     for iDiff = 1:length(condArray1)
%         
%         load([indir num2str(subarray(subj),'%0.2d') '_' condArray1{iDiff} '_ERP']);
%         cond1 = erpdata;
%         load([indir num2str(subarray(subj),'%0.2d') '_' condArray2{iDiff} '_ERP']);
%         cond2 = erpdata;
%         
%         cfg = [];
%         cfg.operation = 'subtract';
%         cfg.parameter = 'avg';
%         
%         raw_effect = ft_math(cfg, cond1, cond2)
%         save([indir num2str(subarray(subj),'%0.2d') '_' condArray1{iDiff} '_' condArray2{iDiff} '_ERP.mat'],'raw_effect');
%         
%     end
% end

%% 3) Perform Grand Averages of ERPs
% Do not keep individuals

%clear
%subarray = [2 3 4 5 6 7 8 9 10 11];
%indir = 'C:\Users\Marta\Documents\BLB_EXP_201706_SGMem\Analysis\';
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
   
% %% 4) ERP PLOTS
% cd(indir)
% load layout
% load GA_CONc_ERP
% load GA_INCc_ERP
% cfg=[];
% cfg.layout=layout;
% cfg.xlim=[-0.1 0.5];
% cfg.ylim=[-4 4];
% cfg.ydir='reverse';
% figure;
% ft_multiplotER(cfg,GA_CONc_ERP, GA_INCc_ERP)


%% 




%% STATISTICS
%%
% ERP POST-HOC ANALYSES (Differences between the conditions for young/old)
% YOUNG
CON=cell(length(subarray),1); % 
INC=cell(length(subarray),1);
YES=cell(length(subarray),1); % 
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

%% Montecarlo 

% Choose the window you are interested in 
window = [0 0.4];

cd(indir)
load layout
load neighbours
cd(scripts_folder)
cfg = [];
cfg.parameter        = 'avg'; 
cfg.latency          = window; %[0 0.4]; [0.050 0.150];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; 
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'parametric'; %added by Nad
cfg.minnbchan        = 1;  %Marta had 1
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.05; 
cfg.numrandomization = 500; % Pump up
cfg.correcttail = 'prob';

cfg.layout = layout; 
cfg.neighbours = neighbours; 

%cfg.design = [ ones(1,length(subarray)) ones(1,length(subarray))*2; 1:length(subarray) 1:length(subarray) ];
%cfg.ivar     = 1; % Row in which the independent variable (group number) is set
%cfg.uvar     = 23; % Row in which the unit of observation variable (subject number) is set

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

CON_INC = ft_timelockstatistics(cfg, CON{:}, INC{:});

% Make different comparisons

% POS vs. NEG
[pos_neg] = ft_timelockstatistics(cfg, POS{:}, NEG{:}); 
save ([indir 'pos_neg.mat'], 'pos_neg')
% 
% % YES vs. NO
% [yes_no] = ft_timelockstatistics(cfg, YES{:}, NO{:}); 
% save ([indir 'yes_no.mat'], 'yes_no')

% %% Stats for retrieval
% cd(indir)
% load layout
% load neighbours
% cd(scripts_folder)
% cfg = [];
% cfg.parameter        = 'avg'; 
% cfg.latency          = [0.050 0.150]; %N1 first
% cfg.method           = 'montecarlo';
% cfg.statistic        = 'depsamplesT';
% cfg.correctm         = 'cluster';
% cfg.clusteralpha     = 0.05; 
% cfg.clusterstatistic = 'maxsum';
% cfg.clusterthreshold = 'parametric'; %added by Nad
% cfg.minnbchan        = 1;  %Marta had 1
% cfg.tail             = 0;
% cfg.clustertail      = 0;
% cfg.alpha            = 0.05; 
% cfg.numrandomization = 10000; % Pump up
% cfg.correcttail = 'prob';
% 
% cfg.layout = layout; 
% cfg.neighbours = neighbours; 
% 
% %cfg.design = [ ones(1,length(subarray)) ones(1,length(subarray))*2; 1:length(subarray) 1:length(subarray) ];
% %cfg.ivar     = 1; % Row in which the independent variable (group number) is set
% %cfg.uvar     = 23; % Row in which the unit of observation variable (subject number) is set
% 
% subj = length(subarray);
% design = zeros(2,2*subj);
% for i = 1:subj
%     design(1,i) = i;
% end
% for i = 1:subj
%     design(1,subj+i)=i;
% end
% design(2,1:subj)=1; design(2,subj+1:2*subj)=2;
% cfg.design = design;
% cfg.ivar=2; cfg.uvar=1;
% % HERE YOU WRITE THE NAME YOU WANT OF THE STATS DATA STRUCTURE
% [ERP_tAvstMA_N1] = ft_timelockstatistics(cfg, C3{:}, C4{:}); 
% save ([indir 'ERP_tAvstMA_N1.mat'], 'ERP_tAvstMA_N1')
% 
% %and P2 
% cfg.latency = [0.17 0.3];
% [ERP_tAvstMA_P2] = ft_timelockstatistics(cfg, C3{:}, C4{:}); 
% save ([indir 'ERP_tAvstMA_P2.mat'], 'ERP_tAvstMA_P2')

% %% STATS PLOTS - Encoding first
% cd(indir)
% load layout
% cfg=[];
% cfg.layout=layout;
% 
% % cfg.highlightsymbolseries=['*','*','*','*','*'];
% cfg.zlim='maxabs';
% cfg.parameter='stat'
% cfg.maskparameter = 'mask';
% cfg.highlight = 'on';
% cfg.alpha=0.05;
% %cfg.colorbar ='yes';
% ft_clusterplot(cfg,pos_neg);
% 
% ft_clusterplot(cfg, ERP_eAvseMA_eM_P2);
%% From these stat plots we know now what time windows we want to plot. in this case, from 650ms to 1000ms in steps of 50ms.
%AGExCOND_T2BvsT1B_YvsO_2sec_18062019.maskedstat=AGExCOND_T2BvsT1B_YvsO_2sec_18062019.stat .* AGExCOND_T2BvsT1B_YvsO_2sec_18062019.mask;
% figure;
% cfg=[];
% cfg.layout=layout;
% %cfg.xlim=[.200 .400];
% cfg.colorbar='WestOutside';
% cfg.parameter='stat';
% cfg.highlightcolorpos = [1 1 1];
% cfg.highlightcolorneg = [0.5 0.5 0.5];
% 
% cfg.zlim= 'maxabs';%'absmax';
% cfg.latency =[.05 .150];
% %cfg.maskparameter = 'mask';
% cfg.highlight ='on';
% cfg.opacitymap='rampup';
% cfg.gridscale=300;
% cfg.colormap = 'jet';
% cfg.shading = 'interp';
% cfg.style='fill';
% cfg.interplimits='electrodes';
% cfg.interpolatenan='no';
% ft_topoplotER(cfg,ERP_eAvseMA_eM_P2)

%% Create plots for encoding
cd(indir); load GA_CONc_ERP; load GA_INCc_ERP;

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_CONvsINC   = ft_math(cfg, GA_CONc_ERP, GA_INCc_ERP);

figure;
% define params for plotting 
timestep = 0.05 %in secs
sampling_rate = 500; %
sample_count = length(pos_neg.time);
% j and M must have the same temporal length
j = [pos_neg.time(1):timestep:pos_neg.time(end)]; %temporal endpoints (in secs) of the ERP avg computed in each subplot
m = [1:timestep*sampling_rate:sample_count]; %temporal endpoints in EEG samples
% get relevant (significant) values for positive
pos_cluster_pvals = [pos_neg.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < pos_neg.cfg.alpha);
pos = ismember(pos_neg.posclusterslabelmat, pos_signif_clust);

% get relevant (significant) values for negative
neg_cluster_pvals = [pos_neg.negclusters(:).prob];
neg_signif_clust = find(neg_cluster_pvals < pos_neg.cfg.alpha);
neg = ismember(pos_neg.negclusterslabelmat, neg_signif_clust);

% First ensure the channels to have the same order in the avg and in the
% stat output. This might not be the case, because ft_math might shuffle
% the order
[i1,i2] = match_str(GA_POSvsNEG.label, pos_neg.label);
    
% plot pliz for positive/negative clusters
ncols = 2
for k = 1:length(m)-1
    subplot((0.4/timestep)/ncols,ncols,k);
    cfg=[];
    cfg.xlim = [j(k) j(k+1)];
   % cfg.zlim = [-5e-14 5e-14];
    cfg.colorbar = 'EastOutside';
    pos_int = zeros(numel(GA_POSvsNEG.label),1);
    pos_int(i1) = all(pos(i2, m(k):m(k+1)),2);
    neg_int = zeros(numel(GA_POSvsNEG.label),1);
    neg_int(i1) = all(neg(i2, m(k):m(k+1)),2);
    cfg.highlight = 'on';
    cfg.highlightchannel = find(pos_int | neg_int);
  
    cfg.comment = 'xlim';
    cfg.maskparameter = 'mask';
    cfg.commentpos = 'title';
    cfg.interactive ='no';
    cfg.layout = layout;
    ft_topoplotER(cfg,GA_POSvsNEG);
end

%% Create plots for retrieval
% Cluster plot first
cd(indir)
load layout
cfg=[];
cfg.layout=layout;
cfg.latency = [0.05 0.3]
cfg.highlightcolorpos = [1 1 1];
cfg.highlightcolorneg = [0.5 0.5 0.5];
% cfg.highlightsymbolseries=['*','*','*','*','*'];
cfg.zlim='maxabs';
cfg.parameter='stat'
cfg.maskparameter = 'mask';
cfg.alpha=0.05;
%cfg.colorbar ='yes';
ft_clusterplot(cfg,ERP_tAvstMA);
% then topoplot
cd(indir); load GA_tAretr_ERP; load GA_tMAretr_ERP;

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_tAretr_vs_tMAretr = ft_math(cfg,GA_tAretr_ERP, GA_tMAretr_ERP)

figure;
% define params for plotting 
timestep = 0.05 %in secs
sampling_rate = 500; %
sample_count = length(ERP_tAvstMA.time);
% j and M must have the same temporal length
j = [0:timestep:ERP_tAvstMA.time(end)]; %temporal endpoints (in secs) of the ERP avg computed in each subplot
m = [1:timestep*sampling_rate:sample_count]; %temporal endpoints in EEG samples
% get relevant (significant) values for positive
pos_cluster_pvals = [ERP_tAvstMA.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < ERP_tAvstMA.cfg.alpha);
pos = ismember(ERP_tAvstMA.posclusterslabelmat, pos_signif_clust);

% get relevant (significant) values for negative (if any)
if ~isempty(ERP_tAvstMA.negclusters)
    neg_cluster_pvals = [ERP_tAvstMA.negclusters(:).prob];
    neg_signif_clust = find(neg_cluster_pvals < ERP_tAvstMA.cfg.alpha);
    neg = ismember(ERP_tAvstMA.negclusterslabelmat, neg_signif_clust);

end
% First ensure the channels to have the same order in the avg and in the
% stat output. This might not be the case, because ft_math might shuffle
% the order
[i1,i2] = match_str(GA_tAretr_vs_tMAretr.label, ERP_tAvstMA.label);
    
% plot pliz for positive/negative clusters
for k = 1:length(m)-1
    subplot(2,5,k);
    cfg=[];
    cfg.xlim = [j(k) j(k+1)];
   % cfg.zlim = [-5e-14 5e-14];
    pos_int = zeros(numel(GA_tAretr_vs_tMAretr.label),1);
    pos_int(i1) = all(pos(i2, m(k):m(k+1)),2);
%     neg_int = zeros(numel(GA_eA_vs_eMA_eM.label),1);
%     neg_int(i1) = all(neg(i2, m(k):m(k+1)),2);
    cfg.highlight = 'on';
    %cfg.highlightchannel = find(pos_int);
    cfg.comment = 'xlim';
    cfg.maskparameter = 'mask';
    cfg.commentpos = 'title';
    %cfg.interactive ='no';
    cfg.layout = layout;
    ft_topoplotER(cfg,GA_tAretr_vs_tMAretr);
end
% %% try other plots
% cfg = [];
% cfg.layout= 'ordered';
% layout2=ft_prepare_layout(cfg,GA_POSc_ERP)
% cfg=[];
% cfg.channel = 'all' % Cz is 28
% ft_singleplotER (cfg, GA_POSc_ERP, GA_NEGc_ERP)
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

%% 
%%% PLOTS BY STEFANIE %%%

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_POSvsNEG    = ft_math(cfg, GA_POSc_ERP, GA_NEGc_ERP);

figure;
% define parameters for plotting
timestep      = 0.05; %(in seconds)
sampling_rate = POS{1}.fsample;
sample_count  = length(pos_neg.time);
j = [0:timestep:1];   % Temporal endpoints (in seconds) of the ERP average computed in each subplot
m = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in M/EEG samples

% get relevant values
pos_cluster_pvals = [pos_neg.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(pos_neg.posclusterslabelmat, pos_clust);

% First ensure the channels to have the same order in the average and in the statistical output.
% This might not be the case, because ft_math might shuffle the order
[i1,i2] = match_str(GA_POSvsNEG.label, pos_neg.label);

% plot
for k = 1:20;
   cfg.figure     = subplot(4,5,k);
   cfg.xlim       = [j(k) j(k+1)];
   cfg.zlim       = [-5e-14 5e-14];
   pos_int        = zeros(numel(GA_POSvsNEG.label),1);
   pos_int(i1)    = all(pos(i2, m(k):m(k+1)), 2);
   cfg.highlight  = 'on';
   cfg.highlightchannel = find(pos_int);
   cfg.comment    = 'xlim';
   cfg.commentpos = 'title';
   cfg.layout     = layout;
   ft_topoplotER(cfg, GA_POSvsNEG);
end

        