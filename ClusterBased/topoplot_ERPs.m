function topoplot_ERPs(stat_clu, data_for_plotting, title)

indir = 'D:\IluAg\ClusterBased\';
load('layout');
sampling_rate = 500; %
timestep=0.05;
sample_count = length(stat_clu.time);
% j and M must have the same temporal length
j = [stat_clu.time(1):timestep:stat_clu.time(end)]; %temporal endpoints (in secs) of the ERP avg computed in each subplot
m = [1:timestep*sampling_rate:sample_count]
% get relevant (significant) values for positive
if isfield(stat_clu,'posclusters')
    pos_cluster_pvals = [stat_clu.posclusters(:).prob];
    pos_signif_clust = find(pos_cluster_pvals < stat_clu.cfg.alpha);
    pos = ismember(stat_clu.posclusterslabelmat, pos_signif_clust);
else
    pos_signif_clust=[];
end
if isfield(stat_clu, 'negclusters')
    % get relevant (significant) values for negative
    neg_cluster_pvals = [stat_clu.negclusters(:).prob];
    neg_signif_clust = find(neg_cluster_pvals < stat_clu.cfg.alpha);
    neg = ismember(stat_clu.negclusterslabelmat, neg_signif_clust);
else
    neg_signif_clust=[];
end
[i1,i2] = match_str(data_for_plotting.label, stat_clu.label);
for k = 1:length(m)-1
    figure;
    %subplot(2,4,k);
    cfg = [];
    cfg.xlim =[j(k) j(k+1)];
    %  cfg.zlim = [-1.0e-13 1.0e-13];
    pos_int = zeros(numel(data_for_plotting.label),1);
    neg_int = zeros(numel(data_for_plotting.label),1);
    neg_int(i1) = all(neg(i2, m(k):m(k+1)), 2);
    pos_int(i1) = all(pos(i2, m(k):m(k+1)), 2);
    cfg.highlight = 'on';
    cfg.highlightchannel = find(pos_int|neg_int);
    cfg.highlightsymbol    = '*'
    cfg.highlightsize      = 4
    cfg.highlightfontsize  = 4
    cfg.comment = 'xlim';
    cfg.commentpos = 'title';
    %         cfg.colorbar='South';
    cfg.markersymbol       = '*'
    cfg.markersize         = 1;
    cfg.layout = layout;
    cfg.style = 'straight_imsat';
    cfg.colorbar='SouthOutSide'
    cfg.fontsize=20;
    ft_topoplotER(cfg, data_for_plotting);
    saveas(gcf,[title '_' num2str(k) '.png'])
end