function [pos, neg] = GetHighlight(comparison)
pos_cluster_pvals = [comparison.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.05);
pos       = ismember(comparison.posclusterslabelmat, pos_clust);
neg_cluster_pvals = [comparison.negclusters(:).prob];
neg_clust = find(neg_cluster_pvals < 0.05);
neg       = ismember(comparison.negclusterslabelmat, neg_clust);
