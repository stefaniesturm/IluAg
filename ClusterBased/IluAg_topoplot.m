function IluAg_topoplot(highlight)



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