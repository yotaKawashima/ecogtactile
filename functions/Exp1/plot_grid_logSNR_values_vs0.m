function plot_grid_logSNR_values_vs0(fig_name, logSNR, trials, bipolar_reref, LeftF0, LeftIM, RightF0, RightIM, faxis)
    % Plot grid longSNR values vs0 (Thomas ver)
    figure;
    set(gcf,'Name',fig_name)
    cmap=colormap('hot');
    cmap=flipud(cmap);
    for nF=1:length(LeftF0)
        subplot(2,4,nF);
        [thisF,idxF]=findclosest(faxis,LeftF0(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('contraF (%g)',LeftF0(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(LeftIM)
        subplot(2,4,nF+2);
        [thisF,idxF]=findclosest(faxis,LeftIM(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('contraIM (%g)',LeftIM(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(RightF0)
        subplot(2,4,nF+4);
        [thisF,idxF]=findclosest(faxis,RightF0(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('ipsiF (%g)',RightF0(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    for nF=1:length(RightIM)
        subplot(2,4,nF+6);
        [thisF,idxF]=findclosest(faxis,RightIM(nF));
        [h, pV, ~, stats]=ttest(squeeze(logSNR(:,idxF,trials)),0,'dim',2);
        [p, FV] = draw_biprref(mean(stats.tstat,3), bipolar_reref, [5 4], [-1 8]);
        title(sprintf('ipsiIM (%g)',RightIM(nF)))
        format_fig;
        colormap(cmap);
        colorbar;
    end
    set(gcf, 'Position', get(0, 'Screensize'));
    
end
