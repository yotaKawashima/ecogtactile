function plot_spectrum(fig_name, logSNR, logPow, trials, bipolar_reref, LeftF0, LeftIM, RightF0, RightIM, faxis)
    % Plot spectrum (Thomas ver)
    figure;
    set(gcf,'Name',fig_name)
    subplot(1,2,1); 
    plot(faxis,mean(logPow(bipolar_reref(:,3)==0,:,trials),3),'Color','k','LineWidth',2)
    hold on;
    plot(faxis,mean(logPow(bipolar_reref(:,3)==1,:,trials),3),'Color',[1 1 1]*0.5,'LineWidth',2)
    xlim([2 45]);
    title({'log(Power) - all somatosens channels','bi-polar (black: vertical; gray: horiz)'})
    format_fig;
    
    subplot(1,2,2); 
    plot(faxis,mean(logSNR(bipolar_reref(:,3)==0,:,trials),3),'Color','k','LineWidth',2)
    hold on;
    plot(faxis,mean(logSNR(bipolar_reref(:,3)==1,:,trials),3),'Color',[1 1 1]*0.5,'LineWidth',2)
    xlim([2 45]);
    title({'log(SNR) - all somatosens channels','bi-polar (black: vertical; gray: horiz)'})
    format_fig;

    % mark fundamentals
    subplot(1,2,2); hold on;
    lgd=[];
    for nF=1:length(LeftF0)
        [thisF,idxF]=findclosest(faxis,LeftF0(nF));
        lgd(1)=scatter(thisF,max(mean(logSNR(:,idxF,trials),3),[],1),'filled','Marker','o','MarkerFaceColor','r','SizeData',100);
    end
    for nF=1:length(RightF0)
        [thisF,idxF]=findclosest(faxis,RightF0(nF));
        lgd(2)=scatter(thisF,max(mean(logSNR(:,idxF,trials),3),[],1),'filled','Marker','o','MarkerFaceColor','b','SizeData',100);
    end
    format_fig;

    % mark IM
    subplot(1,2,2); hold on;
    for nF=1:length(LeftIM)
        [thisF,idxF]=findclosest(faxis,LeftIM(nF));
        lgd(3)=scatter(thisF,max(mean(logSNR(:,idxF,trials),3),[],1),'filled','Marker','*','MarkerEdgeColor','r','SizeData',100,'LineWidth',1.5);
    end
    for nF=1:length(RightIM)
        [thisF,idxF]=findclosest(faxis,RightIM(nF));
        lgd(4)=scatter(thisF,max(mean(logSNR(:,idxF,trials),3),[],1),'filled','Marker','*','MarkerEdgeColor','b','SizeData',100,'LineWidth',1.5);
    end
    legend(lgd,{'contraF','ipsiF','contraIM','ipsiIM'})
    format_fig;
    set(gcf, 'Position', get(0, 'Screensize'));
end
