clear all;

fprintf(' Creating figures .......\n');

load ecogProDataProcess;
load ecogGlobal.mat;

for nversion=1:numOfCondition
    RightF0 = (RightF0_all(:,nversion))';
    LeftF0 =  (LeftF0_all(:,nversion))';
    faxis = faxis_all(:,nversion);
    version=allversions{nversion};
    pow_bych = pow_bych_all(:,:,:,nversion);
    pow_gamma = pow_gamma_all(:,:,:,nversion);
    
    att_1 = pow_gamma(:,1:6,:); %condition 1 in epoch are numbered 1 till 6
    att_3 = pow_gamma(:,6:12,:); %condition 3 in epoch are numbered 6 till 12
    
    figure;
    set(gcf,'position',[-499 801 854 1056],'Name',version)
    subplot(4,2,1)
%     plot(faxis,log(squeeze(mean(mean(pow_bych([find(reref_mat(:,2)==S1channels(1) & reref_mat(:,3)==0) find(reref_mat(:,2)==S1channels(2) & reref_mat(:,3)==0)],:,:),1),2))),'Color','k','LineWidth',2)
    plot(faxis,log(squeeze(mean(mean(pow_bych([find(ismember(reref_mat(:,2),S1channels) & reref_mat(:,3)==0)],:,:),1),2))),'Color','k','LineWidth',2)
    xlim([1 45])
    xlabel('Freq (Hz)')
    ylabel('logPower')
% %     ylabel(sprintf('SNR (ch %g %g) bipolar I',S1channels(1),S1channels(2)))
    ylabel(sprintf('SNR (ch %g - %g) bipolar I',S1channels(1),S1channels(end)))
    
    subplot(4,2,2)
    %plot(faxis,log(squeeze(mean(mean(pow_bych([find(reref_mat(:,2)==S1channels(1) & reref_mat(:,3)==1) find(reref_mat(:,2)==S1channels(2) & reref_mat(:,3)==1)],:,:),1),2))),'Color','k','LineWidth',2)
    plot(faxis,log(squeeze(mean(mean(pow_bych([find(ismember(reref_mat(:,2),S1channels) & reref_mat(:,3)==1)],:,:),1),2))),'Color','k','LineWidth',2)
    xlim([1 45])
    xlabel('Freq (Hz)')
    ylabel('logPower')
% %     ylabel(sprintf('SNR (ch %g %g) bipolar II',S1channels(1),S1channels(2)))
    ylabel(sprintf('SNR (ch %g - %g) bipolar I',S1channels(1),S1channels(end)))
    
    %%
    % extract fundamentals
    fondFreq=[LeftF0 RightF0];
    pow_tag=[];
    pow_tag2=[];
    for nfreq=1:length(fondFreq)
        [~,findfreq]=findclosest(faxis,fondFreq(nfreq));
        pow_tag(:,:,nfreq)=log(pow_bych(:,:,findfreq))-1/2*(log(pow_bych(:,:,findfreq-1)) + log(pow_bych(:,:,findfreq+1)));
        pow_tag2(:,nfreq)=log(mean(pow_bych(:,:,findfreq),2))-1/2*(log(mean(pow_bych(:,:,findfreq-1),2)) + log(mean(pow_bych(:,:,findfreq+1),2)));
    end
    % extract IM
    fondFreq=[abs(LeftF0(2)-LeftF0(1)) abs(RightF0(2)-RightF0(1))];
    pow_IM=[];
    pow_IM2=[];
    for nfreq=1:length(fondFreq)
        [~,findfreq]=findclosest(faxis,fondFreq(nfreq));
        pow_IM(:,:,nfreq)=log(pow_bych(:,:,findfreq))-1/2*(log(pow_bych(:,:,findfreq-1)) + log(pow_bych(:,:,findfreq+1)));
        pow_IM2(:,nfreq)=log(mean(pow_bych(:,:,findfreq),2))-1/2*(log(mean(pow_bych(:,:,findfreq-1),2)) + log(mean(pow_bych(:,:,findfreq+1),2)));
    end
    
    % plot grids
%     posX=repmat(1:elecC,1,6);
    posX=repmat(1:elecC,1,elecR+1);
    posY=[];
    for ind =1:elecR +1
    posY = [posY ind * ones(1,elecC)];
    end %for
    %posY=[ones(1,elecC) 2*ones(1,elecC) 3*ones(1,elecC) 4*ones(1,elecC) 5*ones(1,elecC) 6*ones(1,elecC)];
    
    subplot(4,2,3)
    % titles={sprintf('f0=%g (contra)',LeftF0(1)),sprintf('f0=%g (contra)',LeftF0(2)),sprintf('f0=%g (ipsi)',RightF0(1)),sprintf('f0=%g (ipsi)',RightF0(2))};
    % for nfreq=1:4
    [h, pV, ~, stats]=ttest(squeeze(pow_tag(1:sizeH,:,1:2)),0,'dim',2);
    tp=mean(stats.tstat,3)/3.5*255; %((tp)/maxtp*511+1);
    tp1=tp; tp1(tp1<0)=0;
    tp2=tp; tp2(tp2>0)=0; tp2=abs(tp2);
    scatter(posX,posY,abs(tp),[round(tp1/max(tp1)*255) zeros(length(tp),1) round(tp2/max(tp2)*255)],'filled');
    title('Freq Tag Contra')
%     xlim([0 11])
%     ylim([0 7])
     xlim([0 elecC + 2])
     ylim([0 elecR + 2])

    % end
    
    subplot(4,2,5)
    % titles={sprintf('IM=%g (contra)',LeftF0(1)+LeftF0(2)),sprintf('IM=%g (ipsi)',RightF0(1)+RightF0(2))};
    % for nfreq=1:2
    [h, pV, ~, stats]=ttest(squeeze(pow_IM(1:sizeH,:,1))',0);
    % tp=((tp)/maxtp*511+1);
    tp=stats.tstat'/2*255; %((tp)/maxtp*511+1);
    tp1=tp; tp1(tp1<0)=0;
    tp2=tp; tp2(tp2>0)=0; tp2=abs(tp2);
    scatter(posX,posY,abs(tp),[round(tp1/max(tp1)*255) zeros(length(tp),1) round(tp2/max(tp2)*255)],'filled');
    title('IM Contra')
%     xlim([0 11])
%     ylim([0 7])
     xlim([0 elecC + 2])
     ylim([0 elecR + 2])
    % end
    
%     posX=repmat(1:10,1,5);
%     posY=[ones(1,10) 2*ones(1,10) 3*ones(1,10) 4*ones(1,10) 5*ones(1,10)];
    
    subplot(4,2,4)
    % titles={sprintf('f0=%g (contra)',LeftF0(1)),sprintf('f0=%g (contra)',LeftF0(2)),sprintf('f0=%g (ipsi)',RightF0(1)),sprintf('f0=%g (ipsi)',RightF0(2))};
    % for nfreq=1:4
    [h, pV, ~, stats]=ttest(squeeze(pow_tag(1:sizeH,:,3:4)),0,'dim',2);
    tp=mean(stats.tstat,3)/3.5*255; %((tp)/maxtp*511+1);
    tp1=tp; tp1(tp1<0)=0;
    tp2=tp; tp2(tp2>0)=0; tp2=abs(tp2);
    scatter(posX,posY,abs(tp),[round(tp1/max(tp1)*255) zeros(length(tp),1) round(tp2/max(tp2)*255)],'filled');
    title('Freq Tag Ipsi')
%     xlim([0 11])
%     ylim([0 6])
     xlim([0 elecC + 2])
     ylim([0 elecR + 2])
    % end
    
    subplot(4,2,6)
    % titles={sprintf('IM=%g (contra)',LeftF0(1)+LeftF0(2)),sprintf('IM=%g (ipsi)',RightF0(1)+RightF0(2))};
    % for nfreq=1:2
    [h, pV, ~, stats]=ttest(squeeze(pow_IM(1:sizeH,:,2))',0);
    % tp=((tp)/maxtp*511+1);
    tp=stats.tstat'/2*255; %((tp)/maxtp*511+1);
    tp1=tp; tp1(tp1<0)=0;
    tp2=tp; tp2(tp2>0)=0; tp2=abs(tp2);
    scatter(posX,posY,abs(tp),[round(tp1/max(tp1)*255) zeros(length(tp),1) round(tp2/max(tp2)*255)],'filled');
    title('IM Ipsi')
%     xlim([0 11])
%     ylim([0 6])
     xlim([0 elecC + 2])
     ylim([0 elecR + 2])
    % end
    %%
    %myChans=[find(reref_mat(:,2)==S1channels(1) & reref_mat(:,3)==0) find(reref_mat(:,2)==S1channels(2) & reref_mat(:,3)==0)];
    myChans=find(ismember(reref_mat(:,2),S1channels) & reref_mat(:,3)==0);
    subplot(4,2,7)
    simpleBarPlot(1-0.2,squeeze(mean(mean(pow_tag(myChans,:,1),3),1)),'b',0.38,'r',[],2);
    simpleBarPlot(1+0.2,squeeze(mean(mean(pow_tag(myChans,:,2),3),1)),'b',0.38,'r',[],2);
    
    simpleBarPlot(2-0.2,squeeze(mean(mean(pow_tag(myChans,:,3),3),1)),'k',0.38,'r',[],2);
    simpleBarPlot(2+0.2,squeeze(mean(mean(pow_tag(myChans,:,4),3),1)),'k',0.38,'r',[],2);
    
    
    simpleBarPlot(3-0.2,squeeze(mean(pow_IM(myChans,:,1),1)),'b',0.38,'r',[],2);
    simpleBarPlot(3+0.2,squeeze(mean(pow_IM(myChans,:,2),1)),'k',0.38,'r',[],2);
    
    set(gca,'FontSize',18,'FontWeight','bold')
    xlim([0.5 3.5])
    set(gca,'XTick',1:3,'XTickLabel',{'F0-contra','F0-ipsi','IM'})
    ylabel(sprintf('SNR (ch %g - %g)',S1channels(1),S1channels(end)))
    
% %     myChans=[find(reref_mat(:,2)==S1channels(1) & reref_mat(:,3)==1) find(reref_mat(:,2)==S1channels(2) & reref_mat(:,3)==1)];
    myChans=find(ismember(reref_mat(:,2),S1channels) & reref_mat(:,3)==1);
    subplot(4,2,8)
    simpleBarPlot(1-0.2,squeeze(mean(mean(pow_tag(myChans,:,1),3),1)),'b',0.38,'r',[],2);
    simpleBarPlot(1+0.2,squeeze(mean(mean(pow_tag(myChans,:,2),3),1)),'b',0.38,'r',[],2);
    
    simpleBarPlot(2-0.2,squeeze(mean(mean(pow_tag(myChans,:,3),3),1)),'k',0.38,'r',[],2);
    simpleBarPlot(2+0.2,squeeze(mean(mean(pow_tag(myChans,:,4),3),1)),'k',0.38,'r',[],2);
    
    
    simpleBarPlot(3-0.2,squeeze(mean(pow_IM(myChans,:,1),1)),'b',0.38,'r',[],2);
    simpleBarPlot(3+0.2,squeeze(mean(pow_IM(myChans,:,2),1)),'k',0.38,'r',[],2);
    
    set(gca,'FontSize',18,'FontWeight','bold')
    xlim([0.5 3.5])
    set(gca,'XTick',1:3,'XTickLabel',{'F0-contra','F0-ipsi','IM'})
    ylabel(sprintf('SNR (ch %g - %g)',S1channels(1),S1channels(end)))
end %for nversion