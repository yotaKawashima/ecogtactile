%% Make 2 freq sweep
freqCarrier=300;
duration=4;
samplingRate=44100;
ModRate=1;
snd1=[];
snd2=[];
snd1b=[];
snd2b=[];
for f=3:45
    beep = sin(2*pi*freqCarrier*(0:duration*samplingRate)/samplingRate);
    amplMod1 = (1+ModRate*sin(2*pi*f*(0:duration*samplingRate)/samplingRate));
%     amplMod1(amplMod1<1)=0;
%     amplMod1(amplMod1>=1)=1;
    amplMod2 = (1+ModRate*sin(2*pi*(f+5-2*rem(f,2)-1)*(0:duration*samplingRate)/samplingRate));
%     amplMod2(amplMod2<1)=0;
%     amplMod2(amplMod2>=1)=1;
    
    snd1=[snd1 beep.*amplMod1];
    snd2=[snd2 beep.*amplMod2];
    
     amplMod1 = (1+ModRate*sin(2*pi*f*(0:duration*samplingRate)/samplingRate));
    amplMod1(amplMod1<1)=0;
    amplMod1(amplMod1>=1)=1;
    amplMod2 = (1+ModRate*sin(2*pi*(f+5-2*rem(f,2)-1)*(0:duration*samplingRate)/samplingRate));
    amplMod2(amplMod2<1)=0;
    amplMod2(amplMod2>=1)=1;
    
    snd1b=[snd1b beep.*amplMod1];
    snd2b=[snd2b beep.*amplMod2];
end

%% FFT
param=[];
param.method='fft';
param.mindist=1;
for nW=1:length(3:45)
    data=abs(hilbert(snd1((nW-1)*4*samplingRate+1:(nW)*4*samplingRate)));
    [logSNR, faxis, ~]=get_logSNR(data,samplingRate,param);
    snd1_logSNR(nW,:)=logSNR;
    
    data=abs(hilbert(snd2((nW-1)*4*samplingRate+1:(nW)*4*samplingRate)));
    [logSNR, faxis, ~]=get_logSNR(data,samplingRate,param);
    snd2_logSNR(nW,:)=logSNR;
    
    data=abs(hilbert(snd1((nW-1)*4*samplingRate+1:(nW)*4*samplingRate)))+abs(hilbert(snd2((nW-1)*4*samplingRate+1:(nW)*4*samplingRate)));
    [logSNR, faxis, ~]=get_logSNR(data,samplingRate,param);
    lin_logSNR(nW,:)=logSNR;
    
    data=abs(hilbert(snd1((nW-1)*4*samplingRate+1:(nW)*4*samplingRate))).*abs(hilbert(snd2((nW-1)*4*samplingRate+1:(nW)*4*samplingRate)));
    [logSNR, ~, ~]=get_logSNR(data,samplingRate,param);
    nonlin_logSNR(nW,:)=logSNR;
    
    data=abs(hilbert(snd1b((nW-1)*4*samplingRate+1:(nW)*4*samplingRate)))+abs(hilbert(snd2b((nW-1)*4*samplingRate+1:(nW)*4*samplingRate)));
    [logSNR, faxis, ~]=get_logSNR(data,samplingRate,param);
    lin_logSNRb(nW,:)=logSNR;
    
    data=abs(hilbert(snd1b((nW-1)*4*samplingRate+1:(nW)*4*samplingRate))).*abs(hilbert(snd2b((nW-1)*4*samplingRate+1:(nW)*4*samplingRate)));
    [logSNR, ~, ~]=get_logSNR(data,samplingRate,param);
    nonlin_logSNRb(nW,:)=logSNR;
end

%% Plot
figure;
subplot(1,2,1); 
imagesc(3:45,faxis,lin_logSNR'); ylim([2 55]); set(gca,'YDIR','normal');
xlabel('Sweep Window'); ylabel('Freq (Hz)');
h = colorbar; colormap('hot'); caxis([0 70])
ylabel(h, 'log SNR')
format_fig; format_fig; title('Linear - Sine mod')

subplot(1,2,2); 
imagesc(3:45,faxis,nonlin_logSNR'); ylim([2 55]); set(gca,'YDIR','normal');
xlabel('Sweep Window'); ylabel('Freq (Hz)');
h = colorbar; colormap('hot'); caxis([0 70])
ylabel(h, 'log SNR')
format_fig; title('Non-Linear - Sine mod')

% subplot(2,2,3); 
% imagesc(3:45,faxis,lin_logSNRb'); ylim([2 55]); set(gca,'YDIR','normal');
% xlabel('Sweep Window'); ylabel('Freq (Hz)');
% h = colorbar; colormap('hot'); caxis([0 70])
% ylabel(h, 'log SNR')
% format_fig; title('Linear - Square mod')
% 
% subplot(2,2,4); 
% imagesc(3:45,faxis,nonlin_logSNRb'); ylim([2 55]); set(gca,'YDIR','normal');
% xlabel('Sweep Window'); ylabel('Freq (Hz)');
% h = colorbar; colormap('hot'); caxis([0 70])
% ylabel(h, 'log SNR')
% format_fig;  title('Non-Linear - Square mod')
