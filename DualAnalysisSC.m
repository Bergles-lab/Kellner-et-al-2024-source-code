%Dual color two photon imaging analysis comparing neurons and astrocytes
%in the SC as shown in Figure 2
%Copyright - Vered Kellner April 2024
%load('C:\Users\vered\Documents\MATLAB\dualStructSC_new.mat')
close all
clearvars -except dualstruct 
f=1;
%% load files 
[fn,dname1] = uigetfile('D:\Data\2P\In vivo\SC\*.czi;*.lsm;*.tif','Open gcamp file');
[pathstr, name, ext] = fileparts([dname1 fn]);
bf = bfopen([dname1 fn]);
[m,n] = size(bf{1}{1});
t = size(bf{1},1);
img = zeros(m,n,t,'int16');
for i=1:t
    img(:,:,i) = bf{1}{i};
end
clear bf;

[fn,dname1] = uigetfile('D:\Data\2P\In vivo\SC\*.czi;*.lsm;*.tif','Open RGECO file');
[pathstr, name, ext] = fileparts([dname1 fn]);
bf = bfopen([dname1 fn]);
[m,n] = size(bf{1}{1});
t = size(bf{1},1);
imgRG = zeros(m,n,t,'int16');
for i=1:t
    imgRG(:,:,i) = bf{1}{i};
end
clear bf;

imgCrop=img(13:end,13:end,:);
imgRGCrop=imgRG(13:end,13:end,:);

imgdA=imgCrop;
imgdN=imgRGCrop;

%% spatial analysis
figure('Position',[10,20,1600,700]);  subplot(3,1,1); imagesc(squeeze(mean(imgdN,2)));%caxis([2 15]);
colorbar %axis off; 
title('Neurons')
subplot(3,1,2); imagesc(squeeze(mean(imgdA,2)));%caxis([2 15]);
colorbar %axis off; 
title('Astrocytes')
subplot(3,1,3); plot(sum(squeeze(mean(imgdN))),'k'); xlim([0 size(imgdN,3)]);
hold on; plot(sum(squeeze(mean(imgdA))),'r'); xlim([0 size(imgdA,3)]); 
xlabel('Frames'); ylabel('Intensity'); legend('Neurons','Astrocytes')

%% find peaks over entire image
t=size(imgdA,3);
astTrace=squeeze(mean(mean(imgdA)));
astTraceF=smooth(msbackadj([1:t]',astTrace,'WindowSize',45,'StepSize',45,'Showplot',0));
tN=size(imgdN,3);
nTrace=squeeze(mean(mean(imgdN)));
nTraceF=smooth(msbackadj([1:tN]',nTrace,'WindowSize',45,'StepSize',45,'Showplot',0));

areaInp=inputdlg({'Area?'});
if strcmp(areaInp{:},'SC')>0
    %for SC
    threshAst=median(astTraceF)+mad(astTraceF)/1.5;
    figure('Position',[10,500,1600,350]); findpeaks(astTraceF,'MinPeakHeight',threshAst,'MinPeakProminence',threshAst/2,...
        'MinPeakWidth',5,'MinPeakDistance',20,'WidthReference',...
        'halfheight','Annotate','extents');
    [astPk,astLoc,wAst]=findpeaks(astTraceF,'MinPeakHeight',threshAst,'MinPeakProminence',threshAst/2,...
        'MinPeakWidth',5,'MinPeakDistance',20,'WidthReference',...
        'halfheight','Annotate','extents');    
    %for SC
    threshN=median(nTraceF)+1*mad(nTraceF)/1.5;
    figure('Position',[10,600,1600,350]); findpeaks(nTraceF,'MinPeakHeight',threshN,'MinPeakProminence', threshN/2,...
        'MinPeakWidth',5,'MinPeakDistance',10,'WidthReference',...
        'halfheight','Annotate','extents');
    [nPk,nLoc,wN,pN]=findpeaks(nTraceF,'MinPeakHeight',threshN,'MinPeakProminence', threshN/2,...
        'MinPeakWidth',5,'MinPeakDistance',10,'WidthReference',...
        'halfheight','Annotate','extents');
else
    %for IC
    threshAst=median(astTraceF);
    figure('Position',[10,500,1600,350]);findpeaks(astTraceF,'MinPeakHeight',threshAst,'MinPeakProminence',threshAst/2,...
        'MinPeakWidth',1,'MinPeakDistance',5,'WidthReference',...
        'halfheight','Annotate','extents');
    [astPk,astLoc,wAst]=findpeaks(astTraceF,'MinPeakHeight',threshAst,'MinPeakProminence',threshAst/2,...
        'MinPeakWidth',1,'MinPeakDistance',5,'WidthReference',...
        'halfheight','Annotate','extents');    
    %for IC
    threshN=median(nTraceF);
    figure('Position',[10,500,1600,350]); findpeaks(nTraceF,'MinPeakHeight',threshN,'MinPeakProminence', threshN/2,...
        'MinPeakDistance',1,'WidthReference',...
        'halfheight','Annotate','extents');
    [nPk,nLoc,wN,pN]=findpeaks(nTraceF,'MinPeakHeight',threshN,'MinPeakProminence', threshN/2,...
        'MinPeakDistance',1,'WidthReference',...
        'halfheight','Annotate','extents');
end

%% Find astrocyte signal amplitude 3 frames after neuronal peak
if nLoc(end)+3<=length(astTraceF)
    aAmp=astTraceF(nLoc+3);
    for b=1:100
        inds=randi(t,length(nLoc),1);
        aAmpSim(:,b)=astTraceF(inds);
    end
    figure; scatter(nPk,aAmp,'filled'); hold on; scatter(nPk,mean(aAmpSim,2),'r','MarkerFaceColor','r')
else
    aAmp=astTraceF(nLoc(1:end-1)+3);
    for b=1:100
        inds=randi(t,length(nLoc(1:end-1)),1);
        aAmpSim(:,b)=astTraceF(inds);
    end
    figure; scatter(nPk(1:end-1),aAmp,'filled'); hold on; scatter(nPk(1:end-1),mean(aAmpSim,2),'r','MarkerFaceColor','r')
end

%xlim([0 1]); 
xlabel('Neuronal event amplitude (A.U.)');
ylabel('Astrocyte amplitude 1.5s after neuronal event (A.U.)')
legend('Real data','Randomized data','Location','Northwest')
[R, P]=corrcoef(nPk,aAmp)

%% plot trace
figure('Position',[10,100,1600,700]); 
plot(nTraceF,'k'); hold on; plot(astTraceF,'r')
 xlim([0 length(astTraceF)]); 
ylabel('Fluorescence (A.U.)'); xlabel('Frames')
legend('Neurons','Astrocytes')

figure('Position',[10,100,1600,700]);  plot(nTraceF./max(nTraceF),'k'); hold on; plot(astTraceF./max(astTraceF),'r')
 xlim([0 length(astTraceF)]); 
ylabel('Normalized Fluorescence (A.U.)'); xlabel('Frames')
legend('Neurons','Astrocytes')


figure('Position',[10,100,1600,700]); 
plot(nTraceF./max(nTraceF),'b'); hold on; plot(nLoc,nPk./max(nTraceF),'ok','MarkerFaceColor','k','MarkerSize',7)
plot(astTraceF./max(astTraceF),'r'); plot(astLoc,astPk./max(astTraceF),'or','MarkerFaceColor','r','MarkerSize',7)
 xlim([0 length(astTraceF)]); 
ylabel('Normalized Fluorescence (A.U.)'); xlabel('Frames')
legend('Neurons','Astrocytes')

%% compute the correlation between the astrocyte and neuron signals
rA2N = corr(astTraceF,nTraceF)

%% cross correlations (this is on the whole signal)
[aNcor,aNlag] = xcorr(astTraceF-mean(astTraceF),nTraceF-mean(nTraceF),'normalized'); %need to use zeros here for the movement related events - doesn't work with nans
sampRate=2;
[val,I] = max(abs(aNcor));
lagDiff = aNlag(I)
timeDiff = lagDiff/sampRate

figure
plot(aNlag,aNcor)
xlim([-length(astTraceF) length(astTraceF)])
xlabel('Frames'); ylabel('Normalized correlation')

%% Now gather all the data based on events - astrocytes
lenAtrace=round(max(wAst))+1;
eventATrace=nan(length(astPk),lenAtrace*2+7);
eventATmTrace=nan(length(astPk),lenAtrace*2+7);
 for aa=1:length(astPk)   
    eventAinfo(aa).pkFrm=astLoc(aa);
    eventAinfo(aa).pkAmp=astPk(aa);
    eventAinfo(aa).pkHW=wAst(aa)/sampRate; %hw in seconds
    if astLoc(aa)+round(wAst(aa))+6>length(astTraceF)
        tempInd=astLoc(aa)-round(wAst(aa)):length(astTraceF);
    elseif astLoc(aa)-round(wAst(aa))<=0
        tempInd=1:astLoc(aa)+round(wAst(aa))+6;
    else
        tempInd=astLoc(aa)-round(wAst(aa)):astLoc(aa)+round(wAst(aa))+6;
    end
    startInd=30-find(astLoc(aa)==tempInd);
    if startInd>0        
    eventATrace(aa,startInd:startInd+length(tempInd)-1)=astTraceF(tempInd);  
    else
        eventATrace(aa,1:1+length(tempInd)-1)=astTraceF(tempInd);
    end
    eventATmTrace(aa,1:length(tempInd))=tempInd;
    tempInd=[];
 end
figure; plot(eventATrace')
title('Astrocyte events')

%% Now gather all the data based on events - neurons
lenNtrace=round(max(wN))+1;
eventNTrace=nan(length(nPk),lenNtrace*2+7);
eventNTmTrace=nan(length(nPk),lenNtrace*2+7);
 for aa=1:length(nPk)   
    eventNinfo(aa).pkFrm=nLoc(aa);
    eventNinfo(aa).pkAmp=nPk(aa);
    eventNinfo(aa).pkHW=wN(aa)/sampRate; %hw in seconds
    if nLoc(aa)+round(wN(aa))+6 > length(nTraceF)
        tempIndN=nLoc(aa)-round(wN(aa)):length(nTraceF);
    elseif nLoc(aa)-round(wN(aa))<=0
        tempIndN=1:nLoc(aa)+round(wN(aa))+6;
    else
        tempIndN=nLoc(aa)-round(wN(aa)):nLoc(aa)+round(wN(aa))+6;
    end
    startInd=30-find(nLoc(aa)==tempIndN);
    if startInd>0
    eventNTrace(aa,startInd:startInd+length(tempIndN)-1)=nTraceF(tempIndN);  
    else
        eventNTrace(aa,1:1+length(tempIndN)-1)=nTraceF(tempIndN);
    end
    eventNTmTrace(aa,1:length(tempIndN))=tempIndN;
    tempIndN=[];
 end
 figure; plot(eventNTrace')
title('Neuron events')

%% put into struct
dualstruct(f).Name=name;
dualstruct(f).Area=areaInp{:};
ageList={'6','7','8','9','10','11','12'};
ageAns=listdlg('ListString',ageList,'PromptString','Age (days):');
dualstruct(f).Age=ageList{ageAns};
sexList={'Male','Female','?'};
sexAns=listdlg('ListString',sexList,'PromptString','Sex:');
dualstruct(f).Sex=sexList{sexAns};
dualstruct(f).Pks.AstStats=eventAinfo;
dualstruct(f).Pks.NeurStats=eventNinfo;
dualstruct(f).Traces.eventATrace=eventATrace;
dualstruct(f).Traces.eventATraceTm=eventATmTrace;
dualstruct(f).Traces.eventNTrace=eventNTrace;
dualstruct(f).Traces.eventNTraceTm=eventNTmTrace;
dualstruct(f).Corr.AllFrmAstAmp=aAmp; %3 frames after neuron peak
dualstruct(f).Corr.AllFrmNAmp=nPk;
dualstruct(f).Corr.aAmpSim=aAmpSim;
dualstruct(f).Corr.CorrVal=rA2N;
dualstruct(f).Corr.Timelag=timeDiff; %the lag in seconds
dualstruct(f).Corr.CrossCorr=[aNlag',aNcor]; %need to transform aNlag to plot
dualstruct(f).Frm=size(img,3);
