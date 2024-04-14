%% Group 2P analysis

clearvars -except twoPstruct; 
close all

twoPstruct=twoPROIs;
%% get info
countASC=0;
countAIC=0;
ampASCAll=cell(1,1);
hwASCAll=cell(1,1);
countNSC=0;
countNIC=0;
ampNSCAll=cell(1,1);
hwNSCAll=cell(1,1);

for a=1:size(twoPstruct,2)
    switch twoPstruct(a).Cell
        case 'Astrocyte'
            if strcmp(twoPstruct(a).Area,'SC')
                countASC=countASC+1;                    
                ampASC(countASC)=median([twoPstruct(a).eventInfo(:,1)]);
                hwASC(countASC)=median([twoPstruct(a).eventInfo(:,3)]);
                ampASCAll{countASC}=([twoPstruct(a).eventInfo(:,1)]);
                hwASCAll{countASC}=([twoPstruct(a).eventInfo(:,3)]);
                tempfreq=[twoPstruct(a).eventInfo(:,4)]; 
                freqASC(countASC)=median(accumarray(tempfreq,1))./(twoPstruct(a).frmNum/twoPstruct(a).frmRate);
            else
                countAIC=countAIC+1;
                ampAIC(countAIC)=median([twoPstruct(a).eventInfo(:,1)]);
                hwAIC(countAIC)=median([twoPstruct(a).eventInfo(:,3)]);
                ampAICAll{countAIC}=([twoPstruct(a).eventInfo(:,1)]);
                hwAICAll{countAIC}=([twoPstruct(a).eventInfo(:,3)]);
                tempfreq=[twoPstruct(a).eventInfo(:,4)]; 
                freqAIC(countAIC)=median(accumarray(tempfreq,1))./(twoPstruct(a).frmNum/twoPstruct(a).frmRate);
            end
    end
    tempfreq=[];
end

%% plot frequency
lng=[length(freqASC),length(freqAIC),length(freqNSC),length(freqNIC)];
eventfreqCellMat=nan(max(lng),4);
eventfreqCellMat(1:length(freqASC),1)=freqASC*60;
eventfreqCellMat(1:length(freqAIC),2)=freqAIC*60;
eventfreqCellMat(1:length(freqNSC),3)=freqNSC*60;
eventfreqCellMat(1:length(freqNIC),4)=freqNIC*60;
figure; h=plotSpread(eventfreqCellMat,'distributionColors',[0.5,0.5,0.5]);
hold on;  
er = errorbar([1:4],nanmean(eventfreqCellMat),nanstd(eventfreqCellMat)./sqrt(lng),'ok');    
er.MarkerFaceColor = [0 0 0];                            
er.MarkerSize = 12;  
xlim([0.5 4.5])
ylim([0 9])
set(gca,'XTickLabel',{'A SC','A IC','N SC','N IC'})
ylabel('Frequency (events/min)')

lng=[length(freqASC),length(freqAIC)];
eventfreqCellMat=nan(max(lng),2);
eventfreqCellMat(1:length(freqASC),1)=freqASC*60;
eventfreqCellMat(1:length(freqAIC),2)=freqAIC*60;

figure; h=plotSpread(eventfreqCellMat,'distributionColors',[0.5,0.5,0.5]);
hold on;  
er = errorbar([1:2],nanmean(eventfreqCellMat),nanstd(eventfreqCellMat)./sqrt(lng),'ok');    
er.MarkerFaceColor = [0 0 0];                            
er.MarkerSize = 12;  
xlim([0.5 2.5])
ylim([0 1.2])
set(gca,'XTickLabel',{'A SC','A IC'})
ylabel('Frequency (events/min)')

%% Get data for average event
countNSC=0;
countNIC=0;
eventAllSC=[];
eventAllIC=[];
countASC=0;
countAIC=0;
eventAllASC=[];
eventAllAIC=[];

for a=1:size(twoPstruct,2)
    switch twoPstruct(a).Cell
        case 'Neuron'
            if strcmp(twoPstruct(a).Area,'SC')
                countNSC=countNSC+1;
                avgEventNSC(countNSC,:)=nanmean(twoPstruct(a).Neuropil.eventNTrace);
                stdEventNSC(countNSC,:)=nanstd(twoPstruct(a).Neuropil.eventNTrace);
                numEventNSC(countNSC)=size(twoPstruct(a).Neuropil.eventNTrace,1);
                eventAllSC=[eventAllSC;(twoPstruct(a).Neuropil.eventNTrace)];
                %         eventTmNCtl(a,:)=nanmean(twoPstruct(a).eventNTmNew);
            elseif strcmp(twoPstruct(a).Area,'IC')
                countNIC=countNIC+1;
                avgEventNIC(countNIC,:)=nanmean(twoPstruct(a).Neuropil.eventNTrace);
                stdEventNIC(countNIC,:)=nanstd(twoPstruct(a).Neuropil.eventNTrace);
                numEventNIC(countNIC)=size(twoPstruct(a).Neuropil.eventNTrace,1);
                eventAllIC=[eventAllIC;(twoPstruct(a).Neuropil.eventNTrace)];
            end
        case 'Astrocyte'
            if strcmp(twoPstruct(a).Area,'SC')
                countASC=countASC+1;
                ampASC(countASC)=median([twoPstruct(a).eventInfo(:,1)]);
                hwASC(countASC)=median([twoPstruct(a).eventInfo(:,3)]);
                avgEventASC(countASC,:)=nanmean(twoPstruct(a).eventNTrace);
                stdEventASC(countASC,:)=nanstd(twoPstruct(a).eventNTrace);
                numEventASC(countASC)=size(twoPstruct(a).eventNTrace,1);
                eventAllASC=[eventAllASC;(twoPstruct(a).eventNTrace)];
                %         eventTmNCtl(a,:)=nanmean(twoPstruct(a).eventNTmNew);
            elseif strcmp(twoPstruct(a).Area,'IC')
                countAIC=countAIC+1;
                ampASC(countASC)=median([twoPstruct(a).eventInfo(:,1)]);
                hwASC(countASC)=median([twoPstruct(a).eventInfo(:,3)]);
                avgEventAIC(countAIC,:)=nanmean(twoPstruct(a).eventNTrace);
                stdEventAIC(countAIC,:)=nanstd(twoPstruct(a).eventNTrace);
                numEventAIC(countAIC)=size(twoPstruct(a).eventNTrace,1);
                eventAllAIC=[eventAllAIC;(twoPstruct(a).eventNTrace)];
            end
    end
end

%% plot average event for SC &IC same laser power

countASC=0;
countAIC=0;
eventAllASC=[];
eventAllAIC=[];

for a=1:size(twoPstruct,2)
    switch twoPstruct(a).Cell
        case 'Neuron'
            if strcmp(twoPstruct(a).Area,'SC')
                countNSC=countNSC+1;
                avgEventNSC(countNSC,:)=nanmean(twoPstruct(a).Neuropil.eventNTrace);
                stdEventNSC(countNSC,:)=nanstd(twoPstruct(a).Neuropil.eventNTrace);
                numEventNSC(countNSC)=size(twoPstruct(a).Neuropil.eventNTrace,1);
                eventAllSC=[eventAllSC;(twoPstruct(a).Neuropil.eventNTrace)];
                %         eventTmNCtl(a,:)=nanmean(twoPstruct(a).eventNTmNew);
            elseif strcmp(twoPstruct(a).Area,'IC')
                countNIC=countNIC+1;
                avgEventNIC(countNIC,:)=nanmean(twoPstruct(a).Neuropil.eventNTrace);
                stdEventNIC(countNIC,:)=nanstd(twoPstruct(a).Neuropil.eventNTrace);
                numEventNIC(countNIC)=size(twoPstruct(a).Neuropil.eventNTrace,1);
                eventAllIC=[eventAllIC;(twoPstruct(a).Neuropil.eventNTrace)];
            end
        case 'Astrocyte'
            if strcmp(twoPstruct(a).Area,'SC')
                countASC=countASC+1;
                ampASC(countASC)=median([twoPstruct(a).eventInfo(:,1)]./max([twoPstruct(a).eventInfo(:,1)]));
                hwASC(countASC)=median([twoPstruct(a).eventInfo(:,3)]);
                tempfreq=[twoPstruct(a).eventInfo(:,4)];                
                freqASC(countASC)=median(accumarray(tempfreq,1))./(twoPstruct(a).frmNum/twoPstruct(a).frmRate);
                avgEventASC(countASC,:)=nanmean(twoPstruct(a).eventNTrace);
                stdEventASC(countASC,:)=nanstd(twoPstruct(a).eventNTrace);
                numEventASC(countASC)=size(twoPstruct(a).eventNTrace,1);
                eventAllASC=[eventAllASC;(twoPstruct(a).eventNTrace)];
                intTempSC=cellfun(@trapz,cellfun(@transpose,twoPstruct(a).eventNTraceAll,...
                    'UniformOutput',false),'UniformOutput',false);
                integralSC(countASC)=nanmedian(cellfun(@nanmedian,intTempSC));
                %         eventTmNCtl(a,:)=nanmean(twoPstruct(a).eventNTmNew);
            elseif strcmp(twoPstruct(a).Area,'IC')
                countAIC=countAIC+1;
                ampAIC(countAIC)=median([twoPstruct(a).eventInfo(:,1)]./max([twoPstruct(a).eventInfo(:,1)]));
                hwAIC(countAIC)=median([twoPstruct(a).eventInfo(:,3)]);
                tempfreq=[twoPstruct(a).eventInfo(:,4)];                
                freqAIC(countAIC)=median(accumarray(tempfreq,1))./(twoPstruct(a).frmNum/twoPstruct(a).frmRate);
                avgEventAIC(countAIC,:)=nanmean(twoPstruct(a).eventNTrace);
                stdEventAIC(countAIC,:)=nanstd(twoPstruct(a).eventNTrace);
                numEventAIC(countAIC)=size(twoPstruct(a).eventNTrace,1);
                eventAllAIC=[eventAllAIC;(twoPstruct(a).eventNTrace)];
                intTempIC=cellfun(@trapz,cellfun(@transpose,twoPstruct(a).eventNTraceAll,...
                    'UniformOutput',false),'UniformOutput',false);
                integralIC(countAIC)=nanmedian(cellfun(@nanmedian,intTempIC));
            end
    end
end
%% plot frequency
lng=[length(freqASC),length(freqAIC)];
eventfreqCellMat=nan(max(lng),2);
eventfreqCellMat(1:length(freqASC),1)=freqASC*60;
eventfreqCellMat(1:length(freqAIC),2)=freqAIC*60;

%figure; h=plotSpread(eventfreqCellMat,'distributionColors',[0.5,0.5,0.5]);
figure; plot([0.5:1.5],eventfreqCellMat,'-o','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',10)
hold on;  %bar([1:2],nanmean(eventDurCellMat),0.25,'FaceColor','none')
er = errorbar([0.5:1.5],nanmean(eventfreqCellMat),nanstd(eventfreqCellMat)./sqrt(lng),'-ok');    
er.MarkerFaceColor = [0 0 0];                            
er.MarkerSize = 12;  
% plot([1:2],nanmedian(eventDurCellMat),'dk','MarkerFaceColor','k','MarkerSize',10)
xlim([0 2])
ylim([0 1.2])
set(gca,'XTickLabel',{'','A SC','','A IC',''})
ylabel('Frequency (events/min)')

h=lillietest(eventfreqCellMat(:,1))
h=lillietest(eventfreqCellMat(:,2))
% [h,p,ci,stats] = ttest(eventfreqCellMat(:,1),eventfreqCellMat(:,2))
[p,h,stats] = signrank(eventfreqCellMat(:,1),eventfreqCellMat(:,2))


%plot integral
lng=[length(integralSC),length(integralIC)];
eventIntCellMat=nan(max(lng),2);
eventIntCellMat(1:length(integralSC),1)=integralSC;
eventIntCellMat(1:length(integralIC),2)=integralIC;

figure; plot([0.5:1.5],eventIntCellMat,'-o','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],...
    'MarkerEdgeColor',[0.5 0.5 0.5],'MarkerSize',10)
hold on;  
er = errorbar([0.5:1.5],nanmean(eventIntCellMat),nanstd(eventIntCellMat)./sqrt(lng),'-ok');    
er.MarkerFaceColor = [0 0 0];                            
er.MarkerSize = 12;  
xlim([0 2])
ylim([0 500])
set(gca,'XTickLabel',{'','A SC','','A IC',''})
ylabel('Integral of events')

h=lillietest(eventIntCellMat(:,1))
h=lillietest(eventIntCellMat(:,2))
%paired ttest
% [h,p,ci,stats] = ttest(eventIntCellMat(:,1),eventIntCellMat(:,2))
[p,h,stats] = signrank(eventIntCellMat(:,1),eventIntCellMat(:,2))

%% Plot average event
eventAllNSCSh=eventAllSC(:,10:35);%this is chosen based on the signal
eventAllNSCSh(isnan(eventAllNSCSh))=0;

eventAllNICSh=eventAllIC(:,10:35);
eventAllNICSh(isnan(eventAllNICSh))=0;

eventAllASCSh=eventAllASC(:,10:35);%this is chosen based on the signal
eventAllASCSh(isnan(eventAllASCSh))=0;

eventAllAICSh=eventAllAIC(:,10:35);
eventAllAICSh(isnan(eventAllAICSh))=0;

eventTime=[-20:20];
frameRate=2;
%bring down to 0
normEventNSC=nanmean(avgEventNSC,1);
normEventNSC=normEventNSC-normEventNSC(1);
normEventNIC=nanmean(avgEventNIC,1);
normEventNIC=normEventNIC-normEventNIC(1);

figure; h1=shadedErrorBar(eventTime/frameRate,(normEventNSC),nanmean(stdEventNSC,1)./sqrt(nanmean(numEventNSC)),'lineprops','-k','transparent',1);
hold on; h2=shadedErrorBar(eventTime/frameRate,normEventNIC,nanmean(stdEventNIC,1)./sqrt(nanmean(numEventNIC)),'lineprops','-r','transparent',1);
xlim([-10 10])
% ylim([-0.02 0.20])
legend([h1.mainLine,h2.mainLine],'SC','IC')
xlabel('Time (s)')
ylabel('Amplitude (F.U.)')
title('Neuron')

figure; h1=shadedErrorBar(eventTime/frameRate,(normEventNSC)./max(normEventNSC),nanmean(stdEventNSC,1)./max(normEventNSC)/sqrt(nanmean(numEventNSC)),'lineprops','-k','transparent',1);
hold on; h2=shadedErrorBar(eventTime/frameRate,normEventNIC./max(normEventNIC),nanmean(stdEventNIC,1)./max(normEventNIC)/sqrt(nanmean(numEventNIC)),'lineprops','-r','transparent',1);
xlim([-10 10])
ylim([-0.2 1.2])
legend([h1.mainLine,h2.mainLine],'SC','IC')
xlabel('Time (s)')
ylabel('Normalized to maximal amplitude')
title('Neuron')

%bring down to 0
normEventASC=nanmean(avgEventASC,1);
normEventASC=normEventASC-normEventASC(1);
normEventAIC=nanmean(avgEventAIC,1);
normEventAIC=normEventAIC-normEventAIC(1);

figure; h1=shadedErrorBar(eventTime/frameRate,(normEventASC),nanmean(stdEventASC,1)./sqrt(nanmean(numEventASC)),'lineprops','-k','transparent',1);
hold on; h2=shadedErrorBar(eventTime/frameRate,normEventAIC,nanmean(stdEventAIC,1)./sqrt(nanmean(numEventAIC)),'lineprops','-r','transparent',1);
xlim([-10 10])
% ylim([-0.02 0.20])
legend([h1.mainLine,h2.mainLine],'SC','IC')
xlabel('Time (s)')
ylabel('Amplitude (F.U.)')
title('Astrocyte')

figure; h1=shadedErrorBar(eventTime/frameRate,(normEventASC)./max(normEventASC),nanmean(stdEventASC,1)./max(normEventASC)/sqrt(nanmean(numEventASC)),'lineprops','-k','transparent',1);
hold on; h2=shadedErrorBar(eventTime/frameRate,normEventAIC./max(normEventAIC),nanmean(stdEventAIC,1)./max(normEventAIC)/sqrt(nanmean(numEventAIC)),'lineprops','-r','transparent',1);
xlim([-10 10])
ylim([-0.2 1.4])
legend([h1.mainLine,h2.mainLine],'SC','IC')
xlabel('Time (s)')
ylabel('Normalized to maximal amplitude')
title('Astrocyte')
