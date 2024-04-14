%% Compare over all mice
%gets data from dataStruct generated in the tracking_framework script
%Copyright - Vered Kellner April 2024
%load('SC_tracking_NoGCaMP3_new.mat')
clearvars -except dataStruct
close all

%% interesting parameters to extract
%determine frequency of events in left/right SC
%interval between events?
durEventsAll=[];
distTravelAll=[];
dirAll=[];
biEvent=zeros(1,size(dataStruct,2));
for a=1:size(dataStruct,2)
    countR=0;
    countL=0;
    rightInd=[];
    leftInd=[];
    frmActiveR=[];
    frmActiveL=[];
    if ~isempty(dataStruct(a).Events.WaveInfo)
        frmMatR=nan(dataStruct(a).Events.Frm,size(dataStruct(a).Events.WaveInfo,1));
        frmMatL=nan(dataStruct(a).Events.Frm,size(dataStruct(a).Events.WaveInfo,1));
        for i=1:size(dataStruct(a).Events.WaveInfo,1)
            if strcmp(dataStruct(a).Events.WaveInfo{i}.Side,'Right')
                countR=countR+1;
                rightInd(countR)=i;
                frmMatR(dataStruct(a).Events.WaveInfo{i}.frms,countR)=1;
                frmActiveR=[frmActiveR;[dataStruct(a).Events.WaveInfo{i}.frms(1),dataStruct(a).Events.WaveInfo{i}.frms(end)]];
            else
                countL=countL+1;
                leftInd(countL)=i;
                frmMatL(dataStruct(a).Events.WaveInfo{i}.frms,countL)=1;
                frmActiveL=[frmActiveL;[dataStruct(a).Events.WaveInfo{i}.frms(1),dataStruct(a).Events.WaveInfo{i}.frms(end)]];
            end
            tempIndY=find(~isnan(dataStruct(a).Events.WaveInfo{i}.centroid(1,:)),1,'first');
            tempIndX=find(~isnan(dataStruct(a).Events.WaveInfo{i}.centroid(2,:)),1,'first');
            origP{a}(i,:)=[round(dataStruct(a).Events.WaveInfo{i}.centroid(1,tempIndY)),round(dataStruct(a).Events.WaveInfo{i}.centroid(2,tempIndX))]; % where do the events originate?
            tempIndY_L=find(~isnan(dataStruct(a).Events.WaveInfo{i}.centroid(1,:)),1,'last');
            tempIndX_L=find(~isnan(dataStruct(a).Events.WaveInfo{i}.centroid(2,:)),1,'last');
            endP{a}(i,:)=[round(dataStruct(a).Events.WaveInfo{i}.centroid(1,tempIndY_L)),round(dataStruct(a).Events.WaveInfo{i}.centroid(2,tempIndX_L))];
            %             %distance traveled and information about origin and direction is contained in the following:
            distTravelTemp=origP{a}(i,:)-endP{a}(i,:); %if column 1 is positive- means it went from top to bottom, if column 2 is positive - means it went from right to left
            %distance of the vector sqrt(X^2+Y^2)
            distTravel{a}(i)=sqrt(sum(distTravelTemp.^2));
            dir{a}(i)=rad2deg(atan2(distTravelTemp(2),distTravelTemp(1))); %direction of movement in degrees
            vel{a}(i)=nanmedian(dataStruct(a).Events.WaveInfo{i}.vel); %median velocity
            sdDir{a}(i)=median([nanstd(dataStruct(a).Events.WaveInfo{i}.dir(1,:)),...
                nanstd(dataStruct(a).Events.WaveInfo{i}.dir(2,:))]);
            if isfield(dataStruct(a).Events.WaveInfo{i}, 'amplitude')>0
                amp{a}(i)=dataStruct(a).Events.WaveInfo{i}.amplitude;
            else
                amp{a}=nan;
            end
        end        
        if strcmp(dataStruct(a).use,'Both')
            tempOverlapR=nansum(frmMatR,2);
            overlapR(a)=max(tempOverlapR);% how many events are active on the same side (indicative of branching) at the same time?
            tempOverlapL=nansum(frmMatL,2);
            overlapL(a)=max(tempOverlapL);
            waveFrqL(a)=length(leftInd)/(dataStruct(a).Events.Frm/10)*60;%events per minute
            waveFrqR(a)=length(rightInd)/(dataStruct(a).Events.Frm/10)*60;%events per minute
            avgDurEvents(a)=nanmedian([dataStruct(a).Events.Duration])/10; %median duration in seconds
            mdDurEvents(a)=mad([dataStruct(a).Events.Duration])/10; %mad duration in seconds
            durEventsAll=[durEventsAll;[dataStruct(a).Events.Duration/10,repmat(a,i,1)]];
            distTravelAll=[distTravelAll;[distTravel{a}',repmat(a,i,1)]];
            dirAll=[dirAll;[dir{a}',repmat(a,i,1)]];
            if ~isempty(frmActiveR)
            frmActiveR=sortrows(frmActiveR);
            intrvlR{a}=diff(frmActiveR(:,1));
            end
            if ~isempty(frmActiveL)
            frmActiveL=sortrows(frmActiveL);            
            intrvlL{a}=diff(frmActiveL(:,1));
            end   
            %check bilaterality of events
            if ~isempty(frmActiveR) && ~isempty(frmActiveL)
                lng=[size(frmActiveR,1),size(frmActiveL,1)];
                [M,indM]=max(lng);
                for c=1:M
                    if indM==1
                        if sum(ismember(frmActiveR(c,1)-10:frmActiveR(c,1)+10,frmActiveL(:,1)))>0
                            biEvent(a)=biEvent(a)+1;
                        end
                    elseif indM==2
                        if sum(ismember(frmActiveL(c,1)-10:frmActiveL(c,1)+10,frmActiveR(:,1)))>0
                            biEvent(a)=biEvent(a)+1;
                        end
                    end
                end
                biEventPrcnt(a)=biEvent(a)/M;
            end
        elseif strcmp(dataStruct(a).use,'Right')
            tempOverlapR=nansum(frmMatR,2);
            overlapR(a)=max(tempOverlapR);% how many events are active on the same side (indicative of branching) at the same time?
            waveFrqR(a)=length(rightInd)/(dataStruct(a).Events.Frm/10)*60;%events per minute
            overlapL(a)=nan;
            waveFrqL(a)=nan;
            avgDurEvents(a)=nanmedian([dataStruct(a).Events.Duration(rightInd)])/10; %median duration in seconds
            mdDurEvents(a)=mad([dataStruct(a).Events.Duration(rightInd)])/10; %mad duration in seconds
            origP{a}(leftInd,:)=[];
            endP{a}(leftInd,:)=[];
            distTravel{a}(leftInd)=[];
            dir{a}(leftInd)=[];
            durEventsAll=[durEventsAll;[dataStruct(a).Events.Duration(rightInd)/10,repmat(a,length(rightInd),1)]];
            distTravelAll=[distTravelAll;[distTravel{a}',repmat(a,length(rightInd),1)]];
            dirAll=[dirAll;[dir{a}',repmat(a,length(rightInd),1)]];
            if ~isempty(frmActiveR)
            frmActiveR=sortrows(frmActiveR);
            intrvlR{a}=diff(frmActiveR(:,1));
            end
            intrvlL{a}=nan;
            
        elseif strcmp(dataStruct(a).use,'Left')
            tempOverlapL=nansum(frmMatL,2);
            overlapL(a)=max(tempOverlapL);
            overlapR(a)=nan;
            waveFrqL(a)=length(leftInd)/(dataStruct(a).Events.Frm/10)*60;%events per minute
            waveFrqR(a)=nan;
            avgDurEvents(a)=nanmedian([dataStruct(a).Events.Duration(leftInd)])/10; %median duration in seconds
            mdDurEvents(a)=mad([dataStruct(a).Events.Duration(leftInd)])/10; %mad duration in seconds
            origP{a}(rightInd,:)=[];
            endP{a}(rightInd,:)=[];
            distTravel{a}(rightInd)=[];
            dir{a}(rightInd)=[];
            durEventsAll=[durEventsAll;[dataStruct(a).Events.Duration(leftInd)/10,repmat(a,length(leftInd),1)]];
            distTravelAll=[distTravelAll;[distTravel{a}',repmat(a,length(leftInd),1)]];
            dirAll=[dirAll;[dir{a}',repmat(a,length(leftInd),1)]];
            if ~isempty(frmActiveL)
            frmActiveL=sortrows(frmActiveL);
            intrvlL{a}=diff(frmActiveL(:,1));
            end
            intrvlR{a}=nan;
        end
    else
        overlapL(a)=nan;
        overlapR(a)=nan;
        avgDurEvents(a)=nan;
        mdDurEvents(a)=nan;
        origP{a}=nan;
        endP{a}=nan;
        distTravel{a}=nan;
        dir{a}=nan;
        intrvlR{a}=nan;
        intrvlL{a}=nan;
        amp{a}=nan;
    end
    if ~isempty(dataStruct(a).IC)
        switch dataStruct(a).IC.use
            case 'Both'
                tempL=size(dataStruct(a).IC.LIC,1)/(dataStruct(a).IC.frmNum/10)*60;%events per minute
                tempR=size(dataStruct(a).IC.RIC,1)/(dataStruct(a).IC.frmNum/10)*60;%events per minute
                ICfreqL{a}=tempL;
                ICfreqR{a}=tempR;
                ICfreq{a}=median([tempL tempR]);
            case 'Left'
                tempL=size(dataStruct(a).IC.LIC,1)/(dataStruct(a).IC.frmNum/10)*60;%events per minute                
                ICfreq{a}=size(dataStruct(a).IC.LIC,1)/(dataStruct(a).IC.frmNum/10)*60;%events per minute
                ICfreqL{a}=tempL;
            case 'Right'
                tempR=size(dataStruct(a).IC.RIC,1)/(dataStruct(a).IC.frmNum/10)*60;%events per minute
                ICfreq{a}=size(dataStruct(a).IC.RIC,1)/(dataStruct(a).IC.frmNum/10)*60;%events per minute
                ICfreqR{a}=tempR;
        end
    else
        ICfreq{a}=nan;
    end
end


%% compare neurons to astrocytes
%neurons: SNAP25-GCAMP6s ; Thy1-jRGECO1a
%Astrocytes: Aldh1l1-GCaMP6s/GCaMP6sh; 
neuronIndSN=[];
neuronIndRG=[];
astroIndALGC6=[];
astroIndALGC64HT=[];
for a=1:size(dataStruct,2)
    switch dataStruct(a).Cell
        case 'Neuron'
            if strcmp(dataStruct(a).Promoter,'SNAP25')
                if strcmp(dataStruct(a).Manipulation,'None') || strcmp(dataStruct(a).ManipDetails,'Pre') || strcmp(dataStruct(a).Manipulation,'Pre')...
                    neuronIndSN=[neuronIndSN;a];
                end
            elseif strcmp(dataStruct(a).Promoter,'Thy1')
                if strcmp(dataStruct(a).Manipulation,'None') || strcmp(dataStruct(a).ManipDetails,'Pre') || strcmp(dataStruct(a).Manipulation,'Pre')...
                    neuronIndRG=[neuronIndRG;a];
                end
            end
        case 'Astrocyte'            
            if strcmp(dataStruct(a).Promoter,'Aldh1l1-CreER') && (strcmp(dataStruct(a).Sensor,'GCaMP6s') || strcmp(dataStruct(a).Sensor,'GCaMP6short'))
                if strcmp(dataStruct(a).Manipulation,'None') || strcmp(dataStruct(a).ManipDetails,'Pre') || strcmp(dataStruct(a).Manipulation,'Pre')...
                    if strcmp(dataStruct(a).Tamoxifen,'Twice')    
                    astroIndALGC6=[astroIndALGC6;a];
                    elseif strcmp(dataStruct(a).Tamoxifen,'4HT twice') 
                        astroIndALGC64HT=[astroIndALGC64HT;a];
                    end
                end            
            end            
    end
end
astroIndAll=[astroIndALGC6; astroIndALGC64HT];
%% plot
%frequency
lng=([length(neuronIndSN),length(astroIndALGC6)]);
waveFrq=nan(max(lng),4);
waveFrq(:,1)=waveFrqL(neuronIndSN)';
waveFrq(:,2)=waveFrqR(neuronIndSN)';
waveFrq(1:lng(2),3)=waveFrqL(astroIndALGC6)';
waveFrq(1:lng(2),4)=waveFrqR(astroIndALGC6)';

figure;  h=plotSpread(waveFrq,'distributionColors',{'b','g','r','m'},'xValues',[0.5 0.7 1.5 1.7],'xNames',{'L','R','L','R'});
hold on;
er = errorbar([0.5 0.7 1.5 1.7],nanmean(waveFrq),nanstd(waveFrq)./[sqrt(lng(1)),sqrt(lng(1)),sqrt(lng(2)),sqrt(lng(2))],'ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 10;
ylabel('Wave frequency (events/min)')

%compare neurons to astros without L-R
lng=([length(neuronIndSN),length(astroIndALGC6)]);
waveFrq=nan(max(lng),2);
waveFrq(:,1)=nanmedian([waveFrqL(neuronIndSN);waveFrqR(neuronIndSN)])';
waveFrq(1:lng(2),2)=nanmedian([waveFrqL(astroIndALGC6);waveFrqR(astroIndALGC6)])';

figure;  h=plotSpread(waveFrq,'distributionColors',{'b','m'},'xValues',[0.5 1],'xNames',{'N','A'});
hold on;
er = errorbar([0.5 1],nanmean(waveFrq),nanstd(waveFrq)./[sqrt(lng(1)),sqrt(lng(2))],'ok');
% er.MarkerFaceColor = [0 0 0];
% er.MarkerSize = 10;
% ylim([0 0.1])
ylabel('Wave frequency (events/min)')

h=lillietest(waveFrq(:,3))
h=lillietest(waveFrq(:,4))
if vartest2(waveFrq(:,1),waveFrq(:,2))==0
    [h,p,ci,stats] = ttest2(waveFrq(:,1),waveFrq(:,2),'Vartype','equal')
else
    [h,p,ci,stats] = ttest2(waveFrq(:,1),waveFrq(:,2),'Vartype','unequal')
end
[p,h,stats] = signrank(waveFrq(:,3),waveFrq(:,4))

%compare astros with tamoxifen to astros with 4HT without L-R
lng=([length(astroIndALGC6),length(astroIndALGC64HT)]);
waveFrq=nan(max(lng),2);
waveFrq(:,1)=nanmedian([waveFrqL(astroIndALGC6);waveFrqR(astroIndALGC6)])';
waveFrq(1:lng(2),2)=nanmedian([waveFrqL(astroIndALGC64HT);waveFrqR(astroIndALGC64HT)])';

figure;  h=plotSpread(waveFrq,'distributionColors',{'b','m'},'xValues',[0.5 1],'xNames',{'Tmx','4HT'});
hold on;
er = errorbar([0.5 1],nanmean(waveFrq),nanstd(waveFrq)./[sqrt(lng(1)),sqrt(lng(2))],'ok');
% er.MarkerFaceColor = [0 0 0];
% er.MarkerSize = 10;
ylim([0 12])
ylabel('Wave frequency (events/min)')

%wave interval
avgIntrvlL=cellfun(@median,intrvlL)./10;%for sec
avgIntrvlR=cellfun(@median,intrvlR)./10;
lng=([length(neuronIndSN),length(astroIndALGC6)]);
waveInt=nan(max(lng),4);
waveInt(:,1)=avgIntrvlL(neuronIndSN)';
waveInt(:,2)=avgIntrvlR(neuronIndSN)';
waveInt(1:lng(2),3)=avgIntrvlL(astroIndALGC6)';
waveInt(1:lng(2),4)=avgIntrvlR(astroIndALGC6)';

figure;  h=plotSpread(waveInt,'distributionColors',{'b','g','r','m'},'xValues',[0.5 0.7 1.5 1.7],'xNames',{'L','R','L','R'});
hold on;
er = errorbar([0.5 0.7 1.5 1.7],nanmean(waveInt),nanstd(waveInt)./[sqrt(lng(1)),sqrt(lng(1)),sqrt(lng(2)),sqrt(lng(2))],'ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 10;
ylabel('Wave interval (s)')

%compare neurons to astros without L-R
lng=([length(neuronIndSN),length(astroIndALGC6)]);
waveInt=nan(max(lng),2);
waveInt(:,1)=nanmedian([avgIntrvlL(neuronIndSN);avgIntrvlR(neuronIndSN)])';
waveInt(1:lng(2),2)=nanmedian([avgIntrvlL(astroIndALGC6);avgIntrvlR(astroIndALGC6)])';

figure;  h=plotSpread(waveInt,'distributionColors',{'b','m'},'xValues',[0.5 1],'xNames',{'N','A'});
hold on;
er = errorbar([0.5 1],nanmean(waveInt),nanstd(waveInt)./[sqrt(lng(1)),sqrt(lng(2))],'ok');
% er.MarkerFaceColor = [0 0 0];
% er.MarkerSize = 10;
% ylim([0 0.1])
ylabel('Wave interval (s)')

h=lillietest(waveFrq(:,3))
h=lillietest(waveFrq(:,4))
if vartest2(waveFrq(:,1),waveFrq(:,2))==0
    [h,p,ci,stats] = ttest2(waveFrq(:,1),waveFrq(:,2),'Vartype','equal')
else
    [h,p,ci,stats] = ttest2(waveFrq(:,1),waveFrq(:,2),'Vartype','unequal')
end
[p,h,stats] = signrank(waveFrq(:,3),waveFrq(:,4))


%combine astrocytes with tamoxifen and 4HT
astroAll=[astroIndALGC6; astroIndALGC64HT];
mean(nanmedian([waveFrqL(astroAll);waveFrqR(astroAll)]))
std(nanmedian([waveFrqL(astroAll);waveFrqR(astroAll)]))/sqrt(length(astroAll))
nanmean(waveFrqL(astroAll))
nanmean(waveFrqR(astroAll))
nanstd(waveFrqL(astroAll))./sqrt(length(astroAll))
nanstd(waveFrqR(astroAll))./sqrt(length(astroAll))
[p,h,stats] = signrank(waveFrqL(astroAll),waveFrqR(astroAll))
mean(nanmedian([avgIntrvlL(astroAll);avgIntrvlR(astroAll)]))
std(nanmedian([avgIntrvlL(astroAll);avgIntrvlR(astroAll)]))./sqrt(length(astroAll))

%% bilaterality of peaks
astBothAll=[];
for cc=1:length(astroIndAll)
    if strcmp(dataStruct(astroIndAll(cc)).use,'Both')
        astBothAll=[astBothAll;astroIndAll(cc)];
    end
end
mean(biEventPrcnt(astBothAll)*100)
std(biEventPrcnt(astBothAll)*100)./sqrt(length(astBothAll))

%% look at different ages
ageNind=cell(1,2); %P6-P9,P10-P11
ageAind=cell(1,2);
for a=1:size(dataStruct,2)
    if strcmp(dataStruct(a).Manipulation,'None') || strcmp(dataStruct(a).ManipDetails,'Pre') ||...
            strcmp(dataStruct(a).Manipulation,'Pre') || (strcmp(dataStruct(a).Manipulation,'KOhet')...
            && strcmp(dataStruct(a).ManipDetails,'mGluR5'))
        switch dataStruct(a).Age
            case {'6','7','8','9'}
                if strcmp(dataStruct(a).Cell,'Neuron')
                    ageNind{1}=[ageNind{1};a];
                else
                    ageAind{1}=[ageAind{1};a];
                    
                end
            case {'10','11'}
                if strcmp(dataStruct(a).Cell,'Neuron')
                    ageNind{2}=[ageNind{2};a];
                else
                    ageAind{2}=[ageAind{2};a];
                end
        end
    end
end
tempYng = find(ismember(astroIndAll,cell2mat(ageAind(:,1))));
astroIndY=astroIndAll(tempYng);
tempOld = find(ismember(astroIndAll,cell2mat(ageAind(:,2))));
astroIndO=astroIndAll(tempOld);

%% plot over ages
% plot frequency
lng=([length(ageAind{1}),length(ageAind{2})]);
waveFrq=nan(max(lng),2);
waveFrq(1:lng(1),1)=nanmedian([waveFrqL(ageAind{1});waveFrqR(ageAind{1})])';
waveFrq(1:lng(2),2)=nanmedian([waveFrqL(ageAind{2});waveFrqR(ageAind{2})])';

figure;  h=plotSpread(waveFrq,'xNames',{'P6-9','P10-11'});
hold on;
er = errorbar(nanmean(waveFrq),nanstd(waveFrq)./sqrt(lng),'ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 10;
ylim([0 11])
ylabel('Wave frequency (events/min)')

h=lillietest(waveFrq(:,1))
h=lillietest(waveFrq(:,2))

[p,h,stats] = ranksum(waveFrq(:,1),waveFrq(:,2))
if vartest2(waveFrq(:,1),waveFrq(:,2))==0
    [h,p,ci,stats] = ttest2(waveFrq(:,1),waveFrq(:,2),'Vartype','equal')
else
    [h,p,ci,stats] = ttest2(waveFrq(:,1),waveFrq(:,2),'Vartype','unequal')
end

%plot Inter-event latency
lng=([length(ageAind{1}),length(ageAind{2})]);
waveInt=nan(max(lng),2);
waveInt(1:lng(1),1)=nanmedian([avgIntrvlL(ageAind{1});avgIntrvlR(ageAind{1})])';
waveInt(1:lng(2),2)=nanmedian([avgIntrvlL(ageAind{2});avgIntrvlR(ageAind{2})])';

figure;  h=plotSpread(waveInt,'distributionColors',{'b','m'},'xValues',[0.5 1],'xNames',{'6-9','10-11'});
hold on;
er = errorbar([0.5 1],nanmean(waveInt),nanstd(waveInt)./[sqrt(lng(1)),sqrt(lng(2))],'ok');
% er.MarkerFaceColor = [0 0 0];
% er.MarkerSize = 10;
% ylim([0 0.1])
ylabel('Wave interval (s)')

h=lillietest(waveInt(:,1))
h=lillietest(waveInt(:,2))
if vartest2(waveInt(:,1),waveInt(:,2))==0
    [h,p,ci,stats] = ttest2(waveInt(:,1),waveInt(:,2),'Vartype','equal')
else
    [h,p,ci,stats] = ttest2(waveInt(:,1),waveInt(:,2),'Vartype','unequal')
end
%% MPEP manipulation
mpepN=cell(1,2);
mpepA=cell(1,3);
for a=1:size(dataStruct,2)
    switch dataStruct(a).Cell
        case 'Neuron'
            if strcmp(dataStruct(a).Manipulation,'Pre') && strcmp(dataStruct(a+1).Manipulation,'MPEP')
                mpepN{1}=[mpepN{1};a];
            elseif a==18 || a==21
                mpepN{1}=[mpepN{1};a];
            elseif strcmp(dataStruct(a).Manipulation,'MPEP')
                mpepN{2}=[mpepN{2};a];
            end
        case 'Astrocyte'
            if strcmp(dataStruct(a).Manipulation,'Pre') && strcmp(dataStruct(a+1).Manipulation,'MPEP')
                mpepA{1}=[mpepA{1};a];
            elseif strcmp(dataStruct(a).Manipulation,'MPEP')
                mpepA{2}=[mpepA{2};a];
            elseif strcmp(dataStruct(a).Manipulation,'LY+MPEP')
                mpepA{3}=[mpepA{3};a];
            end
    end
end

%plot astrocyte SC and IC pre and post MPEP
% plot frequency
lng=([cellfun(@length,mpepA(1:2))]);
waveFrq=nan(max(lng),4);
waveFrq(1:lng(1),1)=nanmedian([waveFrqL(mpepA{1});waveFrqR(mpepA{1})])';
waveFrq(1:lng(2),2)=nanmedian([waveFrqL(mpepA{2});waveFrqR(mpepA{2})])';
waveFrq(1:lng(1),3)=cell2mat(ICfreq(mpepA{1}))';
waveFrq(1:lng(2),4)=cell2mat(ICfreq(mpepA{2}))';

figure; plot([1,2],[waveFrq(:,1),waveFrq(:,2)],'-ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
hold on;%plot([3,4],[waveFrq(:,3),waveFrq(:,4)],'-ob','MarkerFaceColor','b','MarkerSize',10)
er = errorbar([1:2],[nanmean(waveFrq(:,1)),nanmean(waveFrq(:,2))],[nanstd(waveFrq(:,1)),...
    nanstd(waveFrq(:,2))]./[sqrt(lng(1:2))],'-ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 12;
xlim([0.5 2.5])
ylim([0 10])
set(gca,'XTickLabel',{'','Pre MPEP','','Post MPEP',''})
ylabel('Frequency (events/min)')
title('SC')
box off
figure;% plot([1,2],[waveFrq(:,1),waveFrq(:,2)],'-ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
plot([1,2],[waveFrq(:,3),waveFrq(:,4)],'-ok','MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10)
hold on;
er = errorbar([1:2],[nanmean(waveFrq(:,3)),nanmean(waveFrq(:,4))],[nanstd(waveFrq(:,3)),...
    nanstd(waveFrq(:,4))]./[sqrt(lng(1:2))],'-ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 12;
xlim([0.5 2.5])
ylim([0 10])
set(gca,'XTickLabel',{'','Pre MPEP','','Post MPEP',''})
ylabel('Frequency (events/min)')
title('IC')
box off

%test for normality
h=lillietest(waveFrq(:,3)-waveFrq(:,4))
h=lillietest(waveFrq(:,1)-waveFrq(:,2))
%paired ttest
[h,p,ci,stats] = ttest(waveFrq(:,3),waveFrq(:,4))
[h,p,ci,stats] = ttest(waveFrq(:,1),waveFrq(:,2))

%% KOs
mG5A=cell(1,8);%control/het/KO/KO+MPEP/KO+LY/Ctl+MPEP/Ctl+MPEP+LY/het+MPEP
mG3A=cell(1,4);

for a=1:size(dataStruct,2)
    if strcmp(dataStruct(a).Cell,'Astrocyte')>0
        switch dataStruct(a).ManipDetails
            case 'mGluR5'
                if strcmp(dataStruct(a).Manipulation,'KOcontrol')
                    mG5A{1}=[mG5A{1};a];
                elseif strcmp(dataStruct(a).Manipulation,'KOhet')
                    mG5A{2}=[mG5A{2};a];
                else
                    mG5A{3}=[mG5A{3};a];
                end
            case 'mGluR5+MPEP'
                if strcmp(dataStruct(a).Manipulation,'KOcontrol')
                    mG5A{6}=[mG5A{6};a];
                elseif strcmp(dataStruct(a).Manipulation,'KOhet')
                    mG5A{8}=[mG5A{8};a];
                else
                    mG5A{4}=[mG5A{4};a];
                end
            case 'mGluR5+MPEP+LY'
                if strcmp(dataStruct(a).Manipulation,'KOcontrol')
                    mG5A{7}=[mG5A{7};a];
                end
            case 'mGluR5+LY'
                mG5A{5}=[mG5A{5};a];
            case 'Grm3cKO'
                if strcmp(dataStruct(a).Manipulation,'KOcontrol')
                    mG3A{1}=[mG3A{1};a];
                else
                    mG3A{2}=[mG3A{2};a];
                end
            case 'mGluR3+LY'
                mG3A{3}=[mG3A{3};a];
            case 'mGluR3+MPEP'
                mG3A{4}=[mG3A{4};a];
        end
    end
end


%% PLOTS FOR KOs
% look at astrocyte mglur5 only (with IC info)
mG5CtlA=[astroIndALGC64HT;mG5A{1};mG5A{2}];
lng=[length(mG5CtlA),length(mG5A{3}), length(mG5CtlA),length(mG5A{3})];
waveFrq=nan(max(lng),4);
waveFrq(1:lng(1),1)=nanmedian([waveFrqL(mG5CtlA);waveFrqR(mG5CtlA)])'; %mglur5 controls
waveFrq(1:lng(2),2)=nanmedian([waveFrqL(mG5A{3});waveFrqR(mG5A{3})])'; %mglur5 KOs
waveFrq(1:lng(3),3)=[ICfreq{mG5CtlA}]'; %mglur5 Ctls IC
waveFrq(1:lng(4),4)=[ICfreq{mG5A{3}}]'; %mglur5 KOs IC

figure;  h=plotSpread(waveFrq,'distributionColors',{'k','r','k','r'},...
    'xNames',{'Ctl-SC','mG5 KO-SC','Ctl-IC','mG5 KO-IC'});
hold on;
er = errorbar(nanmean(waveFrq),nanstd(waveFrq)./sqrt(lng),'ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 10;
ylabel('Wave frequency (events/min)')


h=lillietest(waveFrq(:,1))
h=lillietest(waveFrq(:,2))

if vartest2(waveFrq(:,1),waveFrq(:,2))==0
    [h,p,ci,stats] = ttest2(waveFrq(:,1),waveFrq(:,2),'Vartype','equal')
else
    [h,p,ci,stats] = ttest2(waveFrq(:,1),waveFrq(:,2),'Vartype','unequal')
end

h=lillietest(waveFrq(:,3))
h=lillietest(waveFrq(:,4))

if vartest2(waveFrq(:,3),waveFrq(:,4))==0
    [h,p,ci,stats] = ttest2(waveFrq(:,3),waveFrq(:,4),'Vartype','equal')
else
    [h,p,ci,stats] = ttest2(waveFrq(:,3),waveFrq(:,4),'Vartype','unequal')
end

%Astrocye mglur3 with IC info
mG3CtlA=[astroIndALGC64HT;mG3A{1}];
lng=[length(mG3CtlA),length(mG3A{2}),length(mG3A{4}),...
    length(mG3CtlA),length(mG3A{2}),length(mG3A{4})];
waveFrq=nan(max(lng),6);
waveFrq(1:lng(1),1)=nanmedian([waveFrqL(mG3CtlA);waveFrqR(mG3CtlA)])'; %mglur3 controls
waveFrq(1:lng(2),2)=nanmedian([waveFrqL(mG3A{2});waveFrqR(mG3A{2})])'; %mglur3 KOs
waveFrq(1:lng(3),3)=nanmedian([waveFrqL(mG3A{4});waveFrqR(mG3A{4})])'; %mglur3 KOs + MPEP
waveFrq(1:lng(4),4)=[ICfreq{mG3CtlA}]'; %mglur3 Ctls IC
waveFrq(1:lng(5),5)=[ICfreq{mG3A{2}}]'; %mglur3 KOs IC
waveFrq(1:lng(6),6)=[ICfreq{mG3A{4}}]'; %mglur3 KOs IC + MPEP

figure;  h=plotSpread(waveFrq,'distributionColors',{'k','r','m','k','r','m'},...
    'xNames',{'Ctl-SC','mG3 KO-SC','mG3 KO-SC MPEP','Ctl-IC',...
    'mG3 KO-IC','mG3 KO-IC MPEP'});
hold on;
er = errorbar(nanmean(waveFrq),nanstd(waveFrq)./sqrt(lng),'ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 10;
ylabel('Wave frequency (events/min)')

%Prepare data for a 3 way nested anova
%put data into a vector (you will then label each variable in the table)
respTemp=[waveFrq(find(~isnan(waveFrq(:,1))),1);...
    waveFrq(find(~isnan(waveFrq(:,2))),2);waveFrq(find(~isnan(waveFrq(:,3))),3);...
    waveFrq(find(~isnan(waveFrq(:,4))),4);waveFrq(find(~isnan(waveFrq(:,5))),5);...
    waveFrq(find(~isnan(waveFrq(:,6))),6)];

Ctrl5TblMM=table(respTemp,[ones(length(find(~isnan(waveFrq(:,1)))),1);...
    ones(length(find(~isnan(waveFrq(:,2)))),1)*2;...
    ones(length(find(~isnan(waveFrq(:,3)))),1)*2;...
    ones(length(find(~isnan(waveFrq(:,4)))),1)*1;...
    ones(length(find(~isnan(waveFrq(:,5)))),1)*2;...
    ones(length(find(~isnan(waveFrq(:,6)))),1)*2],...
    [ones(length(find(~isnan(waveFrq(:,1)))),1)*1;...
    ones(length(find(~isnan(waveFrq(:,2)))),1)*1;...
    ones(length(find(~isnan(waveFrq(:,3)))),1)*1;...
    ones(length(find(~isnan(waveFrq(:,4)))),1)*2;...
    ones(length(find(~isnan(waveFrq(:,5)))),1)*2;...
    ones(length(find(~isnan(waveFrq(:,6)))),1)*2],...
    [ones(length(find(~isnan(waveFrq(:,1)))),1);...
    ones(length(find(~isnan(waveFrq(:,2)))),1)*1;...
    ones(length(find(~isnan(waveFrq(:,3)))),1)*2;...
    ones(length(find(~isnan(waveFrq(:,4)))),1)*1;...
    ones(length(find(~isnan(waveFrq(:,5)))),1)*1;...
    ones(length(find(~isnan(waveFrq(:,6)))),1)*2],...
    [(1:length(find(~isnan(waveFrq(:,1)))))';...
    (length(find(~isnan(waveFrq(:,1))))+1:length(find(~isnan(waveFrq(:,1))))...
    +length(find(~isnan(waveFrq(:,2)))))';...
    (11:14)';(2:length(find(~isnan(waveFrq(:,1)))))';...
    (length(find(~isnan(waveFrq(:,1))))+1:length(find(~isnan(waveFrq(:,1))))...
    +length(find(~isnan(waveFrq(:,2)))))';(11:14)'],...
    'VariableNames',{'Responses','Genotype','Area','Treatment','Mice'});
Ctrl5TblMM.Genotype=categorical(Ctrl5TblMM.Genotype);
Ctrl5TblMM.Area=categorical(Ctrl5TblMM.Area);
Ctrl5TblMM.Mice=categorical(Ctrl5TblMM.Mice);

%%3 way nested ANOVA
[p,tbl,stats,terms]=anovan(respTemp,{Ctrl5TblMM.Genotype Ctrl5TblMM.Area Ctrl5TblMM.Treatment},'model',...
    'interaction','varnames',{'Gentoype','Area','Treatment'},'nested',[0 0 0;0 0 0;1 0 0]);

%% eye ablations
countABL=0;
for a=1:size(dataStruct,2)
    if strcmp(dataStruct(a).Manipulation,'Enucleation')
         if strcmp(dataStruct(a).ManipDetails,'Left')
             countABL=countABL+1;
            contraSC(countABL)=waveFrqR(a);
            contraIC(countABL)=ICfreqR{a};
            ipsiSC(countABL)=waveFrqL(a);
            ipsiIC(countABL)=ICfreqL{a};
         elseif strcmp(dataStruct(a).ManipDetails,'Right')
             countABL=countABL+1;
            contraSC(countABL)=waveFrqL(a);
            contraIC(countABL)=ICfreqL{a};
            ipsiSC(countABL)=waveFrqR(a);
            ipsiIC(countABL)=ICfreqR{a};
         end
    end
end

%% eye ablations  plot

%frequency
waveFrq=nan(3,4);
waveFrq(:,1)=ipsiSC';%SC Ipsi
waveFrq(:,2)=contraSC';%sc contra
waveFrq(:,3)=ipsiIC';%IC Ipsi
waveFrq(:,4)=contraIC';%IC contra

figure;  plot([1:2],waveFrq(:,1:2),'-o','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5,0.5,0.5],...
    'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerSize',10)
% h=plotSpread(waveFrq,'distributionColors',{'b','g','r','m'},'xValues',[0.5 0.7 1.5 1.7],'xNames',{'L','R','L','R'});
hold on;
plot([3:4],waveFrq(:,3:4),'-o','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5,0.5,0.5],...
    'MarkerEdgeColor',[0.5,0.5,0.5],'MarkerSize',10)
er = errorbar([1:4],nanmean(waveFrq),nanstd(waveFrq)./[sqrt(3)],'ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 12;
xlim([0.5 4.5])
ylabel('Wave frequency (events/min)')
set(gca,'XTickLabel',{'','Ipsi-SC','','Contra-SC','','Ipsi-IC','','Contra-IC'})


%paired ttest
[h,p,ci,stats] = ttest(waveFrq(:,1),waveFrq(:,2))
[h,p,ci,stats] = ttest(waveFrq(:,3),waveFrq(:,4))

%compare ipsiSC to control SC (astro)
lng=([length(astroIndY),length(ipsiSC)]);
waveFrq=nan(max(lng),2);
waveFrq(1:lng(1),1)=nanmedian([waveFrqL(astroIndY);waveFrqR(astroIndY)])'; %compare to controls
waveFrq(1:lng(2),2)=ipsiSC;
nanmean(waveFrq)
nanstd(waveFrq)./sqrt(lng)

% h=lillietest(waveFrq(:,1))
% h=lillietest(waveFrq(:,2))

[p,h,stats] = ranksum(waveFrq(:,1),waveFrq(:,2))

%compare both left and right IC to ablation IC
lng=([length(astroIndY),length(ipsiIC)]);
waveFrq=nan(max(lng),2);
waveFrq(1:lng(1),1)=cell2mat([ICfreq(astroIndY)])'; %compare to controls
waveFrq(1:lng(2),2)=nanmedian([ipsiIC; contraIC])';
nanmean(waveFrq)
nanstd(waveFrq)./sqrt(lng)

% h=lillietest(waveFrq(:,1))
% h=lillietest(waveFrq(:,2))
[p,h,stats] = ranksum(waveFrq(:,1),waveFrq(:,2))
