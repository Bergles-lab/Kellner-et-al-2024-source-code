%% group analysis of SC dual imaging 
%load('C:\Users\vered\Documents\MATLAB\dualStructSC_new.mat')
close all
clearvars -except dualstruct

%% get info
 frmRate=2;
 aFrqSC=nan(size(dualstruct,2),1);
 aIntSC=cell(size(dualstruct,2),1);
 aDurSC=nan(size(dualstruct,2),1);
 nFrqSC=nan(size(dualstruct,2),1);
 nIntSC=cell(size(dualstruct,2),1);
 nDurSC=nan(size(dualstruct,2),1);
 aNcorrSC=nan(size(dualstruct,2),1);
 aNlagSC=nan(size(dualstruct,2),1);
 countSC=0;
 aFrqIC=nan(size(dualstruct,2),1);
 aIntIC=cell(size(dualstruct,2),1);
 aDurIC=nan(size(dualstruct,2),1);
 nFrqIC=nan(size(dualstruct,2),1);
 nIntIC=cell(size(dualstruct,2),1);
 nDurIC=nan(size(dualstruct,2),1);
 aNcorrIC=nan(size(dualstruct,2),1);
 aNlagIC=nan(size(dualstruct,2),1);
 allFrmAAmpSC=[];
 allFrmNAmpSC=[];
 allFrmSimSC=[];
 allFrmAAmpIC=[];
 allFrmNAmpIC=[];
 allFrmSimIC=[];
 countIC=0;
 for aa=1:size(dualstruct,2)
     if strcmp(dualstruct(aa).Area,'SC')
         countSC=countSC+1;
         frms=dualstruct(aa).Frm;
         aFrqSC(aa,1)=size(dualstruct(aa).Pks.AstStats,2)/(frms/frmRate)*60; %events per minute
         aIntSC{countSC,1}=diff([dualstruct(aa).Pks.AstStats.pkFrm]); %interval between events
         aDurSC(aa,1)=median([dualstruct(aa).Pks.AstStats.pkHW])/frmRate; %half height in secs
         nFrqSC(aa,1)=size(dualstruct(aa).Pks.NeurStats,2)/(frms/frmRate)*60; %events per minute
         nIntSC{countSC,1}=diff([dualstruct(aa).Pks.NeurStats.pkFrm]); %interval between events
         nDurSC(aa,1)=median([dualstruct(aa).Pks.NeurStats.pkHW])/frmRate; %half height in secs
         aNcorrSC(aa,1)=dualstruct(aa).Corr.CorrVal;
         aNlagSC(aa,1)=dualstruct(aa).Corr.Timelag; %in seconds 
         allFrmAAmpSC=[allFrmAAmpSC;dualstruct(aa).Corr.AllFrmAstAmp];
         allFrmNAmpSC=[allFrmNAmpSC;dualstruct(aa).Corr.AllFrmNAmp];
         allFrmSimSC=[allFrmSimSC;mean(dualstruct(aa).Corr.aAmpSim,2)];
     else
         countIC=countIC+1;
         frms=dualstruct(aa).Frm;
         aFrqIC(aa,1)=size(dualstruct(aa).Pks.AstStats,2)/(frms/frmRate)*60; %events per minute
         aIntIC{countIC,1}=diff([dualstruct(aa).Pks.AstStats.pkFrm]); %interval between events
         aDurIC(aa,1)=median([dualstruct(aa).Pks.AstStats.pkHW])/frmRate; %half height in secs
         nFrqIC(aa,1)=size(dualstruct(aa).Pks.NeurStats,2)/(frms/frmRate)*60; %events per minute
         nIntIC{countIC,1}=diff([dualstruct(aa).Pks.NeurStats.pkFrm]); %interval between events
         nDurIC(aa,1)=median([dualstruct(aa).Pks.NeurStats.pkHW])/frmRate; %half height in secs
         aNcorrIC(aa,1)=dualstruct(aa).Corr.CorrVal;
         aNlagIC(aa,1)=dualstruct(aa).Corr.Timelag; %in seconds
         allFrmAAmpIC=[allFrmAAmpIC;dualstruct(aa).Corr.AllFrmAstAmp];
         allFrmNAmpIC=[allFrmNAmpIC;dualstruct(aa).Corr.AllFrmNAmp];
         allFrmSimIC=[allFrmSimIC;mean(dualstruct(aa).Corr.aAmpSim,2)];
     end
 end
scICpairInds=[7,9;10,11;12,13;15,16;17,18];
 %% calculate actual time lag 
 for aa=1:size(dualstruct,2)
     astFrmInds=[dualstruct(aa).Pks.AstStats.pkFrm];
     neurFrmInds=[dualstruct(aa).Pks.NeurStats.pkFrm];
     pkLocMat=zeros(1,2);
     indMat=zeros(1,2);
     for a=1:length(astFrmInds)
         if strcmp(dualstruct(aa).Area,'IC')
             tempInd=astFrmInds(a)-6:astFrmInds(a)+1;
         else
             tempInd=astFrmInds(a)-10:astFrmInds(a)+1;
         end
         tempN=find(ismember(neurFrmInds,tempInd),1,'Last');
         if ~isempty(tempN)
             pkLocMat(a,:)=[astFrmInds(a),neurFrmInds(tempN)];
             indMat(a,:)=[a,tempN];
         else
             pkLocMat(a,:)=[astFrmInds(a),nan];
             indMat(a,:)=[a,nan];
         end
     end
     dualstruct(aa).Pks.IndMat=indMat;
     astNeurdiff=pkLocMat(:,1)-pkLocMat(:,2);
     if strcmp(dualstruct(aa).Area,'SC')
         avgDelaySC(aa)=nanmedian(astNeurdiff)/2; %in seconds
         avgDelayIC(aa)=nan;
     else
         avgDelayIC(aa)=nanmedian(astNeurdiff)/2; %in seconds
         avgDelaySC(aa)=nan;
     end
 end
 
 %average lag time (already in seconds)
 nanmean(avgDelaySC(scICpairInds(:,1)))
 nanstd(avgDelaySC(scICpairInds(:,1)))./sqrt(sum(~isnan(avgDelaySC(scICpairInds(:,1)))))
 
 nanmean(avgDelayIC(scICpairInds(:,2)))
 nanstd(avgDelayIC(scICpairInds(:,2)))./sqrt(sum(~isnan(avgDelayIC(scICpairInds(:,2)))))
 
 h=lillietest(avgDelaySC(scICpairInds(:,1)))
h=lillietest(avgDelayIC(scICpairInds(:,2)))
[p,h,stats] = ranksum(avgDelaySC(scICpairInds(:,1)),avgDelayIC(scICpairInds(:,2)))

%% plot frequency
figure; plot(1:2,[aFrqSC(scICpairInds(:,1)),nFrqSC(scICpairInds(:,1))],'-o','Color','b','MarkerFaceColor','b','MarkerSize',10);
hold on;
er = errorbar([1 2],nanmean([aFrqSC(scICpairInds(:,1)),nFrqSC(scICpairInds(:,1))]),nanstd([aFrqSC(scICpairInds(:,1)),nFrqSC(scICpairInds(:,1))])./...
    sqrt([size(scICpairInds,1),size(scICpairInds,1)]),'-ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 12;
xlim([0.5 2.5])
ylim([0 10])
ylabel('Wave frequency (events/min)')
set(gca,'XTickLabel',{'','Astrocytes','','Neurons',''})
box off

h=lillietest(aFrqSC(scICpairInds(:,1)))
h=lillietest(nFrqSC(scICpairInds(:,1)))

%paired ttest
[h,p,ci,stats] = ttest(aFrqSC(scICpairInds(:,1)),nFrqSC(scICpairInds(:,1)))

figure; plot(1:2,[aFrqIC(scICpairInds(:,2)),nFrqIC(scICpairInds(:,2))],'-o','Color','b','MarkerFaceColor','b','MarkerSize',10);
hold on;
er = errorbar([1 2],nanmean([aFrqIC(scICpairInds(:,2)),nFrqIC(scICpairInds(:,2))]),nanstd([aFrqIC(scICpairInds(:,2)),nFrqIC(scICpairInds(:,2))])./...
    sqrt([size(scICpairInds,1),size(scICpairInds,1)]),'-ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 12;
xlim([0.5 2.5])
ylim([0 10])
ylabel('Wave frequency (events/min)')
set(gca,'XTickLabel',{'','Astrocytes','','Neurons',''})
box off

h=lillietest(aFrqIC(scICpairInds(:,2)))
h=lillietest(nFrqIC(scICpairInds(:,2)))
%paired ttest
[h,p,ci,stats] = ttest(aFrqIC(scICpairInds(:,2)),nFrqIC(scICpairInds(:,2)))

%compare IC to SC
h=lillietest(aFrqSC(scICpairInds(:,1)))
h=lillietest(aFrqIC(scICpairInds(:,2)))
%paired ttest
[h,p,ci,stats] = ttest(aFrqSC(scICpairInds(:,1)),aFrqIC(scICpairInds(:,2)))

h=lillietest(nFrqSC(scICpairInds(:,1)))
h=lillietest(nFrqIC(scICpairInds(:,2)))
%paired ttest
[h,p,ci,stats] = ttest(nFrqSC(scICpairInds(:,1)),nFrqIC(scICpairInds(:,2)))

%% look at different ages
ageSCind=cell(1,2); %P6-P7-P8-P9,P10-P11
ageICind=cell(1,2);
for a=1:size(dualstruct,2)
        switch dualstruct(a).Age            
            case {'6','7','8','9'} 
                if strcmp(dualstruct(a).Area,'SC')
                    ageSCind{1}=[ageSCind{1};a];
                else
                    ageICind{1}=[ageICind{1};a];
                end
            case {'10','11'}
                if strcmp(dualstruct(a).Area,'SC')
                    ageSCind{2}=[ageSCind{2};a];
                else
                    ageICind{2}=[ageICind{2};a];
                end            
        end    
end

%% plot frequency according to age
lng=([cellfun(@length,ageSCind)]);
waveFrq=nan(max(lng),4);
waveFrq(1:lng(1),1)=aFrqSC(ageSCind{1}); %P6-P7-8-9
waveFrq(1:lng(1),2)=nFrqSC(ageSCind{1});
waveFrq(1:lng(2),3)=aFrqSC(ageSCind{2});%%P10-11
waveFrq(1:lng(2),4)=nFrqSC(ageSCind{2});

figure;  h=plot([0.5 1],waveFrq(:,1:2)','-o','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10);
hold on; plot([1.5 2],waveFrq(:,3:4)','-o','Color',[0.5 0.5 0.5],'MarkerFaceColor',[0.5 0.5 0.5],'MarkerSize',10);
er = errorbar([0.5 1 1.5 2],nanmean(waveFrq),nanstd(waveFrq)./sqrt([lng(1),lng(1),lng(2),lng(2)]),'ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 12;
ylabel('Wave frequency (events/min)')
ylim([0 3])
xlim([0 2.5])
box off
title('SC')

%prepare data for anova
y=[aFrqSC(ageSCind{1})',nFrqSC(ageSCind{1})',aFrqSC(ageSCind{2})',nFrqSC(ageSCind{2})'];
g1={'A','A','A','A','A','A','A','A','A',...
'N','N','N','N','N','N','N','N','N','A','A','A','A','N','N','N','N'};
g2=[repmat(1,1,18),repmat(2,1,8)];
[p,tbl,stats,terms]=anovan(y,{g1,g2},'model','interaction');
figure; results = multcompare(stats,'Dimension',[1 2])

%% plot correlations
figure; plot(1:2,[aNcorrIC(scICpairInds(:,2)),aNcorrSC(scICpairInds(:,1))],'o','Color','b','MarkerFaceColor','b','MarkerSize',10);
hold on;
er = errorbar([1 2],nanmean([aNcorrIC(scICpairInds(:,2)),aNcorrSC(scICpairInds(:,1))]),nanstd([aNcorrIC(scICpairInds(:,2)),aNcorrSC(scICpairInds(:,1))])./...
    sqrt([size(scICpairInds,1),size(scICpairInds,1)]),'ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 12;
xlim([0.5 2.5])
ylim([0 1])
ylabel('Correlation between astrocytes and neurons')
set(gca,'XTickLabel',{'','IC','','SC',''})
box off

h=lillietest(aNcorrIC(scICpairInds(:,2)))
h=lillietest(aNcorrSC(scICpairInds(:,1)))
if vartest2(aNcorrIC(scICpairInds(:,2)),aNcorrSC(scICpairInds(:,1)))==0
    [h,p,ci,stats] = ttest2(aNcorrIC(scICpairInds(:,2)),aNcorrSC(scICpairInds(:,1)),'Vartype','equal')
else
    [h,p,ci,stats] = ttest2(aNcorrIC(scICpairInds(:,2)),aNcorrSC(scICpairInds(:,1)),'Vartype','unequal')
end

avgLagSC=nanmean(aNlagSC)
semLagSC=nanstd(aNlagSC)./sqrt(sum(~isnan(aNlagSC)))
avgLagIC=nanmean(aNlagIC)
semLagIC=nanstd(aNlagIC)./sqrt(sum(~isnan(aNlagIC)))

h=lillietest(aNlagSC)
h=lillietest(aNlagIC)
[p,h,stats] = ranksum(aNlagSC,aNlagIC)

figure; plot(1:2,[aNcorrSC,aNlagSC],'o','Color','b','MarkerFaceColor','b');
hold on;
er = errorbar([1 2],nanmean([aNcorrSC,aNlagSC]),nanstd([aNcorrSC,aNlagSC])./...
    sqrt([sum(~isnan(aNcorrSC)),sum(~isnan(aNlagSC))]),'ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 10;
xlim([0.5 2.5])
ylim([0 3])
ylabel('Correlation between astrocytes and neurons')
set(gca,'XTickLabel',{'','Correlation','','Lag',''})

%plot correlations according to age
lng=([cellfun(@length,ageSCind)]);
corrMat=nan(max(lng),2);
corrMat(1:lng(1),1)=aNcorrSC(ageSCind{1});%P6-P7-8-9
corrMat(1:lng(2),2)=aNcorrSC(ageSCind{2});%P10-11
figure;  h=plotSpread(corrMat,'xNames',{'P6-9','P10-11'});
hold on;
er = errorbar(nanmean(corrMat),nanstd(corrMat)./sqrt([lng(1),lng(2)]),'ok');
er.MarkerFaceColor = [0 0 0];
er.MarkerSize = 12;
ylabel('Correlation')
ylim([0 1])
box off
h=lillietest(corrMat(:,1))
h=lillietest(corrMat(:,2))
if vartest2(corrMat(:,1),corrMat(:,2))==0
    [h,p,ci,stats] = ttest2(corrMat(:,1),corrMat(:,2),'Vartype','equal')
else
    [h,p,ci,stats] = ttest2(corrMat(:,1),corrMat(:,2),'Vartype','unequal')
end
%% plot intervals
figure;
for a=1:length(aIntSC)-sum(cellfun(@isempty,aIntSC))
    hold on;
    [N edges]=histcounts(aIntSC{a}/frmRate,[0:20:200]);
    plot(edges(1:end-1),N,'b')
end

figure;
for a=1:length(aIntIC)-sum(cellfun(@isempty,aIntIC))
    hold on;
    [N edges]=histcounts(aIntIC{a}/frmRate,[0:20:200]);
    plot(edges(1:end-1),N,'r')
end

%% plot average event with time - 
% calculate actual time lag
countSC=0;
for aa=1:size(dualstruct,2)
    if strcmp(dualstruct(aa).Area,'SC')
        countSC=countSC+1;
        astFrmInds=[dualstruct(aa).Pks.AstStats.pkFrm];
        neurFrmInds=[dualstruct(aa).Pks.NeurStats.pkFrm];
        pkLocMat=zeros(1,2);
        indMat=zeros(1,2);
        for bb=1:length(astFrmInds)
            if astFrmInds(bb)-10<0
                tempInd=1:astFrmInds(bb);
            else
                tempInd=astFrmInds(bb)-10:astFrmInds(bb);
            end
            tempN=find(ismember(neurFrmInds,tempInd),1,'Last');
            if ~isempty(tempN)
                pkLocMat(bb,:)=[astFrmInds(bb),neurFrmInds(tempN)];
                indMat(bb,:)=[bb,tempN];
            else
                pkLocMat(bb,:)=[nan,nan];
                indMat(bb,:)=[nan,nan];
            end
        end
        
        astNeurdiff{countSC}=pkLocMat(:,1)-pkLocMat(:,2);
        avgDelay(countSC)=nanmean(astNeurdiff{countSC})/2 %in seconds
        
        indMatnew=indMat(~isnan(indMat(:,1)),:);
        eventNyes{countSC}=[dualstruct(aa).Traces.eventNTrace(indMatnew(:,2),:)];
        eventTmN{countSC}=[dualstruct(aa).Traces.eventNTraceTm(indMatnew(:,2),:)];
        eventA{countSC}=[dualstruct(aa).Traces.eventATrace(indMatnew(:,1),:)];
        eventTmA{countSC}=[dualstruct(aa).Traces.eventATraceTm(indMatnew(:,1),:)];
        
        
        %find the difference between the times of events
        normTmN=[];
        normTmA=[];
        for a=1:size(eventTmA{countSC},1)
            lng=[length(eventTmN{countSC}(a,:)),length(eventTmA{countSC}(a,:))];
            newTm=nan(2,max(lng));
            newTm(1,1:length(eventTmN{countSC}(a,:)))=eventTmN{countSC}(a,:);%neur
            newTm(2,1:length(eventTmA{countSC}(a,:)))=eventTmA{countSC}(a,:);%ast
            tmDffAN=newTm(1,:)-newTm(2,:);
            normTmN(a,:)=newTm(1,:)-newTm(1,1)+1;
            normTmA(a,:)=normTmN(a,:)+abs(tmDffAN(:,1));
            %         tempFirstdff=sum(isnan(normTmN(1,:)))-sum(isnan(eventTmA{a}(1,:)));
            %         normTmA(1,:)=normTmA(1,:)+tempFirstdff;
            
        end
        avgTmN{countSC}=nanmean(normTmN,1);
        avgTmA{countSC}=nanmean(normTmA,1);
        eventN{countSC}=[dualstruct(aa).Traces.eventNTrace(indMatnew(:,2),:)];
        eventA{countSC}=[dualstruct(aa).Traces.eventATrace(indMatnew(:,1),:)];
    end
end

%% plot neuron astrocyte average trace - stopped here

%bring down to 0
tempN=cellfun(@nanmean,eventN,'UniformOutput',false);
mxN=max(cellfun(@length,tempN));
avgEventNMat=nan(length(tempN),mxN);
for b=1:length(tempN)
    avgEventNMat(b,1:length(tempN{b}))=tempN{b};
end
normEventN=nanmean(avgEventNMat,1);
normEventN=normEventN-normEventN(1);
tempA=cellfun(@nanmean,eventA,'UniformOutput',false);
mxA=max(cellfun(@length,tempA));
avgEventAMat=nan(length(tempA),mxA);
for b=1:length(tempA)
    avgEventAMat(b,1:length(tempA{b}))=tempA{b};
end
normEventA=nanmean(avgEventAMat,1);
normEventA=normEventA-normEventA(1);

%find sems for events
semEventN=nanstd(avgEventNMat)./sqrt(size(avgEventNMat,1));
semEventA=nanstd(avgEventAMat)./sqrt(size(avgEventAMat,1));

%get average time
mxNT=max(cellfun(@length,avgTmN));
avgEventTmNMat=nan(length(avgTmN),mxNT);
for b=1:length(avgTmN)
    avgEventTmNMat(b,1:length(avgTmN{b}))=avgTmN{b};
end

mxAT=max(cellfun(@length,avgTmA));
avgEventTmAMat=nan(length(avgTmA),mxAT);
for b=1:length(avgTmA)
    avgEventTmAMat(b,1:length(avgTmA{b}))=avgTmA{b};
end
if size(avgEventTmAMat,2)~=length(normEventA)
    sz=min([size(avgEventTmAMat,2),length(normEventA)]);
    avgEventTmAMatnew=avgEventTmAMat(:,1:sz);
end
figure; h1=shadedErrorBar(round(nanmean(avgEventTmNMat))./2,smooth(normEventN)./max(smooth(normEventN)),(semEventN),'lineprops','-m','transparent',1);
hold on; h2=shadedErrorBar(round(nanmean(avgEventTmAMatnew))./2,smooth(normEventA)./max(smooth(normEventA)),(semEventA),'lineprops','-g','transparent',1);
xlim([0 35])
% ylim([-0.02 0.20])
legend([h1.mainLine,h2.mainLine],'Neurons','Astrocytes')
xlabel('Relative time (s)')
ylabel('Normalized amplitude')

%% Figure 2E - astro amp after neuronal peak
figure; 
scatter(allFrmNAmpSC,allFrmAAmpSC,'filled','MarkerFaceColor','k','MarkerEdgeColor','k')
xlabel('Neuronal event amplitude (A.U.)');
ylabel('Astrocyte amplitude 1.5s after neuronal event (A.U.)')
%legend('Randomized data','Real data','Location','Northwest')
lsline
[R, P]=corrcoef(allFrmNAmpSC,allFrmAAmpSC)
