%% Two Photon imaging analysis using ROIs
%as in Figure 3G-H
%Copyright - Vered Kellner 2024
% load('C:\Users\vered\Documents\MATLAB\SCIC_2P_30ROIs.mat')
clearvars -except twoPROIs; 
close all
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

cellAns=questdlg('Cell type?','','Neuron','Astrocyte','Astrocyte');
areaInp=questdlg('Area?','','SC','IC','SC');

%% Try drawing ROIS
%first crop area of interest
figure; h=imagesc(mean(img,3)/mean(max(mean(img,3))));colormap gray; caxis([0 1])
brect = imrect(gca,[0,0,159,159]);%100

setResizable(brect,1);
wait(brect);
pos = getPosition(brect);
pos = int16(round(pos));
crimg=img((pos(2)):(pos(2)+pos(4)),pos(1):(pos(1)+pos(3)),:);
imgCropL=crimg; clear crimg

figure; h=imagesc(mean(imgCropL,3)/mean(max(mean(imgCropL,3))));colormap gray; caxis([0 1])
nROI = imrect(gca,[0,0,20,20]);%100
wait(nROI);
Nmask = createMask(nROI, h);
for a=2:30    
    nROI = imrect(gca,[0,0,20,20]);%100
    wait(nROI);
    Nmask = Nmask+createMask(nROI, h);
end

Lrois=bwlabel(Nmask);
    
%% average signal for each roi
t=size(img,3);
for ii=1:max(max(Lrois))
    [tempR, tempC]=find(Lrois==ii);
    tempN=squeeze(mean(mean(img(tempR,tempC,:))));
    roisNeur(:,ii) = (msbackadj([1:t]',tempN,'WindowSize',45,'StepSize',45,'Showplot',0));
end

%% find peaks
dataNeurMat=[]; %pks, frms,hw,roi
t=size(img,3);
for cc=1:size(roisNeur,2) %go through each ROI
    if strcmp(cellAns,'Astrocyte')
        threshNeur=median(roisNeur(:,cc))+3*mad(roisNeur(:,cc));
    [pksNeur,locsNeur,wNeur] = findpeaks(roisNeur(:,cc),1:t,'MinPeakHeight',threshNeur,'MinPeakProminence',threshNeur,...
        'MinPeakWidth',3,'MaxPeakWidth',15,'WidthReference','halfheight','annotate','extents');
    figure; findpeaks(roisNeur(:,cc),1:t,'MinPeakHeight',threshNeur,'MinPeakProminence',threshNeur,...
        'MinPeakWidth',3,'MaxPeakWidth',15,'WidthReference','halfheight','annotate','extents');
    else
        if strcmp(areaInp,'IC')
        threshNeur=median(roisNeur(:,cc))+1*mad(roisNeur(:,cc));
        else
            threshNeur=median(roisNeur(:,cc))+3*mad(roisNeur(:,cc));
        end
        [pksNeur,locsNeur,wNeur] = findpeaks(roisNeur(:,cc),1:t,'MinPeakHeight',threshNeur,'MinPeakProminence',threshNeur,...
        'WidthReference','halfheight','annotate','extents');
    figure; findpeaks(roisNeur(:,cc),1:t,'MinPeakHeight',threshNeur,'MinPeakProminence',threshNeur,...
        'WidthReference','halfheight','annotate','extents');
    end
        dataNeurMat=[dataNeurMat;[pksNeur,locsNeur',wNeur',repmat(cc,length(locsNeur),1)]];
end

%% plot

[t,n] = size(roisNeur);
numToShow = n;
figure('Position',[200,150,1200,800]); plot(roisNeur(:,1:numToShow) + 5*repmat(1:numToShow,t,1),'Color','k');

%% Now gather all the data based on those events
[t,n] = size(roisNeur);
i=10; %seconds
frmRate=2;
eventT=[-i*frmRate:i*frmRate]; %5 seconds
for aa=1:n %over all ROIS
    tempInd=find(dataNeurMat(:,4)==aa);
    tempEventTrace=nan(length(tempInd),length(eventT));
    for bb=1:length(tempInd) %go through all the peaks
         if dataNeurMat(tempInd(bb),2)-i*frmRate<=0           
        elseif dataNeurMat(tempInd(bb),2)+i*frmRate>size(img,3)          
        else
            tempEventTrace(bb,1:end)=roisNeur(dataNeurMat(tempInd(bb),2)-i*frmRate:dataNeurMat(tempInd(bb),2)+i*frmRate,aa);
        end

    end
    allEvents{aa}=tempEventTrace;
    hiEventInd=find(dataNeurMat(tempInd,1)>=prctile(dataNeurMat(tempInd,1),75));
    eventHi(aa,:)=nanmean(tempEventTrace(hiEventInd,:),1);
    eventNTrace(aa,:)=nanmean(tempEventTrace,1);    
end
figure; plot(eventNTrace')
figure; plot(eventHi')
figure; plot(allEvents{3}')
%% put into struct
twoPROIs(f).Name=name;
twoPROIs(f).Cell=cellAns;
twoPROIs(f).Area=areaInp;
twoPROIs(f).frmRate=frmRate;
twoPROIs(f).frmNum=size(img,3);
twoPROIs(f).eventInfo=dataNeurMat;
twoPROIs(f).eventNTrace=eventNTrace;
twoPROIs(f).eventNTraceHi=eventHi;
twoPROIs(f).eventNTraceAll=allEvents;
twoPROIs(f).ROIs=Lrois;
