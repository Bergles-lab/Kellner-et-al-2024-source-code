%% This code was used to generate the raster plots in Figure 1D,E,G
%Copyright - Vered Kellner April 2024
%%
close all
clearvars

%% load file
[fn,dname1] = uigetfile('D:\Data\RAW WF data\*.czi','Open raw movie');
[pathstr, name, ext] = fileparts([dname1 fn]);
switch ext
    case '.czi'
        bf = bfopen([dname1 fn]);
        tic;
        [m,n] = size(bf{1}{1});
        t = size(bf{1},1);
        img = zeros(m,n,t,'int16');
        for i=1:t
            img(:,:,i) = bf{1}{i};
        end
        clear bf;
        %         img = imrotate(img,180); %rotate if needed
        toc;
  
        %% bleach correction
        sampRate = 10; %sampling rate in Hz
        [imgBC] = bleachCorrect(img, sampRate);
        imgBC=imrotate(imgBC,180);
    case '.tif'
        imgBC = loadTif([dname1 fn],16);
        sampRate = 10; %sampling rate in Hz
        figure;
        imagesc(mean(imgBC,3)); axis image; colormap gray
        clear img
end
%% crop image for grid analysis
%run script separately for left and right SC
figure;
imagesc(mean(imgBC,3)); axis image; colormap gray; caxis([500 5000])
brect = imrect(gca,[0,0,96,84]);%100

setResizable(brect,1);
wait(brect);
pos = getPosition(brect);
pos = int16(round(pos));
crimg=imgBC((pos(2)):(pos(2)+pos(4)),pos(1):(pos(1)+pos(3)),:);
close
imgCrop=crimg; clear crimg
LRstr={'LSC','RSC'};
LRans=listdlg('ListString',LRstr);

%% get dFF
%currently normalizing to median, can be changed in normalizeImg to
%percentile
[dFoF,Fo]=normalizeImg(double(imgCrop));
sigTrace=squeeze(nanmean(nanmean(dFoF),2));
sampRate=10; %manually change as needed
win=5; %check that the window is correct for your data
t=size(dFoF,3);
T=[0:1/sampRate:(t/sampRate)-1/sampRate];
sigTraceF=msbackadj(T',smooth(sigTrace),'SHOWPLOT',0,'WindowSize',win,'StepSize',win);

%% Grid ROIS

figure; imshow(mean(imgCrop,3)/mean(max(mean(imgCrop,3))));
widthImg = size(imgCrop,2);
heightImg = size(imgCrop,1);
sizeSq = 10; %adjust manually as needed

[positiveIndices] = getGrid(widthImg,heightImg,sizeSq);
T=size(imgCrop,3);
for i=1:size(positiveIndices,1)
    hold on;
    tempX=positiveIndices(i,1:2:end);
    tempY=positiveIndices(i,2:2:end);
        pgon = polyshape(tempX,tempY);
    plot(tempX,tempY,'Color','g');
    if positiveIndices(i,6)<heightImg
        xs = [positiveIndices(i,2):positiveIndices(i,6)];
    else
        xs=[positiveIndices(i,2):heightImg];
    end
    if positiveIndices(i,3)<widthImg
        ys = [positiveIndices(i,1):positiveIndices(i,3)];
    else
        ys=[positiveIndices(i,1):widthImg];
    end
%     plot(xs,ys,'Color','g');
    
    temp = squeeze(mean(mean(dFoF(xs,ys,:),2),1));
    rois(:,i) = smooth(msbackadj([1:T]',temp,'WindowSize',20,'StepSize',20,'Showplot',0));
%if want to see the numbers uncomment both lines
%         [cX,cY]=centroid(pgon);
%         text(cX,cY,num2str(i))
end
%% crop larger image for movement analysis
figure;
imagesc(mean(imgBC,3)); axis image; colormap gray; caxis([500 5000])
brect = imrect(gca,[0,0,349,149]);%250

setResizable(brect,1);
wait(brect);
pos = getPosition(brect);
pos = int16(round(pos));
crimg=imgBC((pos(2)):(pos(2)+pos(4)),pos(1):(pos(1)+pos(3)),:);
close
imgCropLG=crimg; clear crimg

%% find movement times 
imgCropLG=double(imgCropLG);
norms = [];
regto = fft2(mean((imgCropLG),3));
t = size(imgCropLG,3);
parfor j = 1:t
    [output, Greg] = dftregistration(regto,fft2(imgCropLG(:,:,j)),2) %the last number will determine how sensitive this is
    norms(j) = norm([1 0 output(3); 0 1 output(4); 0 0 1],2);
end

hit = zeros(size(norms));
for j = 1:t
    if norms(j) > 1
        %look 30 ahead
        if j>5 && j < t-30
            if sum(norms(j-1:j+30)-1) > 0
                hit(j-5:j+30) = 1;            
            end
        elseif j<5
            if sum(norms(j-1:j+30)-1) > 0
                hit(1:j+30) = 1;            
            end
        elseif j>t-30
            hit(j-5:end) = 1;
        end
    end
end
figure; plot(norms-1); hold on; plot(hit)
mvmInd=find(hit>0);

%% plot grids
[t,n] = size(rois);
numToShow = n; %choose a number here that is not too large ~50
figure('Position',[200,150,1200,800]); plot(rois(:,1:numToShow) - .1*repmat(1:numToShow,t,1),'Color','k'); 

%% remove movement
roisNMVM=rois;
roisNMVM(mvmInd,:)=0; %for some cases need to zero 
roisNMVM2=rois;
roisNMVM2(mvmInd,:)=nan; %generally better to use nan but not for everything
sigTraceF2=sigTraceF;
sigTraceF2(mvmInd)=nan;

%% find peaks
[t,n] = size(roisNMVM);
dataMat=[]; %pks, frms,hw,roi
figure('Position',[200,150,1200,800]);
for cc=1:size(roisNMVM,2) %go through each ROI
    threshAst(cc)=nanmedian(roisNMVM2(:,cc))+3*mad(roisNMVM2(:,cc));
    [pksAst,locsAst,wAst,pAst] = findpeaks((double(roisNMVM(:,cc))),1:t,'MinPeakHeight',threshAst(cc),...
        'MinPeakWidth',15,'MinPeakProminence', threshAst(cc)/2,'WidthReference','halfheight');
    plot(roisNMVM(:,cc) - .1*repmat(cc,t,1),'Color','k');
    hold on; plot(locsAst-.1*cc,pksAst-.1*cc,'*r')
    dataMat=[dataMat;[pksAst,locsAst',wAst',pAst,repmat(cc,length(locsAst),1)]];    
end
dataMatSrt=sortrows(dataMat,2); %sort according to peak location (frames)
tempDiff=diff(dataMatSrt(:,2));
tempDiff=[0;tempDiff]; %to make this the correct n

%% Remove movement related events
[reps, vals]=groupcounts(dataMatSrt(:,2));
repsBadInd=find(reps>=5);%was 5 for astros
valsBadInd=vals(repsBadInd);
dataMatSrtNew=dataMatSrt;
for n=1:length(valsBadInd)    
   tempInd=find(dataMatSrtNew(:,2)==valsBadInd(n));
   dataMatSrtNew(tempInd,:)=[];
   tempInd=[];
end

figure('Position',[200,150,1200,800]);
for cc=1:size(dataMatSrtNew,1) %go through each ROI
    tempInd=dataMatSrtNew(cc,5);
    plot(rois(:,tempInd) - .1*repmat(tempInd,t,1),'Color','k');
    hold on; plot(dataMatSrtNew(cc,2)-.1*tempInd,dataMatSrtNew(cc,1)-.1*tempInd,'*r')    
end
%% make raster plot of detected peaks
%first arrange the data
roiMat=sortrows(dataMatSrtNew,5);
roiCell=cell(size(rois,2),1); %need to treat each roi as a trial
for nn=1:size(rois,2)
    tempInd=find(roiMat(:,5)==nn);
    roiCell(nn,1)={roiMat(tempInd,2)'./10}; %event times in seconds            
end

LineFormat = struct();
LineFormat.Color = 'k';
LineFormat.LineWidth = 1.5;
figure;  plotSpikeRaster(roiCell,'PlotType','vertline','LineFormat',LineFormat);
