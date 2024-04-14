%% Tracking of SC waves
%Code was developed by Xuelong Mi in the lab of Guoqiang Yu and adapted by
%Vered Kellner 
close all
clearvars -except dataStruct
f=1;
%% setting
startup;    % initialize
preset = 1; % default
opts = util.parseParam(preset,1);

cellAns=questdlg('Cell type?','','Neuron','Astrocyte','Neuron');
useList={'Left','Right','Both'};
useAns=listdlg('ListString',useList,'PromptString','Use?');
%% large-impact paramters
if strcmp(cellAns,'Astrocyte')
    opts.lenFilter = 20;     % the minimum duration of signal wave. Will filter out short signal waves - was 3 for neurons
    opts.maxJump = 20;       % allowed max jump (due to merging, some detection will not be found) - was 8 for neurons
%     opts.lenFilter = 10;     % the minimum duration of signal wave. Will filter out short signal waves - was 3 for neurons
%     opts.maxJump = 10;       % allowed max jump (due to merging, some detection will not be found) - was 8 for neurons
% %     % For example, if one signal wave appear at frame 1 and frame 5, then we at
    % least set this parameter to 4 to link them.
    
    opts.jumpProb = 0.9;    % Set from 0 to 1. It's the false negative probability - was 0.5 for neurons
    % (the probablity I didn't detect the signal in one frame.) But due to the
    % merging, I set it to 0.5. Large this parameter, more easily to link
    % detections for large jump.
    %% other parameters
    opts.thrARScl = 0.2;    % bottom-up to check detection, start from (opts.thrARScl x noise) was 3 for neurons, 0.2 for astro
    opts.smoXY = 5;         % smooth parameter for select region was 1 for neurons
    opts.minCir = 0.05;     % minimum circularity for one detection was 0.05
    opts.stepRatio = 0.1;   % the step for bottom-up. Should be small if SNR is low. Now is 0.1 x noise one step
    opts.cut = 30;          % parameter for calculating dF. If the F0 change strongly, should be small.
    opts.movAvgWin = 10;    % parameter for calculating dF. If the F0 change strongly, should be small.
    opts.growIoU = 2;       % parameter for calculating IoU. grow the region then judge IoU.
    
else
    opts.lenFilter = 20; %was 20 for neurons
    opts.maxJump = 20;
    opts.jumpProb = 0.9;
    %% other parameters
    opts.thrARScl = 0.1;    % bottom-up to check detection, start from (opts.thrARScl x noise) was 1 for neurons then 2 then 0.5
    opts.smoXY = 3;         % smooth parameter for select region was 2 for neurons
    opts.minCir = 0.05;     % minimum circularity for one detection - was 0.05
    opts.stepRatio = 0.1;   % the step for bottom-up. Should be small if SNR is low. Now is 0.1 x noise one step
    opts.cut = 30;          % parameter for calculating dF. If the F0 change strongly, should be small.
    opts.movAvgWin = 10;    % parameter for calculating dF. If the F0 change strongly, should be small.
    opts.growIoU = 2;       % parameter for calculating IoU. grow the region then judge IoU.
end
opts.minSize = 100;     % minimum number of pixels for one detection in each frame
opts.maxSize = 5000;    % maximize number of pixels for one detection in each frame

% linking parameters
opts.IoULimit = 0.01;    % IoU limitation for linking together
opts.zThr = 0.2;          % filter small-score - was 3 for neurons was 0.2 for dim  was 0.5


%% file path
p0 = 'D:\Data\RAW WF data\';  % folder name
f0 = uigetfile([p0,'*.czi']);  % file name

%% load data
load('random_Seed');
rng(s);
[folder, name, ext] = fileparts(strcat(p0,'\',f0));
path0 = [p0,name,'\'];
% if ~exist(path0,'dir') && ~isempty(path0)
%     mkdir(path0);    
% end
[datOrg,opts] = burst.prep1(p0,f0,[],opts);  % read data

%% 1. registration
% datReg = ui.algo.registrate_CC(datOrg);

%% 2. background
[dFOrg] = burst.actTopTracking(datOrg,opts);  % read data

%% find movement times (from Travis / Vered)
Y=datOrg;
norms = [];
regto = fft2(mean((Y),3));
t = size(Y,3);
parfor j = 1:t
[output, Greg] = dftregistration(regto,fft2(Y(:,:,j)),1) %the last number will determine how sensitive this is
norms(j) = norm([1 0 output(3); 0 1 output(4); 0 0 1],2);
end
hit = zeros(size(norms));
for j = 1:t
if norms(j) > 1
    %look 50 ahead
    if j>5 && j < t-30
        if sum(norms(j+1:j+30)-1) > 0
            hit(j-5:j+30) = 1;
        elseif hit(j-5) == 1
            %                 hit(i:i+5) = 1;
            %                 hit(i+6:i+56)=0;
            hit(j:j+1) = 1;
            hit(j+2:j+30)=0;
        else
            hit(j) = 0;
        end
    end
end
end
figure; plot(norms-1); hold on; plot(hit)
mvmInd=find(hit>0);

%% 3. detection
tic
zScoreMap = burst.detect(dFOrg,opts);
toc

%% 4. linking
tic
[evtLst,~,evtLen] = burst.linking(zScoreMap,opts);
toc

%% 5. filter
evtLen2=evtLen;
evtLst2=evtLst;
evtLst = evtLst(evtLen>opts.lenFilter);
evtLen = evtLen(evtLen>opts.lenFilter);

%% feature
[H,W,T] = size(datOrg);
movFea = cell(numel(evtLst),1);
% determine if events are in left or right SC
figure; imagesc(squeeze(mean(datOrg,3)))
lh = imline(gca);
[s] = round(lh.getPosition); x = s(1); 
% % Original scale dF, without stabilization (we use square root to stabilize
% % data), and smooth filter to denoise.
dF = imgaussfilt(dFOrg.*(2*datOrg-dFOrg)*opts.maxValueDat,3);
for i = 1:numel(evtLst)

vel = nan(1,T);
% first row of dir: positive - north, negative - south
% second row of dir: positive - east, negative - west
dir = nan(2,T);     

[ih,iw,it] = ind2sub([H,W,T],evtLst{i});
t0 = min(it);
t1 = max(it);
ts = unique(it); %active frames for this event
centroid = nan(2,T);
for j = 1:numel(ts)
  curIt = ts(j); 
  curIh = ih(it==curIt);
  curIw = iw(it==curIt);
  curPix = sub2ind([H,W],curIh,curIw);
  weight = dFOrg(:,:,curIt);
  weight = weight(curPix);
  weight = weight/sum(weight);
  centroidX = weight'*curIh;
  centroidY = weight'*curIw;
  centroid(1,curIt) = H+1-centroidX;
  centroid(2,curIt) = centroidY;
end

vel(t0) = 0;
dir(:,t0) = 0;
pret = t0;
t = t0+1;
while(t<=t1)
  if(~isnan(centroid(1,t)))
      prePos = centroid(:,pret);
      curPos = centroid(:,t);
      move = curPos-prePos;
      dur = t-pret;
      vel(pret+1:t) = sqrt(sum(move.^2))/dur;
      dir(:,pret+1:t) = repmat(move/sqrt(sum(move.^2)),1,dur);  
      pret = t;
  end
  t = t+1;
end
if useAns==1
   evtLoc='Left';
elseif useAns==2
   evtLoc='Right';
elseif useAns==3
   if nanmedian(centroid(2,:))>x
       evtLoc='Right';
   else
       evtLoc='Left';
   end
end
movFea{i}.centroid = centroid;
movFea{i}.vel = vel;
movFea{i}.dir = dir;
movFea{i}.frms=ts; %active frames
movFea{i}.Side=evtLoc; %left or right
movFea{i}.amplitude = median(dF(evtLst{i}));
end

%% show
tic
ov = zeros(3*H,W,3,T);
ov(1:H,:,1,:) = datOrg;
ov(1:H,:,2,:) = datOrg;
ov(1:H,:,3,:) = datOrg;
ov(H+1:2*H,:,:,:) = double(plt.regionMapWithData_label(evtLst,movFea,datOrg*0.5,0.5))/255; 
ov(2*H+1:end,:,1,:) = dFOrg/max(dFOrg(:));
ov(2*H+1:end,:,2,:) = ov(2*H+1:end,:,1,:);
ov(2*H+1:end,:,3,:) = ov(2*H+1:end,:,1,:);
zzshow(ov);
toc
%% manually filter out events - VK added Jan 2022
%need to compare the above video with raw video, write down in notepad
%which roi numbers to keep and input here
userAns=inputdlg('Which events to keep?');
keepInd=str2num(userAns{:});
keepInd=unique(keepInd);

%% create structure for each mouse
dataStruct(f).Cell=cellAns;
sensorList={'GCaMP3','GCaMP6s','RGECO','GCaMP6short'};
sensorAns=listdlg('ListString',sensorList,'PromptString','Sensor:');
dataStruct(f).Sensor=sensorList{sensorAns};
promoterList={'GLAST-CreER','SNAP25','Aldh1l1-CreER','Thy1'};
promoterAns=listdlg('ListString',promoterList,'PromptString','Promoter:');
dataStruct(f).Promoter=promoterList{promoterAns};
dataStruct(f).Name=name;
ageList={'5','6','7','8','9','10','11','12','13','14','15'};
ageAns=listdlg('ListString',ageList,'PromptString','Age (days):');
dataStruct(f).Age=ageList{ageAns};
sexList={'Male','Female','?'};
sexAns=listdlg('ListString',sexList,'PromptString','Sex:');
dataStruct(f).Sex=sexList{sexAns};
tmxList={'None','Once','Twice','4HT twice'};
tmxAns=listdlg('ListString',tmxList,'PromptString','Tamoxifen injection:');
dataStruct(f).Tamoxifen=tmxList{tmxAns};
manipList={'None','LY','Vehicle','MPEP','KOfull','KOhet','KOcontrol','Pre','LY+MPEP'};
manipAns=listdlg('ListString',manipList,'PromptString','Manipulation:');
dataStruct(f).Manipulation=manipList{manipAns};
switch manipAns
    case 1
        dataStruct(f).ManipDetails='None';
    case 2       
        lyList={'Pre','immediate','10min','30min','1hr','3hr','5.5hr'};
        lyAns=listdlg('ListString',lyList);
        dataStruct(f).ManipDetails=lyList{lyAns};
    case 3
        vehList={'Pre','24hr CTEP','Saline 30min','Saline immediate','Saline 3hr','ddH2O immediate','ddH2O 30min'};
        vehAns=listdlg('ListString',vehList);
        dataStruct(f).ManipDetails=vehList{vehAns};
    case 4
        mpepList={'Pre','immediate','10min','30 min','1-2hr','4-6hr','8-9hr','19-21hr','@5min','@3min','3x & @5min'};
        mpepAns=listdlg('ListString',mpepList);
        dataStruct(f).ManipDetails=mpepList{mpepAns};        
    case 5
        KOList={'FMR1','FMR1cKO','IP3R2','mGluR5','Grm3KO','Grm3cKO','mGluR5+MPEP','mGluR5+LY','mGluR5+MPEP+LY','mGluR3+MPEP','mGluR3+LY','mGluR3+LY+MPEP'};
        KOAns=listdlg('ListString',KOList);
        dataStruct(f).ManipDetails=KOList{KOAns};
    case 6
        KOList={'FMR1','FMR1cKO','IP3R2','mGluR5','Grm3KO','Grm3cKO','mGluR5+MPEP','mGluR5+LY','mGluR5+MPEP+LY','mGluR3+MPEP','mGluR3+LY','mGluR3+LY+MPEP'};
        KOAns=listdlg('ListString',KOList);
        dataStruct(f).ManipDetails=KOList{KOAns};
    case 7
        KOList={'FMR1','FMR1cKO','IP3R2','mGluR5','Grm3KO','Grm3cKO','mGluR5+MPEP','mGluR5+LY','mGluR5+MPEP+LY','mGluR3+MPEP','mGluR3+LY','mGluR3+LY+MPEP'};
        KOAns=listdlg('ListString',KOList);
        dataStruct(f).ManipDetails=KOList{KOAns};
    case 8
        dataStruct(f).ManipDetails='None';
    case 9
        lyMPList={'Pre','immediate','10min','30min','1hr','3hr','5.5hr'};
        lyAns=listdlg('ListString',lyMPList);
        dataStruct(f).ManipDetails=lyMPList{lyAns};
end
cmntList={'Good','OK','Slightly cloudy','Pretty cloudy','Bad','1024x1024','No/low activity','Animal moving a lot','pinworms'};
cmntAns=listdlg('ListString',cmntList,'PromptString','Brain state');
dataStruct(f).comment=cmntList{cmntAns};
% useList={'Left','Right','Both'};
% useAns=listdlg('ListString',useList,'PromptString','Use?');
dataStruct(f).use=useList{useAns};

dataStruct(f).Events.Duration=evtLen(keepInd); %frames
dataStruct(f).Events.Area=evtLst(keepInd); %pixels in 3D use [x,y,t]=ind2sub([H,W,T],evtlst{i}) to extract coordinates of pixel
dataStruct(f).Events.WaveInfo=movFea(keepInd); %centroid, direction, velocity, active frames, side, amp
dataStruct(f).Events.Frm=T;
dataStruct(f).Events.MvmInd=mvmInd;
