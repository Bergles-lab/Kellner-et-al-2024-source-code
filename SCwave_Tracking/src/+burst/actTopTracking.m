function [dFOrg] = actTopTracking(datOrg,opts,evtSpatialMask,ff)
    
    opts.smoT = 0;
    opts.pixNoise = 0;
    opts.smoXY = 0.5;

    % valid region
    if ~isfield(opts,'fgFluo')
        opts.fgFluo = 0;
    end
    [H,W,T] = size(datOrg);
    datOrgMean = mean(datOrg,3);
    msk000 = var(datOrg,0,3)>1e-8;
    
    if exist('evtSpatialMask','var') && ~isempty(evtSpatialMask)
        evtSpatialMask = evtSpatialMask.*msk000;
    else
        evtSpatialMask = msk000;
    end
    noiseEstMask = evtSpatialMask.*(datOrgMean>opts.fgFluo);
        
    nonvalid = datOrg==0;
    se = strel(true(2*opts.smoXY+1,2*opts.smoXY+1,5));
    nonvalid = imdilate(nonvalid,se);
    
    % smooth the data
    dat = datOrg;
    if opts.smoXY>0
        for tt=1:size(dat,3)
            dat(:,:,tt) = imgaussfilt(dat(:,:,tt),opts.smoXY);
        end
    end
    
    % noise for smoothed data
    gap = 1;
    xx = dat;
    xx(nonvalid) = nan;
    xx = (xx(:,:,gap+1:end) - xx(:,:,1:end-gap)).^2;
    stdMap = sqrt(nanmedian(xx,3)/0.9133);
    
    if opts.smoT>0
%          dat = imgaussfilt3(dat,[1e-3,1e-3,opts.smoT]);
%          datOrg = imgaussfilt3(datOrg,[1e-3,1e-3,opts.smoT]);
         dat = movmean(dat,opts.smoT*2+1,3);
         datOrg = movmean(datOrg,opts.smoT*2+1,3);
    end
    

    %% Bayesian Deal With noise
    if(opts.pixNoise)
        stdMapGau = burst.noiseEstimate(dat,nonvalid,gap,evtSpatialMask,opts);
    else
        xx = (dat(:,:,gap+1:end) - dat(:,:,1:end-gap)).^2;
        stdN = sqrt(nanmedian(xx,3)/0.9133);  % use median to estimate the noise, avoid the signal impact
        stdN = median(stdN(:));
        stdMapGau = ones(H,W)*stdN;
    end
    
    %% linear estimation of F0 new
    tic;
    disp('Baseline removing');
    opts.cut = min(opts.cut,T);
    datIn = dat./repmat(stdMapGau,1,1,T);   % stablization, make large SNR pixel has more weight for F0. Also for parameter setting.
    [dF,xBias] = burst.getDfBlk_linear_removeNaN2(datIn,opts.cut,opts.movAvgWin,nonvalid);
    dFOrg = datOrg - dat + dF.*repmat(stdMapGau,1,1,T); % dFOrg - F0
    toc;
    opts.stdMap = stdMapGau;
    opts.xBias = xBias;
    
    if exist('ff','var')
        waitbar(1,ff);
    end
    
end

