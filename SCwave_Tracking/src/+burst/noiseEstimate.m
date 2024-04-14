function stdMapGau = noiseEstimate(xx,nonvalid,gap,evtSpatialMask,opts)
    [H,W,T] = size(xx);
    
    %% median to estimate the noise
    xx(nonvalid) = nan;
%     stdMap = zeros(H,W);
%     for x = 1:H
%         for y = 1:W
%             dif = (xx(x,y,gap+1:end) - xx(x,y,1:end-gap)).^2;
%             stdMap(x,y) = sqrt(nanmedian(dif)/0.9133);
%         end
%     end
    xx = (xx(:,:,gap+1:end) - xx(:,:,1:end-gap)).^2;
    stdMap = sqrt(nanmedian(xx,3)/0.9133);  % use median to estimate the noise, avoid the signal impact
    stdMap(~evtSpatialMask) = nan;  % only in selected part
    stdMap(stdMap==0) = nan;    % the part with no variation shouldn't be used
    mu0 = nanmean(stdMap(:));
    nForEachPixel = T-sum(isnan(xx),3);  % # for estimated std
    clear xx;
%     if(~opts.isAlignmentNanData)   % if the part data is nonvalid for a large part, use all valid pixels to estimate only one std
%         stdMap = ones(H,W) * mu0;
%     end
    
    %% Bayesian Deal With noise
    % post probability to estimate the variance
    % E[y|x] = (sigma_0^2 x + sigma_1^2 mu)/(sigma_0^2 + sigma_1^2)
    sig0square = nanstd(stdMap(:)).^2;
    stdMap(isnan(stdMap)) = mu0;
    sig1square = stdMap.^2/2./nForEachPixel;
    stdMap_afterDeal = (stdMap*sig0square + mu0*sig1square)./(sig0square+sig1square);
    stdMapGau = double(imgaussfilt(stdMap_afterDeal)) + 1e-6;
    stdMapGau(stdMapGau<max(stdMapGau(:))/30) = max(stdMapGau(:))/30;   % avoid too small std

end