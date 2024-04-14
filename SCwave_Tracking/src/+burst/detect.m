function zscoreMaps = detect(datOrg,opts)
    %% stablization
    [H,W,T] = size(datOrg);
    gap = 1;
    xx = (datOrg(:,:,gap+1:end) - datOrg(:,:,1:end-gap)).^2;
    stdN = sqrt(nanmedian(xx,3)/0.9133);  % use median to estimate the noise, avoid the signal impact
    stdN = median(stdN(:));
    datOrg = datOrg/stdN;

    %% Use the smoothed data
    datSmo = imgaussfilt(datOrg,opts.smoXY);
    x_direction = zeros(1,9);
    y_direction = zeros(1,9);
    ind_dir = 1;
    for xx = -1:1
        for yy = -1:1
            x_direction(ind_dir) = xx;
            y_direction(ind_dir) = yy;
            ind_dir = ind_dir + 1;
        end
    end
    zscoreMaps = zeros(H,W,T);
    thresholds = [opts.thrARScl:opts.stepRatio:max(datSmo(:))];
    
    parfor t = 1:T
        zscore_map = zeros(H,W);
        %% from low to high
        for kk = 1:numel(thresholds)
            thre = thresholds(kk);
            detection = datSmo(:,:,t)>=thre;
%             detection = imerode(detection,strel('disk',3));
%             detection = imdilate(detection,strel('disk',3));
            components = bwconncomp(detection);
            pixLst = components.PixelIdxList;
            area = regionprops(components,'area'); area = [area.Area];
            perimeter = regionprops(components,'perimeter'); perimeter = [perimeter.Perimeter];
            circularity  = (4*pi*area)./perimeter.^2;

            sz = cellfun(@numel,pixLst);
            pixLst = pixLst(sz>=opts.minSize & sz<=opts.maxSize & circularity > opts.minCir);
            components_len = numel(pixLst);
            zscores = zeros(components_len,1);
            for com_ind = 1:components_len
                com_element = pixLst{com_ind};
                zscores(com_ind) = ui.algo.calRegionScore(com_element,datOrg(:,:,t),datSmo(:,:,t));
                zscore_map(com_element) = max(zscores(com_ind),zscore_map(com_element));
            end

        end
        foreground = imregionalmax(zscore_map);
        zscore_map(~foreground) = 0;
        zscoreMaps(:,:,t) = zscore_map;
    end

end