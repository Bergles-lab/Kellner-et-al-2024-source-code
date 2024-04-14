function datReg = registrate_CC(data,ff)
    %% register
    tic;
    [H,W,T] = size(data);
    scl = 3;
    dataScl = imresize(data(:,:,1),1/scl);
    [H0,W0] = size(dataScl);
    dataScl = zeros(H0,W0,T);
    for t = 1:T
        dataScl(:,:,t) = imresize(data(:,:,t),1/scl);
    end
    
    ref = mean(dataScl(:,:,1:10),3);
    refOrg = mean(data(:,:,1:10),3);
    tforms = cell(T,1);
    
    if exist('ff','var')
        waitbar(0.1,ff);
    end
        
    %% resized data
    WExtend = 0;
    HExtend = 0;
    
    for t = 1:T
        if exist('ff','var') & mod(t,100)==0
            waitbar(t/T*0.8 + 0.1,ff);
        end
        matrix = ui.algo.calCC(dataScl(:,:,t),ref);
        [~,id] = max(matrix(:));
        [hShift,wShift] = ind2sub(size(matrix),id); 
        [hShift,wShift] = ui.algo.searchBestShift(data(:,:,t),refOrg,hShift-H0,wShift-W0,scl);
        tform = eye(3);
        tform(3,1) = -wShift;
        tform(3,2) = -hShift;
        WExtend = max(WExtend,abs(wShift));
        HExtend = max(HExtend,abs(hShift));
        tforms{t} = affine2d(tform);
    end
    
%     %%  original scale
%     for t = 1:T
%         matrix = test(data(:,:,t),refOrg);
%         [~,id] = max(matrix(:));
%         [hShift,wShift] = ind2sub(size(matrix),id); 
%         tform = eye(3);
%         tform(3,1) = W-wShift;
%         tform(3,2) = H-hShift;
%         WExtend = max(WExtend,abs(wShift));
%         HExtend = max(HExtend,abs(hShift));
%         tforms{t} = affine2d(tform);
%     end
    
    datReg = zeros(H+2*HExtend,W+2*WExtend,T,'single');
    datReg(HExtend+1:HExtend+H,WExtend+1:WExtend+W,:) = data;
    valid_region = false(H+2*HExtend,W+2*WExtend,T);
    valid_region(HExtend+1:HExtend+H,WExtend+1:WExtend+W,1) = true;
    for i = 1:T
        tform = tforms{i};
        datReg(:,:,i) = imwarp(datReg(:,:,i),tform,'OutputView',imref2d([H+2*HExtend,W+2*WExtend]));
        valid_region(:,:,i) = imwarp(valid_region(:,:,1),tform,'OutputView',imref2d([H+2*HExtend,W+2*WExtend]));
    end
    Mask = sum(valid_region,3);
    Mask = Mask==T;
    ihw = find(Mask);
    [ih,iw] = ind2sub(size(Mask),ihw);
    minH = min(ih);maxH = max(ih);
    minW = min(iw);maxW = max(iw);
    datReg = datReg(minH:maxH,minW:maxW,:);
    
    if exist('ff','var')
        waitbar(1,ff);
    end
    toc;
end