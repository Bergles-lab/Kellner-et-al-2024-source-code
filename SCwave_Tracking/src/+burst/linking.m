function [traLst,traScore,traLen] = linking(zscoreMaps,opts)
    %% tracking setting
    zscoreMaps(zscoreMaps<opts.zThr) = 0;
    reward = -log(opts.IoULimit);
    C_en = reward*0.5;
    C_ex = reward*0.5;
    maxJump = opts.maxJump;                 % depends on the data
    [H,W,T] = size(zscoreMaps);
    zscoreMaps(zscoreMaps<0) = 0;   % filter
    jumpProb = opts.jumpProb;
    %% since the moving speed is slow, use overlap as the measure.
    nodeMap = zeros(H*W,T);
    disp('Linking - Constructing the graph - Node');
    d_m = [];
    cntNode = 0;
    nodeLst = cell(100000,1);
    nodeFrame = zeros(100000,1);
    for t = 1:T
        % Node Property
        curMap = zscoreMaps(:,:,t);
        scoreList = unique(curMap);
        scoreList = scoreList(2:end);
        for i = 1:numel(scoreList)
           pix = find(curMap==scoreList(i));
           nodeMap(pix,t) = i + cntNode;
           nodeLst{i + cntNode} = pix;
        end
        
        %% Node in graph
        d_m0 = zeros(numel(scoreList),4);
        d_m0(:,1) = [1:numel(scoreList)] + cntNode;
        d_m0(:,2) = C_en;
        d_m0(:,3) = C_ex;
        d_m0(:,4) = -reward;
        d_m = [d_m;d_m0];
        nodeFrame([1:numel(scoreList)] + cntNode) = t;
        cntNode = cntNode + numel(scoreList);
    end
    
    nodeLst = nodeLst(1:cntNode);
    nodeFrame = nodeFrame(1:cntNode);
    
    %% Graph construction
    disp('Linking - Constructing the graph - linking');
    % Each 2D detection is one node. Use overlap
    s_t_Origin = [];
    curJump = 1;
    while curJump<=maxJump
        for t = 1:T
            if t+curJump<=T
                map1 = nodeMap(:,t);
                map2 = nodeMap(:,t+curJump);
                overLapMap = (map1>0 & map2>0);
                label1 = map1(overLapMap);
                label2 = map2(overLapMap);
                pair = unique(label1*100000 + label2);
                label1 = floor(pair/100000);
                label2 = pair - label1*100000;
                
                s_t0 = zeros(numel(pair),4);
                s_t0(:,1) = label1;
                s_t0(:,2) = label2;
                s_t0(:,4) = curJump-1;
                
                for k = 1:numel(pair)
                   s_t0(k,3) = probCal(nodeLst{label1(k)},nodeLst{label2(k)},H,W,opts); % probability of such overlap
                end
                s_t_Origin = [s_t_Origin;s_t0];
            end
        end
        curJump = curJump + 1;
    end

    %% minimum circulation for tracking - CINDA
    disp('Linking - Link nodes');
    
    s_t = zeros(size(s_t_Origin,1),3);
    s_t(:,1:2) = s_t_Origin(:,1:2);
    % calculate cost
    s_t(:,3) = -log(jumpProb) * s_t_Origin(:,4) - log(s_t_Origin(:,3));
    [trajectories, costs] = mcc4mot(d_m,s_t);
       
    %% show results
    labelMap = zeros(H,W,T);
    traLst = cell(numel(trajectories),1);
    traScore = zeros(numel(trajectories),1);
    traLen = zeros(numel(trajectories),1);
    for i = 1:numel(trajectories)
        idInCurTra = trajectories{i};
        pix = [];
        for k = 1:numel(idInCurTra)
            curPix = nodeLst{idInCurTra(k)} + H*W*(nodeFrame(idInCurTra(k))-1);
            pix = [pix;curPix];
        end
        traLst{i} = pix;
        traScore(i) = -costs(i);
        traLen(i) = numel(idInCurTra);
    end
    traLst = traLst(traScore>1e-7);
    traScore = traScore(traScore>1e-7);
    traLen = traLen(traScore>1e-7);
%     traScore = traScore + C_en + C_ex;
end
function p = probCal(pix1,pix2,H,W,opts)

    x_direction = [-1,0,1,-1,1,-1,0,1];
    y_direction = [-1,-1,-1,0,0,1,1,1];
    
    % grow pix1
    new_add = pix1;
    for i = 1:opts.growIoU
        curGrow = [];
        [x_position,y_position] = ind2sub([H,W],new_add);
        for dir_ind = 1:numel(x_direction)
            x_neighbor = min(max(x_position + x_direction(dir_ind),1),H); 
            y_neighbor = min(max(y_position + y_direction(dir_ind),1),W); 
            curGrow = [curGrow;sub2ind([H,W],x_neighbor,y_neighbor)];
        end
        new_add = setdiff(curGrow,pix1);
        pix1 = [pix1;new_add];
    end
    
    % grow pix2
    new_add = pix2;
    for i = 1:opts.growIoU
        curGrow = [];
        [x_position,y_position] = ind2sub([H,W],new_add);
        for dir_ind = 1:numel(x_direction)
            x_neighbor = min(max(x_position + x_direction(dir_ind),1),H); 
            y_neighbor = min(max(y_position + y_direction(dir_ind),1),W); 
            curGrow = [curGrow;sub2ind([H,W],x_neighbor,y_neighbor)];
        end
        new_add = setdiff(curGrow,pix2);
        pix2 = [pix2;new_add];
    end
    
    nInter = numel(intersect(pix1,pix2));
    nUnion = numel(union(pix1,pix2));
%     p1 = binocdf(nInter,numel(pix2),numel(pix1)/nUnion);
%     p2 = binocdf(nInter,numel(pix1),numel(pix2)/nUnion);
%     p = max(p1,p2);
    p = nInter/nUnion;
end