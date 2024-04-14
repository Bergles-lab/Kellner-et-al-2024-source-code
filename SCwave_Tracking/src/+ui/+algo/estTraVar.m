function traVar = estTraVar(tra,nodeFrameID,centerNode,opts)
    
    T = numel(tra);
    
    traVar = ones(T,2)*(opts.minSigma)^2;    
    if(T<=3)
        return;
    end
    % forward 
    sumVar = 0;
    for t = 3:T
        curS = tra(t);
        s1 = tra(max(t-opts.vecNodes,1));
        s2 = tra(t-1);
        v = (centerNode(:,s2) - centerNode(:,s1))/(nodeFrameID(s2)-nodeFrameID(s1));
        predict = centerNode(:,s2) + v*(nodeFrameID(curS)-nodeFrameID(s2));
        varF(t) = sum((centerNode(:,curS) - predict).^2);
    end
    %% backward
    sumVar = 0;
    for t = T-2:-1:1
        curS = tra(t);
        s1 = tra(min(t+opts.vecNodes,T));
        s2 = tra(t+1);
        v = (centerNode(:,s2) - centerNode(:,s1))/(nodeFrameID(s2)-nodeFrameID(s1));
        predict = centerNode(:,s2) + v*(nodeFrameID(curS)-nodeFrameID(s2));
        varB(t) = sum((centerNode(:,curS) - predict).^2);
    end
    varF(1:2) = varB(1:2);
    varB(T-1:T) = varF(T-1:T);
    for t = 1:T
        t0 = max(t-opts.vecNodes,1);
        traVar(t,1) = mean(varF(t0:t));
    end
    for t = 1:T
        t1 = min(t+opts.vecNodes,T);
        traVar(t,2) = mean(varF(t:t1));
    end
    traVar = min(traVar,1000);
end