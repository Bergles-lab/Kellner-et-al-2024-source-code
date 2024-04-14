function [newCost] = costUpdate(id1,id2,jump,centerNode,NodeMap,trajectories,nodeFrameID,traVar,moveCost,opts)
    %% first trajectory and second trajectory
    if(NodeMap(id1,1)==0)
        tra1 = id1;
        traVar1 = (opts.minSigma)^2;
    else
        tra1 = trajectories{NodeMap(id1,1)};
        tra1 = tra1(max(NodeMap(id1,2)-opts.vecNodes,1):NodeMap(id1,2));
        traVar1 = traVar{NodeMap(id1,1)};
        traVar1 = traVar1(NodeMap(id1,2),1);
    end
    if(NodeMap(id2,1)==0)
        tra2 = id2;
        traVar2 = (opts.minSigma)^2;
    else
        tra2 = trajectories{NodeMap(id2,1)};
        tra2 = tra2(NodeMap(id2,2):min(NodeMap(id2,2)+opts.vecNodes,numel(tra2)));
        traVar2 = traVar{NodeMap(id2,1)};
        traVar2 = traVar2(NodeMap(id2,2),2);
    end

    %% predict position
    x1 = centerNode(:,id1);
    x2 = centerNode(:,id2);

    if(numel(tra1)==1)
        cost1 = sqrt(sum((x2-x1).^2))*moveCost;
    else
        motion1 = centerNode(:,tra1(end)) - centerNode(:,tra1(1));
        motionDur1 = nodeFrameID(tra1(end)) - nodeFrameID(tra1(1));
        predict1 = x1 + motion1/motionDur1*(1+jump);
        dist1 = sqrt(sum((x2 - predict1).^2));
        sigma1 = max(sqrt(traVar1 + 2*(traVar1)/(motionDur1)^2*(1+jump)^2),1);
        cost1 = dist1/sigma1;
    end
    if(numel(tra2)==1)
        cost2 = sqrt(sum((x2-x1).^2))*moveCost;
    else
        motion2 = centerNode(:,tra2(end)) - centerNode(:,tra2(1));
        motionDur2 = nodeFrameID(tra2(end)) - nodeFrameID(tra2(1));
        predict2 = x2 - motion2/motionDur2*(1+jump);
        dist2 = sqrt(sum((x1 - predict2).^2));
        sigma2 = max(sqrt(traVar2 + 2*(traVar2)/(motionDur2)^2*(1+jump)^2),1);
        cost2 = dist2/sigma2;
    end
    newCost = sqrt((cost1^2 + cost2^2)/2);
end