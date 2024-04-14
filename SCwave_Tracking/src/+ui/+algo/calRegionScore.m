function score = calRegionScore(com_element,datOrg,datSmo)
    x_direction = [-1,-1,-1,0,0,1,1,1];
    y_direction = [-1,0,1,-1,1,-1,0,1];

    [H,W] = size(datSmo);
    com_neighbor = [];
    com_grow = com_element;
    new_add = com_element;
    growCnt = 0;
    while(numel(com_neighbor)<numel(com_element) && growCnt<100)
        pix1 = [];
        growCnt = growCnt+ 1;
        [x_position,y_position] = ind2sub([H,W],new_add);
        for dir_ind = 1:numel(x_direction)
            x_neighbor = min(max(x_position + x_direction(dir_ind),1),H); 
            y_neighbor = min(max(y_position + y_direction(dir_ind),1),W); 
            pix1 = [pix1;sub2ind([H,W],x_neighbor,y_neighbor)];
        end
        new_add = setdiff(pix1,com_grow);
        com_grow = [com_grow;new_add];
        com_neighbor = [com_neighbor;new_add];
    end
%     com_neighbor = setdiff(com_neighbor,com_element);
    
    fg = datOrg(com_element);
    bg = datOrg(com_neighbor);
    L = mean(fg)-mean(bg);
    [mu, sigma] = burst.ksegments_orderstatistics_fin(fg, bg);
    score = (L-mu)/sigma;
end