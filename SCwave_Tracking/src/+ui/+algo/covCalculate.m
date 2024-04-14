function covMatrix = covCalculate(dF,sizeH,sizeW)
    [H,W,T] = size(dF);
    lH = (sizeH-1)/2;
    lW = (sizeW-1)/2;
    num = 1e5;
    covMatrix = zeros(sizeH,sizeW);
    
    %% selected samples
    pix = randi([1,H*W*T],num,1);
    [ih,iw,it] = ind2sub([H,W,T],pix);
    selected = (ih>=lH+1) & (iw>=lW+1) & (ih <= H-lH) & (iw<=W-lW);
    ih = ih(selected);
    iw = iw(selected);
    it = it(selected);
    pixSelected = sub2ind([H,W,T],ih,iw,it);
    % dF is already normalized, noise variance is 1
    for x = 0:lH
       for y = 0:lW
           if(x==0 && y==0)
                covMatrix(lH+1,lW+1) = 1; 
           else
                pixCur = sub2ind([H,W,T],ih+x,iw+y,it);
                rho = (2 - mean((dF(pixSelected)-dF(pixCur)).^2))/2;
                covMatrix(lH+1+x,lW+1+y) = rho;
                covMatrix(lH+1-x,lW+1-y) = rho;
                covMatrix(lH+1-x,lW+1+y) = rho;
                covMatrix(lH+1+x,lW+1-y) = rho;
           end
       end
    end
    covMatrix(covMatrix<0) = 0;


end