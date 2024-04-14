function [hShift,wShift] = searchBestShift(moving,ref,hShift,wShift,scl)
%     tic;
    hShift = hShift*scl;
    wShift = wShift*scl;
    [H,W] = size(moving);
    if(hShift>=0)
        ih1 = 1+hShift:H; ih2 = 1:H-hShift;
    else
        ih1 = 1:H+hShift; ih2 = 1-hShift:H;
    end
    
    if(wShift>=0)
        iw1 = 1+wShift:W; iw2 = 1:W-wShift;
    else
        iw1 = 1:W+wShift; iw2 = 1-wShift:W;
    end

    radiu = floor(scl/2);
    len = radiu*2+1;
    
    ih = ih1(1)+radiu:ih1(end)-radiu;
    iw = iw1(1)+radiu:iw1(end)-radiu;
    ref = ref(ih2(1)+radiu:ih2(end)-radiu,iw2(1)+radiu:iw2(end)-radiu);
    ref = ref(:);
    corr = zeros(len^2,1);
    
    cnt = 1;
    
    for y = -radiu:1:radiu
        for x = -radiu:1:radiu
            selectData = moving(ih+x,iw+y);
            corr(cnt) = (selectData(:))'*ref;
            cnt = cnt + 1;
        end
    end
    
    [~,id] = max(corr);
    [hchange,wchange] = ind2sub([len,len],id);
    wShift = wShift + wchange - radiu - 1;
    hShift = hShift + hchange - radiu - 1;
%     toc;
end