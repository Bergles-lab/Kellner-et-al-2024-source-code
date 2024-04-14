function score = fitModel(zscore,fitModel2)
    if(zscore<=7)
        score = log(2*(1-normcdf(zscore)));
    else
        score = fitModel2(zscore);
    end

end