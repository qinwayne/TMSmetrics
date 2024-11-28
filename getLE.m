
function biomarkers = getLE(response,d, biomarkers, nonartefact_length)
    response = response(:,nonartefact_length);%nonartefact_length(end)
    nonprob_data = zscore(mean(response,1));% normalise

    % wolf

    best_tau = 1;
    best_m = 3;
    
    datcnt = length(nonprob_data);
    tau = best_tau;
    ndim = best_m;
    ires = 10;
    maxbox = 6000;
    db = basgen_new(nonprob_data, tau, ndim, ires, datcnt, maxbox);

    %
    dt = 1;
    evolve = 2; %20
    dismin = 0.001;
    dismax = 0.3;
    thmax = 30;
    
    [out, SUM] = fet(db, dt, evolve, dismin, dismax, thmax);
    biomarkers.LE_est(:,d) = out(:,4);

end
