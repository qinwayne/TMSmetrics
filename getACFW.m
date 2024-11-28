function biomarkers = getACFW(response,d, biomarkers, nonartefact_length)
    data = response(:,nonartefact_length(1):end);

% for the passive biomarkers
    nonprob_detrend = detrend(response(:,nonartefact_length(end):end));
    nonprob_detrend_norm = (nonprob_detrend-mean(nonprob_detrend,2))./max(abs(nonprob_detrend-mean(nonprob_detrend,2)),[],2);
    arE(d) = 0;    
    for i = 1:size(data,1)

        %% autocorrelation
        autocorrE = autocorr(nonprob_detrend_norm(i,:),round(size(nonprob_detrend_norm,2)/2));%sum((Ehalf(1:end-1)-mean(Ehalf)).*(Ehalf(2:end)-mean(Ehalf)))/varE(j);
        arE(d) = arE(d) + sum(autocorrE>0.6)/size(data,1);

    end
    biomarkers.arE(d) = arE(d);
end