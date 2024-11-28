
function biomarkers = getPILI(fs,response,d, biomarkers, nonartefact_length)   
    data = response(:,nonartefact_length(1):end);
    
    %% data for use
    baseline = data(:,end-end/10:end);
    data_no_dc = data(1:length(nonartefact_length)) - mean(baseline,2);

    data_norm = data_no_dc./max(abs(data_no_dc),[],2);
    data_norm_mean = mean(data_norm,1);

    prob_no_dc = response(:,nonartefact_length)-mean(baseline,2);
    prob_norm = prob_no_dc./max(abs(prob_no_dc),[],2);
    prob_norm_mean = mean(prob_norm,1);
    prob_mean = mean(prob_no_dc,1);

   
    est_aucs = 0;
    est_aucs_norm = 0;
    
    %% fit nonlinear model to find decay term
    tspan = 0:1/fs:size(data_norm_mean,2)/fs-1/fs;
    data_for_fitting = prob_norm_mean;

    for i = 1:size(data_norm_mean,1)
        %% PILI
        est_aucs_norm   =   est_aucs_norm   +   trapz(tspan,(data_norm_mean(i,:)))/size(data_norm_mean,1); %est_aucs(d)+trapz(tspan,abs(prob_data_norm(i,:)))/size(data,1);
        est_aucs        =   est_aucs        +   trapz(tspan,abs(data_norm_mean(i,:)))/size(data_norm_mean,1); %est_aucs(d)+trapz(tspan,abs(prob_data_norm(i,:)))/size(data,1);

    end
    biomarkers.pili(d) = est_aucs_norm;
    biomarkers.pili2(d) = est_aucs;
end