%%
fqbd = {'Delta','Theta','Alpha','Beta','Gamma'};
for fqbdi = 1:length(fqbd)
    clearvars -except fqbdi fqbd
    clc
    close all

    folders = dir();
    nonartefact_length = 1:550;
    fs = 1000;
    Delta = [1 4];
    Theta = [4 8];
    Alpha = [8 12];
    Beta = [12 30];
    Gamma = [30 80];

    freqband = fqbd{fqbdi};
    for fid = 3:length(folders)
        cd(folders(fid).name)
        matfiles = dir('*.mat');
        for fileid = 1:length(matfiles)
            clear bioM
            load(matfiles(fileid).name)
            data = eval(matfiles(fileid).name(1:end-4));
            bioM.pili = 0;
            for chid = 1:size(data,1)
                norm_response = data(chid,:);
                bioM = getPILI(fs,norm_response,chid, bioM, nonartefact_length);
                bioM = getACFW(norm_response,chid, bioM, nonartefact_length);
                bioM = getLE(norm_response,chid, bioM, nonartefact_length);
            end
        end
    end
end

