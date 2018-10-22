% 2017 11 28  Rename files to form better sequence

clear

addpath ~/internal_2tb/ooi_sonar/ooi_sonar_code/

data_path = '/media/wu-jung/wjlee_apl_2/ooi_zplsc_palm_nmf/';

cd(data_path)  % has to switch folder to ensure reading h5 files correctly

n_comp = 3;
files = dir(fullfile(data_path,sprintf('palm_nmf_r%d_beta*.mat',n_comp)));


for iF=1:length(files)
    s = strsplit(files(iF).name,'.mat');

    % distinguish if beta or eta param first in filename
    if strcmp(s{1}(13),'e')
        eta = str2double(s{1}(16:end));
        newname = sprintf('palm_nmf_r%d_eta%06.2f.mat',n_comp,eta);
    else
        betaW = str2double(s{1}(18:20));
        betaH = str2double(s{1}(27:29));
        eta = str2double(s{1}(34:end));
        newname = sprintf('palm_nmf_r4_betaW0.1_betaH0.1_eta%06.2f.mat',eta);
    end

    disp(['oldname:', files(iF).name]);
    disp(['newname:', newname]);

    if ~strcmp(files(iF).name,newname)
        movefile(files(iF).name,newname);
    end
end
