% 2017 11 28  Rename files to form better sequence

clear

addpath ~/internal_2tb/ooi_sonar/ooi_sonar_code/

data_path = '/media/wu-jung/wjlee_apl_2/ooi_zplsc_new/';

cd(data_path)  % has to switch folder to ensure reading h5 files correctly

files = dir(fullfile(data_path,'*.mat'));


for iF=1:length(files)
    s = strsplit(files(iF).name,'.mat');
    num = str2num(s{1}(13:end));
    newname = sprintf('palm_nmf_eta%010.1f.mat',num);
    if ~strcmp(files(iF).name,newname)
        movefile(files(iF).name,newname);
    end
end
