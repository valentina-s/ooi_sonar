% Generate parameter file with param combinations

% Path
save_path = '/scratch/ch153/wjl/';
save_file = sprintf('params_20191105.txt');

% Parameter combination
rank_all = 2:10;
sm_all = [5e6];  % smoothness
sp_all = [0.1, 0.2, 0.5, 10, 20];  % sparsity
rep_all = 1:10;   % repetitions
betaHW = 0.1;
max_iter = 2e4;    % max iteration

% Open text file
fileID = fopen(fullfile(save_path, save_file),'w');
fprintf(fileID,'%4s %6s %6s %10s %10s %10s %7s\n','rank','betaH','betaW','sparsity','smoothness','max_iter','rng_num');


for r = rank_all
    for sm = sm_all
        for sp = sp_all
            for rep = rep_all
                fprintf(fileID,'%4d %6.2f %6.2f %10.2e %10.2e %10.2e %7d\n',r, betaHW, betaHW, sp, sm, max_iter, rep);
            end
        end
    end
end

fclose(fileID);


