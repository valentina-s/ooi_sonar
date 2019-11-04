% Test writing and reaiding formatted txt files

save_path = '/scratch/ch153/wjl/';

rank_all = 2:10;   % rank
betaHW = reshape(repmat([0, 0.1, 1, 100, 1000, 10000],2,1),1,[]);
sm = 5e6;        % smoothness
sp = 2;          % sparsity (lambda)
max_iter = 5e3;  % max iteration

for irank=1:length(rank_all)

    r = rank_all(irank);
    save_file = sprintf('params_20191027_betaHW_rank%02d.txt',r);

    % write text file
    fileID = fopen(fullfile(save_path, save_file),'w');
    fprintf(fileID,'%4s %5s %5s %9s %9s %9s\n',...
            'rank','betaH','betaW','sparsity','smoothness','max_iter');
    
    for ibeta = 1:length(betaHW)
        betaH = betaHW(ibeta);
        betaW = betaHW(ibeta);
    
	fprintf(fileID,'%4d %05.2f %05.2f %09.2e %09.2e %09.2e\n',...
                r, betaH, betaW, sp, sm, max_iter);
    end

    fclose(fileID);
    
end
