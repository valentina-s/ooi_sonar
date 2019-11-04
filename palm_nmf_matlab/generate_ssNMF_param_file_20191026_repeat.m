% Test writing and reaiding formatted txt files

save_path = '/scratch/ch153/wjl/';

repeat_num = 100;
rank_all = 2:10;   % rank

for irank=1:length(rank_all)

    r = rank_all(irank);
    save_file = sprintf('params_20191026_rank%02d.txt',r);
    
    betaH = 0.1;
    betaW = 0.1;
    sm = 5e6;        % smoothness
    sp = 1;          % sparsity (lambda)
    max_iter = 2e4;  % max iteration
    
    % write text file
    fileID = fopen(fullfile(save_path, save_file),'w');
    fprintf(fileID,'%4s %5s %5s %9s %9s %9s\n','rank','betaH','betaW','sparsity','smoothness','max_iter');

    for repeat=1:repeat_num
	fprintf(fileID,'%4d %05.2f %05.2f %09.2e %09.2e %09.2e\n',r, betaH, betaW, sp, sm, max_iter);
    end
    fclose(fileID);
    
end
