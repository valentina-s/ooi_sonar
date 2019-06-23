% Test writing and reaiding formatted txt files

save_path = '/scratch/ch153/wjl/';

rank_all = 3:10;   % rank

for irank=1:length(rank_all)
    
    r = rank_all(irank);
    save_file = sprintf('params_20190622_rank%02d.txt',r);
    
    betaH = 0.1;
    betaW = 0.1;
    sp = [1,2,5,10,20,50];  % sparsity (lambda)
    max_iter = 2e4;    % max iteration
    
    % smoothness
    sm_order = [5,6,7,8];
    sm = repmat([1,2,5],length(sm_order),1);
    for iorder = 1:length(sm_order)
        sm(iorder,:) = sm(iorder,:)*10^sm_order(iorder);
    end
    sm = sm';
    sm = sm(:);
    
    % write text file
    fileID = fopen(fullfile(save_path, save_file),'w');
    fprintf(fileID,'%4s %5s %5s %9s %9s %9s\n','rank','betaH','betaW','sparsity','smoothness','max_iter');
    for isp=1:length(sp)
        for ism=1:length(sm)
            fprintf(fileID,'%4d %05.2f %05.2f %09.2e %09.2e %09.2e\n',r, betaH, betaW, sp(isp), sm(ism), max_iter);
        end
    end
    fclose(fileID);
    
end
