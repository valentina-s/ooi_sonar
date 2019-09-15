% Writing and reading formatted txt files
% New version on 2019/08/02 to just do sparsity=0

save_path = '/scratch/ch153/wjl/';

save_file = sprintf('params_20190802.txt');
fileID = fopen(fullfile(save_path, save_file),'w');
fprintf(fileID,'%4s %5s %5s %9s %9s %9s\n','rank','betaH','betaW','sparsity','smoothness','max_iter');

rank_all = 3:10;   % rank

for irank=1:length(rank_all)
    
    r = rank_all(irank);
    
    betaH = 0.1;
    betaW = 0.1;
    sp = 0;  % sparsity (lambda)
    max_iter = 2e4;    % max iteration
    
    % smoothness
    % sm_order = [5,6,7,8];  % runs sent on 20190623
    sm_order = [1,2,3,4];   % runs sent on 20190626
    sm = repmat([1,2,5],length(sm_order),1);
    for iorder = 1:length(sm_order)
        sm(iorder,:) = sm(iorder,:)*10^sm_order(iorder);
    end
    sm = sm';
    sm = sm(:);
    
    % write text file
    for isp=1:length(sp)
        for ism=1:length(sm)
            fprintf(fileID,'%4d %05.2f %05.2f %09.2e %09.2e %09.2e\n',r, betaH, betaW, sp(isp), sm(ism), max_iter);
        end
    end
    
end

fclose(fileID);
