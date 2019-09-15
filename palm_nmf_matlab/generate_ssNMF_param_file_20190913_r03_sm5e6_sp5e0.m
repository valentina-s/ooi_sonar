% Writing and reading formatted txt files
% New version on 2019/09/13 to do repeated run on the same parameter combination

save_path = '/scratch/ch153/wjl/';

save_file = sprintf('params_20190913_r03_sm5e6_sp5e0.txt');

fileID = fopen(fullfile(save_path, save_file),'w');
fprintf(fileID,'%4s %5s %5s %9s %9s %9s\n','rank','betaH','betaW','sparsity','smoothness','max_iter');

repeat_num = 50;
r = 3;

betaH = 0.1;
betaW = 0.1;
sm = 5e6;   % smoothness
sp = 5;  % sparsity (lambda)
max_iter = 2e4;    % max iteration


for repeat=1:repeat_num
    
    % write text file
    fprintf(fileID,'%4d %05.2f %05.2f %09.2e %09.2e %09.2e\n',...
            r, betaH, betaW, sp, sm, max_iter);
    
end

fclose(fileID);
