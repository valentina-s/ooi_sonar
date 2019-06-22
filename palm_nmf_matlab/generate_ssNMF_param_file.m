% Test writing and reaiding formatted txt files

save_path = '~/Downloads/';
save_file = 'test.txt';

rank = 3;          % rank
betaH = 0.1;
betaW = 0.1;
sp = [2,5,10,20];  % sparsity (lambda)

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
fprintf(fileID,'%4s %5s %5s %9s %9s\n','rank','betaH','betaW','sparsity','smoothness');
for isp=1:length(sp)
    for ism=1:length(sm)
        fprintf(fileID,'%4d %05.2f %05.2f %09.2e %09.2e\n',rank, betaH, betaW, sp(isp), sm(ism));
    end
end
fclose(fileID);

% % read text file
% fileID = fopen(fullfile(save_path,save_file));
% C = textscan(fileID,'%d %f %f %f %f','HeaderLines',1);
% fclose(fileID);
