function seq = find_match_factor_seq(rho,rank)
% Find matching component sequence based on correlation coefficients
%
% Inputs
%   rho   correlation coefficients (based on H or W, H seems to give better results)
%   rank  rank of the decomposition
%
% Outputs
%   seq   matching sequence
%
% Wu-Jung Lee | leewujung@gmail.com
% 2019 05 03

icomp = 1;
while icomp<=rank
    if icomp==1
        seq = nan(rank,1);
    end
    [~,max_ind] = max(rho(:));
    [x_ind,y_ind] = ind2sub(size(rho),max_ind);
    seq(x_ind) = y_ind;
    rho(x_ind,:) = ones(size(rho,2),1)*-2;
    rho(:,y_ind) = ones(1,size(rho,1))*-2;
    icomp = icomp+1;
end
