function [H] = NMF_proj(V,W,k)


% V: m*n
% W m*k
% H k*n
N_iter=1000;
rng(0,'twister');
[m,n]=size(V);
H=rand(k,n)+eps;

for i=1:N_iter
    % update H
    WV=W'*V;
    WWH=W'*W*H;
    H=H.*(WV./WWH);
    H=H./nansum(H,1); % Additional normalization step
    


end