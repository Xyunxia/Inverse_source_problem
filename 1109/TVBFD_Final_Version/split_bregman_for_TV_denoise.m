function [uk,uk_hat] = split_bregman_for_TV_denoise(u0,P,Ksi1,Ksi2,Ksi_norm2,gm_hat,lambda,beta,D1)
% purpose: 
% split bregman iteration for solve  the optimal problem
% $$
% u^{m+1}=\arg \min J_{\mathrm{TV}}(u)+\frac{\lambda}{2}\left\|P^{-1} \hat{u}-\hat{g}^{m}\right\|^{2}
% $$

eps = 5*1e-6;
[width,height] = size(u0);
d_tidle_x = zeros(width,height);
d_tidle_y = zeros(width,height);
dx = zeros(width,height);
dy = zeros(width,height);
% persistent C;  % 静态变量，函数调用结束后这个值仍然存在，可供下一次调用
% if isempty(C)  % 这里不适合用静态变量，因为lambda是在变的
   C = (lambda*P.^2 + beta * Ksi_norm2 ); 
% end

% err = [];
% err = [err norm(P .* u0_hat - gm_hat)];
for k = 0 : 5*1e+4
    %% solve d-subproblem
%     vector_shrink(u0,width,height,beta);
[d_tidle_x_new,d_tidle_y_new,dx_new,dy_new]...
    =solveDsubproblem(d_tidle_x,d_tidle_y,dx,dy,u0,beta,D1);
    %% solve u-subproblem
    w1_hat = fft2(d_tidle_x_new);
    w2_hat = fft2(d_tidle_y_new);
    temp = Ksi1.*w1_hat + Ksi2.*w2_hat;
    uk_hat = (lambda * P.*gm_hat - 1i*beta * temp)./C;
    uk_hat(isnan(uk_hat)) = 0;
    uk_hat(isinf(uk_hat)) = 0;
    uk = abs(ifft2(uk_hat)); 
    %% caculate err
    normu = norm(uk);
    if normu>1e-3
        diffNorm = norm(uk-u0)/normu;
    else
        fprintf('Bregman iteration has err!');
        break; 
    end
    
    if(diffNorm <eps ) 
        break;
    end
    %% update
    u0 = uk;
    dx = dx_new;
    dy = dy_new;
    d_tidle_x = d_tidle_x_new;
    d_tidle_y = d_tidle_y_new;
end
% fprintf('split Bregman iter = %d , u relative err = %0.9g \n',k,diffNorm);
end