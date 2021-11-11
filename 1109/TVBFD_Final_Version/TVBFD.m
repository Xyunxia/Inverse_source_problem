function [um,um_hat,TErr] = TVBFD(g0_hat,P,Ksi1,Ksi2,Ksi_norm2,delta,lambda,beta,D1)
% Popusure: use Bregman iteration solve the TV model
% \underset{f \in \mathcal{A}}{\arg \min } ~J_{\mathrm{TV}}\left(f\right)+
% \frac{\lambda}{2}\left\|P\cdot \widehat{f}(\xi)-\widehat{g}^{\delta}\right\|_{L^{2}
%     \left(\mathbb{R}^{2}\right)}^{2}\label{gs11}

[width,height] = size(g0_hat);
% um = zeros(width,height);
u0_hat = zeros(width,height);
u0 = zeros(width,height);
gm_hat = g0_hat;

tau = 1.01*44.3925; 
% 给tau乘44.3925是因为用 g_delta = g + delta*randn(r,w);生成的含噪观测，
% 其误差的范数||g_delta - g|| 约等于 44.3925*delta;

% TErr = zeros(80,1);

m = 1;
Cd = sqrt(width*height); % 二维保范离散傅里叶变换的范数与fft2范数之间相差的常数
tic
while(norm(P.*u0_hat - g0_hat)/Cd >= tau*delta)
% while(1)   % 展示迭代效果时用
    [um,um_hat] = split_bregman_for_TV_denoise(u0,P,Ksi1,Ksi2,Ksi_norm2,gm_hat,lambda,beta,D1);
    
   % TErr(m) = norm(P.*um_hat - g0_hat)/Cd;

    if(m > 80 )  
        break;
    end
    %% update
    gm_hat = gm_hat + (g0_hat - P .* um_hat) ; 
    u0_hat = um_hat;
    u0 = um;
    m = m + 1;
end
toc
% index = (TErr==0);
% TErr(index) =[]; 
% fprintf('iter = %d \n',m);
% TErr
%     imagesc(real(um));
%     imshow(um);
%     pause(1);
end