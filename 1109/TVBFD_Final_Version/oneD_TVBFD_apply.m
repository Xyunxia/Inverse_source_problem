function [um,um_hat,Err] = oneD_TVBFD_apply(lambda,beta,g_delta_hat,P,ksi1,delta,D)
% purpose: solve the TV model
% \underset{f \in \mathcal{A}}{\arg \min } ~J_{\mathrm{TV}}\left(f\right)+
% \frac{\lambda}{2}\left\|P\cdot \widehat{f}(\xi)-\widehat{g}^{\delta}\right\|_{L^{2}
%     \left(\mathbb{R}^{2}\right)}^{2}\label{gs11}

[width,height] = size(g_delta_hat);

% initialize
gm_hat = g_delta_hat;
um = zeros(width,height);
um_hat = zeros(width,height);
tau = 1.001*55.07;   
% 给tau乘55.07是因为用 g_delta = g+delta*randn(r,w);生成的含噪观测，
% 其误差的范数||g_delta - g|| 约等于 55.07*delta;
m=1; 
eps = 5*1e-6;
C = (lambda*P.^2 + beta * ksi1.^2 ); 
thresh = 1.0/beta;
Err = zeros(80,1);

% ferr = zeros(80,1);% 测试迭代过程时用

 tic
 % fft不保范，它与保范离散傅里叶变换之间相差一个常数sqrt(N)
while(norm(P.*um_hat - g_delta_hat)/sqrt(width) >= tau*delta)
% while(1)  % 测试迭代过程时用
    um0 = um;
   
    %% split bregman method
    d_tidle_x = zeros(width,height);
    dx = zeros(width,height);
    for iter =0:1:5*1e+4
        %% solve d-subproblem
        grad_x = D*um0;
        dx_new = zeros(width,height);
        dx(1:end-1)=dx(1:end-1)+ grad_x(1:width-1)- d_tidle_x(1:width-1);   
        abs_dx = abs(dx(1:end-1));
        index1 = find(abs_dx>1e-5);
        index2 = find(abs_dx<1e-5);
        temp = max(abs_dx - thresh ,0);
        if ~isempty(index1)
            dx_new(index1) = temp(index1).*dx(index1)./abs_dx(index1);
        end
        if ~isempty(index2)
            dx_new(index2) = 0;
        end
       % update dtidle 
        d_tidle_x(1:width-1) = 2*dx_new(1:width-1) - dx(1:width-1);
        dx = dx_new;
       clear abs_dx;clear temp;clear index1;clear index2;
        %% solve u-subproblem
        w_x_hat = fft(d_tidle_x);
        umk_hat = (lambda * P.*gm_hat - 1i*beta*ksi1.*w_x_hat)./C;
        umk_hat(isnan(umk_hat)) = 0;
        umk_hat(isinf(umk_hat)) = 0;
        umk = abs(ifft(umk_hat));
        
        udiffer = norm(umk-um0);
        umnorm = norm(umk);
        %% caculate err
        if umnorm > 1e-3
            diffNorm = udiffer/umnorm;
        else
            diffNorm = 0;
            printf('Bregman iteration has err!');
            break; 
        end
       
        if(diffNorm <eps )  % || err(end)>err(end-1)
            break;
        end     
        um0 = umk;
    end
    
   fprintf('Bregman iter = %d, u-subproblem err = %0.9g\n',iter,diffNorm);    
%    ferr(m) =  norm(f - umk);
   Err(m) =  norm(P.*umk_hat - g_delta_hat)/sqrt(width);
   
   if(m > 80 ) 
        break;
   end
    %% update
    um = umk;
    um_hat = umk_hat;
    gm_hat = gm_hat + (g_delta_hat - P .* um_hat) ; 
    m = m + 1;
end
toc
index = (Err==0);
Err(index) =[]; 
Err
%  figure;
%  plot(1:size(Err,2),Err);
%  xlabel('m');
%  ylabel('|| I\wedge{f}^m - {\wedge}{g}^{\delta}||');
%  hold on;
%  plot([1,size(Err,2)],[tau*delta tau*delta],'r-');
%  title(['\delta =',num2str(delta),'\lambda = ',num2str(lambda)])
%  figure;
%  plot(1:size(ferr,2),ferr);
%  xlabel('m');
%  ylabel('|| fm - f||');
%  title(['\delta =',num2str(delta),'\lambda = ',num2str(lambda)])
end