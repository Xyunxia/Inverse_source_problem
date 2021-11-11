function [relative_err,lambda_opt]=auto_turing_lambda(f,J,ksi1,delta,beta,D)
g_delta_hat = generate_noise_measure(f,delta,J);
 % 采用ROF TV去噪文章中的调lambda的算法
lambda0 = 0.7079/delta + 0.6849/(delta^2);
lamd = zeros(1,10);
SpatialErr = zeros(1,10);
for tr = 1:10
  fm = oneD_TVBFD_apply(lambda0,beta,g_delta_hat,J,ksi1,delta,D);
  lamd(tr) = lambda0;
  SpatialErr(tr) = norm(fm-f) ;
  lambda0 = lambda0*SpatialErr(tr)/delta;
end
relative_err = SpatialErr/norm(f);
% 画图时剔除误差非常大的\lambda
index=find(relative_err>1.3);
SpatialErr(index) = [];
lamd(index) = [];
relative_err(index) = [];
if ~isempty(lamd)
    [a,index] = min(relative_err);
    plot(lamd,relative_err,'-o');
    lambda_opt = lamd(index);
    text(lambda_opt,a+0.02,['(',num2str(lambda_opt),',',num2str(a),')'],'color','r');
else    
    lambda_opt =0;
end
end