function  [lambda_opt,b] = turing_lambda(f,J,x,ksi1,delta,beta,flag,D)
% purpose : given a noise level delta, from g_delta = g+delta*randn(N,1);
% generate noise measurement, then use TVBFD algorithm to solve the unknow
% source function f form the noise measurement g_delta

g_delta_hat = generate_noise_measure(f,delta,J);
%% Manually set possible lambda values randomly
% 1e-6,0.5*1e-6,1e-5,0.5*1e-5,1e-4,0.5*1e-4,1e-3,0.5*1e-3,0.5*1e-2,...
% 0.01,0.5*1e-1,0.5,1  too small is bad!
% first try:
% lambda = [5,50,150,250,500,1e+3,2.5*1e+3,5*1e+3,1e+4,1.5*1e+4];
% find that the opt lambda is betwen 5*1e+3~1e+4

% sencond try:
lambda0 = 0.7079/delta + 0.6849/(delta^2);
lambda = [lambda0,1e+3,2.5*1e+3,5*1e+3,1e+4,1.5*1e+4,2.5*1e+4,3.5*1e+4,...
    4.5*1e+4,5*1e+4,5.15*1e+4,5.25*1e+4];
    % ,1e+5,1.25*1e+5,2.5*1e+5];
% small = [5:15:150];
% mid = [155:25:475 500:125:1e+3];
% big = [1050:250:1e+4 1.25*1e+4,2.5*1e+4,5*1e+4,8*1e+4,1e+5,1.25*1e+5,...
%     2.5*1e+5,5*1e+5,8*1e+5,1e+6,2.5*1e+6,5*1e+6,8*1e+6,1e+7];
% lambda = [small mid big];
lambda = sort(lambda);  % 排序，从小到大排列
        % lamda取得过大可能导致Bregman迭代过程一步就满足终止条件了 
si = size(lambda,2);
spatial_err = zeros(1,si);

for b = 1:si
    fm = oneD_TVBFD_apply(lambda(b),beta,g_delta_hat,J,ksi1,delta,D);
    spatial_err(b) = norm(fm-f) ;
end
relative_err = spatial_err/norm(f);  % 相对误差是（近似值-真实值）/真实值
% 画图时剔除误差非常大的\lambda
index=find(relative_err>1.3);
spatial_err(index) = [];
lambda(index) = [];
relative_err(index) = [];
if flag == 1
    source_type = 'Piecewise linear';
else
    source_type = 'smooth';
end
% figure name里可以有tex字符
% figurenamestr = [source_type,',','\delta =',num2str(delta),',','use spatial err turing \lambda '];
% figure(1)
% plot(lambda,spatial_err,'g.');
% xlabel('\lambda');
% ylabel('|| f_{TV}^m - f||');
% title(figurenamestr);
% hold on;
% [a,index] = min(spatial_err);
% lambda_opt = lambda(index);
% plot(lambda_opt,a,'r*');
% text(lambda_opt,a+0.02,['(',num2str(lambda_opt),',',num2str(a),')'],'color','r');
% hold off
% % file name里不可以有tex字符
% filename = [source_type,',','delta =',num2str(delta),',','use spatial err turing lambda.png'];
% print(1,'-dpng',filename)   % 保存为png格式的图片

figure(2)
plot(lambda,relative_err,'g.');
xlabel('\lambda');
ylabel('|| f_{TV}^m - f||/||f||');
titlename = [source_type,',','\delta =',num2str(delta),',','use relative err turing \lambda '];
title(titlename);
hold on;
[b,index] = min(relative_err);
lambda_opt = lambda(index);
plot(lambda_opt,b,'r*');
text(lambda_opt,b+0.02,['(',num2str(lambda_opt),',',num2str(b),')'],'color','r');
hold off
filen = [source_type,',','delta =',num2str(delta),',','use relative err turing lambda.png'];
print(2,'-dpng',filen)   % 保存为png格式的图片

fm = oneD_TVBFD_apply(lambda_opt,beta,g_delta_hat,J,ksi1,delta,D);
figure(3)
plot(x,fm,'b:',x,f,'r--','linewidth',1.5);
xlabel('x');
ylabel('f(x)');
legend('f^m','f exact');
finalname=[source_type,',','\delta =',num2str(delta),',', 'relative err = ',num2str(b)];
title(finalname);
finame = [source_type,',','delta =',num2str(delta),',', 'relative err = ',num2str(b),'.png'];
print(3,'-dpng',finame)   % 保存为png格式的图片
%% 观察迭代过程中的变化
% 需要对oneD_TVBFD_apply的输入参数进行改变
%% best \lambda and max iter M
% [~,~,Err1,ferr1] = oneD_TVBFD_apply(f,2.5*1e+4,beta,g_delta_hat,J,ksi1,delta,D);
% [~,~,Err2,ferr2] = oneD_TVBFD_apply(f,1e+5,beta,g_delta_hat,J,ksi1,delta,D);
% [~,~,Err3,ferr3] = oneD_TVBFD_apply(f,5*1e+3,beta,g_delta_hat,J,ksi1,delta,D);
% [a1,b1]=min(Err1);
% [a2,b2]= min(Err2);
% [a3,b3]=min(Err3);
% tau = 1.001*55.07;  
% figure(2);
% xl = 1:size(Err1,1);
%  plot(xl,Err1,'g.',xl,Err2,'b.',xl,Err3,'k.');
%  xlabel('m');
%  ylabel('|| I\wedge{f}^m - {\wedge}{g}^{\delta}||');
%  
%  hold on;
%  plot([1,size(Err1,1)],[tau*delta tau*delta],'r-');
%  legend('\lambda = 2.5*1e+4','\lambda = 1e+5','\lambda = 5*1e+3','\tau\delta');
%  plot([xl(b1) xl(b2) xl(b3)],[a1 a2 a3],'ro');
% text(xl(b1), a1, ['(',num2str(xl(b1)),',',num2str(a1),')'],'color','g');
% text(xl(b2), a2, ['(',num2str(xl(b2)),',',num2str(a2),')'],'color','b');
% text(xl(b3), a3, ['(',num2str(xl(b3)),',',num2str(a3),')'],'color','k');
%  title(['\delta =',num2str(delta)])
%  hold off;
%  finame = ['delta =',num2str(delta),',iterF.png'];
%  print(2,'-dpng',finame)  
%  fx = 1:size(ferr1,1);
%  
%  [a1,b1]=min(ferr1);
% [a2,b2]= min(ferr2);
% [a3,b3]=min(ferr3);
%  figure(3);
%  plot(fx,ferr1,'g.',fx,ferr2,'b.',fx,ferr3,'k.');
%  xlabel('m');
%  ylabel('|| fm - f||');
%  title(['\delta =',num2str(delta)])
%  legend('\lambda = 2.5*1e+4','\lambda = 1e+5','\lambda = 5*1e+3');
%   hold on;
%  plot([xl(b1) xl(b2) xl(b3)],[a1 a2 a3],'ro');
% text(fx(b1), a1, ['(',num2str(fx(b1)),',',num2str(a1),')'],'color','g');
% text(fx(b2), a2, ['(',num2str(fx(b2)),',',num2str(a2),')'],'color','b');
% text(fx(b3), a3, ['(',num2str(fx(b3)),',',num2str(a3),')'],'color','k');
% hold off;
%  finame = ['delta =',num2str(delta),',iterS.png'];
%  print(3,'-dpng',finame)  
end