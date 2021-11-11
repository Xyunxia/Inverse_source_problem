function [lambda_opt,b] = turing_lambda2D(f,J,X,Y,Ksi1,Ksi2,Ksi_norm2,delta,beta,flag,D1)
% purpose : given a noise level delta, from g_delta = g+delta*randn(N,1);
% generate noise measurement, then use TVBFD algorithm to solve the unknow
% source function f form the noise measurement g_delta by different lambda.
% Finally, compare the results and find the optimal parameter output.

%% generate the noise mesurement of the u_hat(ksi,T)
[g_delta_hat,g,g_delta] = generate_noise_measure(f,delta,J);
%% Manually set possible lambda values randomly
% 1e-6,0.5*1e-6,1e-5,0.5*1e-5,1e-4,0.5*1e-4,1e-3,0.5*1e-3, 二维问题lambda太小的话效果也不太好
%  
% small = [50:5:145];
% mid = [ 1050:250:1e+4];
% 500和800效果还不错
% small = [0.5*1e-2,0.01,0.5*1e-1,0.5,1,50:50:250 500:250:1e+3,1e+3,2.5*1e+3];
% mid = [5*1e+3:250:5*1e+4];
% big = [8*1e+4,1e+5,1.25*1e+5,2.5*1e+5,3*1e+5,5*1e+5,8*1e+5,1e+6,1e+7];
% lambda = [small mid big];

lambda0 = 0.7079/delta + 0.6849/(delta^2);
% first try:
lambda = [lambda0,5*1e+3,1.5*1e+3,2.5*1e+3,1.5*1e+4];

% sencond try:
% lambda = [lambda0,1e+3,2.5*1e+3,5*1e+3,1e+4,1.5*1e+4,2.5*1e+4,3.5*1e+4,...
%     4.5*1e+4,5*1e+4,5.15*1e+4,5.25*1e+4];

lambda = sort(lambda);  % 排序，从小到大排列
        % lamda取得过大可能导致Bregman迭代过程一步就满足终止条件了
si = size(lambda,2);
spatial_err = zeros(1,si);

for k = 1:si
    fm = TVBFD(g_delta_hat,J,Ksi1,Ksi2,Ksi_norm2,delta,...
        lambda(k),beta,D1);
    spatial_err(k) = norm(fm-f) ;
end
relative_err = spatial_err/norm(f);  % 相对误差是（近似值-真实值）/真实值
% 画图时剔除误差非常大的\lambda
index = find(relative_err>1.3);
spatial_err(index) = [];
lambda(index) = [];
relative_err(index) = [];
clear index;
if flag == 1
    source_type = 'Piecewise constant';
else
    source_type = 'continuous';
end
if isempty(lambda)==0
    % figure name里可以有tex字符
%     figure(1)
%     plot(lambda,spatial_err,'g.');
%     xlabel('\lambda');
%     ylabel('|| f_{TV}^m - f||');
%     figurenamestr = [source_type,',','\delta =',num2str(delta),',','use spatial err turing \lambda '];
%     title(figurenamestr);
%     hold on;
%     [a,index] = min(spatial_err);
%     lambda_opt = lambda(index);
%     plot(lambda_opt,a,'r*');
%     text(lambda_opt,a,['(',num2str(lambda_opt),',',num2str(a),')'],'color','r');
%     
%     % file name里不可以有tex字符
%     filename = [source_type,',','delta =',num2str(delta),',','use spatial err turing lambda.png'];
%     print(1,'-dpng',filename)   % 保存为png格式的图片

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
    text(lambda_opt,b,['(',num2str(lambda_opt),',',num2str(b),')'],'color','r');
    hold off
    filen = [source_type,',','delta =',num2str(delta),',','use relative err turing lambda.png'];
    print(2,'-dpng',filen)   % 保存为png格式的图片

    fm = TVBFD(g_delta_hat,J,Ksi1,Ksi2,Ksi_norm2,delta,...
            lambda_opt,beta,D1);
    %% plot 平面图
    figure(3)
    subplot(2,2,1)
    imagesc(f);
    colorbar;
    title('f(x,y)');
    subplot(2,2,2)
    imagesc(g);
    colorbar;
    title('g(x)=u(x,T)');
    subplot(2,2,3)
    imagesc(g_delta);
    colorbar;
    title(['g^{\delta}(x),\delta=',num2str(delta)])
    subplot(2,2,4)
    imagesc(fm);
    colorbar;
    title(['f_{TV}^{\delta}(x,y),RelErr =',num2str(b),'\lambda=',num2str(lambda_opt)]);
    finame = [source_type,',','delta =',num2str(delta),',', 'relative err = ',num2str(b),'2D','.png'];
    print(3,'-dpng',finame)   % 保存为png格式的图片
    %% 单个的图
    figure(5)
    mesh(X,Y,fm);
    xlabel('x');
    xlabel('y');
    xlabel('f_{TV}^m(x,y)');
    title(['delta = ',num2str(delta),'relative err = ',num2str(b),'.png'])
    finame = ([source_type,',','delta =',num2str(delta),',', 'relative err =',num2str(b),'.png'];
    print(5,'-dpng',finame)
    %% plot 3D图
    figure(4);
    subplot(2,2,1)
    % surf(X,Y,f);
    mesh(X,Y,f);
%     colorbar;
    title('f(x,y)');
    subplot(2,2,2)
    % surf(X,Y,g);
    mesh(X,Y,g);
%     colorbar;
    title('g(x)=u(x,T)');
    subplot(2,2,3)
    % surf(X,Y,g_delta);
    mesh(X,Y,g_delta);
%     colorbar;
    title(['g^{\delta}(x),\delta=',num2str(delta)])
    subplot(2,2,4)
    % surf(X,Y,fm);
    mesh(X,Y,fm);
%     colorbar;
    title(['f_{TV}^{\delta}(x,y),RelErr =',num2str(b),'\lambda=',num2str(lambda_opt)]);

    finalname=[source_type,',','delta =',num2str(delta),',', 'relative err = ',num2str(b),',3D.png'];
    % title(finalname);
    print(4,'-dpng',finalname)   % 保存为png格式的图片
else
    fprintf('The code has err!');
    lambda_opt = 0;b=0;
end
end