function [prio_relative_err,poster_relative_err] = direct_solve_method(delta,f,J,Ksi,N1,X,flag)
% purpose: 
% relazation the direct method for solving the inverse source
% problem, which is proposed by Nguyen Van Thang in the artical
% Identifying an unknown source term in a time-space fractional parabolic equation

% Input：
% f : the value of exact source function in the discrete mesh point
% in one dimensional case, f must be a vector;
% in two dimensional case, f must be a matrix;
% delta : noise level
% J: Numerical approximation of integral values associated with ...
%    the frequency domain grid, which satisfied that 
%    \hat{f} = \hat{g}./J;
%    the size of J is equal to the size of J
% X: the discrete mesh grid in spacial, for One Dimensional case, X is a vector;
%    for two Dimensional case, X is a matrix, which is first col is
%    the x coordinate , the second col is the y coordinate.
% Ksi: the discrete mesh grid in frequency .
% flag: indicate the type of source function, has 4 different values
% flag = 0   1D smooth source;
% flag = 1   1D piecewise linear source
% flag = 2   2D continues source;
% flag = 3   2D piecewise constent source
[r,w] = size(f);
E = 120;  
p=3;
gamma=1;

v = (E/delta)^(1/(p+gamma));

if w == 1
    %% 1D
     %% generate the noise mesurement of the u_hat(ksi,T);
    g_delta_hat = generate_noise_measure(f,delta,J);
%     f_hat = fft(f);
%     g_hat = f_hat.*J;
%     g = abs(ifft(g_hat));
%      % method 1 : norm(eta) = delta
%     g_delta = g + delta*randn(r,w)/55.07;  
%     tau = 1.001;  

%    % method 2 : norm(eta) = delta*55.07
%     g_delta = g + delta*randn(r,w); % norm(eta) = delta*55.07
%    g_delta_hat = fft(g_delta);
    tau = 1.001*55.07; 
    %% use a priori parameter choice rule
    xvalid = find(abs(Ksi)<v);
    Sv_hat = zeros(r,w);
    Sv_hat(xvalid) = g_delta_hat(xvalid);
    fv_hat = Sv_hat./J;
    fv_prio = abs(ifft(fv_hat));
    
    %% use a posteriori parameter choice rule
    % G(v_{\delta}) = ||S_{v_{\delta}}-g^{\delta}||_{L_2} = \tau \delta
    Sv_opt_hat = dichotomizing_search_the_posteriori_para(g_delta_hat,delta,Ksi,N1,flag,tau) ;
    fv_poster_hat = Sv_opt_hat./J;
    fv_poster = abs(ifft(fv_poster_hat));
else
    %% 2D
    %% generate the noise mesurement of the u_hat(ksi,T);
    [g_delta_hat,~,g_delta] = generate_noise_measure(f,delta,J);
%         f_hat = fft2(f);
%     g_hat = f_hat.*J;
%     g = abs(ifft2(g_hat));
%      % method 1 : norm(eta) = delta
%     g_delta = g + delta*randn(r,w)/44.2725;
%     tau = 1.01;  %*44.2725

%     % method 2 : norm(eta) = delta*44.2725
%     g_delta = g + delta*randn(r,w);
%     g_delta_hat = fft2(g_delta);
    tau = 1.01*44.2725;
    %% use a priori parameter choice rule
%     [Ksi1,Ksi2] = meshgrid(ksi1,ksi1);
    temp1 = Ksi(:,1);temp2 = Ksi(:,2);
    ind1 = find(abs(temp1)<v);
    ind2 = find(abs(temp2)<v);
    index = intersect(ind1,ind2'); % 找两个数组共同包含的元素
    vector_g_delta_hat = g_delta_hat(:);
    Valid_temp = zeros(r*w,1);
    Valid_temp(index) = vector_g_delta_hat(index);
    Sv_hat = reshape(Valid_temp,r,w);
    fv_hat = Sv_hat./J;
    fv_prio = abs(ifft2(fv_hat));
    
      %% use a posteriori parameter choice rule
    % G(v_{\delta}) = ||S_{v_{\delta}}-g^{\delta}||_{L_2} = \tau \delta
    Sv_opt_hat = dichotomizing_search_the_posteriori_para(g_delta_hat,delta,Ksi,N1,flag,tau) ;
    fv_poster_hat = Sv_opt_hat./J;
    fv_poster = abs(ifft2(fv_poster_hat));
end

%% calculate err
err_prio = norm(fv_prio-f);
normf = norm(f);
prio_relative_err = err_prio/normf;   
err_poster = norm(fv_poster - f);
poster_relative_err = err_poster/normf;
%% plot
if flag == 0   
    type = 'smooth' ;
elseif flag == 1  
    type = 'piecewise linear'; 
elseif flag == 2
    type = 'continues' ;
elseif flag == 3 
    type = 'piecewise constent';
else
    fprintf('Please pass in the correct parameter value!')
end
if flag ==0 || flag == 1  
    figure(5)
    plot(X,f,'r--',X,fv_prio,'b:','linewidth',1.5);
    legend('Exact','Reconstructed');
    xlabel('x');
    ylabel('f^{v}(x)');
    title([type,',','\delta =',num2str(delta),',','priori relative err =',num2str(prio_relative_err)]);
    filename = [type,',','delta =',num2str(delta),',','priori solve relative err =',num2str(prio_relative_err),'.png'];
    print(5,'-dpng',filename);
    figure(6)
    plot(X,f,'r--',X,fv_poster,'b:','linewidth',1.5);
    legend('Exact','Reconstructed');
    xlabel('x');
    ylabel('f^{v}(x)');
    title([type,',','\delta =',num2str(delta),',','poster relative err =',num2str(poster_relative_err)]);
    filename = [type,',','delta =',num2str(delta),',','poster solve relative err =',num2str(poster_relative_err),'.png'];
    print(6,'-dpng',filename);
 % 出来的图是高瘦的，不好看
%     figure(8)
%     subplot(1,2,1)
%     plot(X,f,'r--',X,fv_prio,'b:');
% %     axis( [min(X) max(X)]) 
%     legend('Exact','Reconstructed');
%     xlabel('x');
%     ylabel('f^{v}(x)');
%     title(['\delta =',num2str(delta),',','priori ',num2str(prio_relative_err)]);
%     subplot(1,2,2)
%     plot(X,f,'r--',X,fv_poster,'b:');
%     legend('Exact','Reconstructed');
%     xlabel('x');
%     ylabel('f^{v}(x)');
%     title(['\delta =',num2str(delta),',','poster ',num2str(poster_relative_err)]);
%     filename = [type,',','delta =',num2str(delta),',','direct2.png'];
%      print(8,'-dpng',filename);
    %% 两个图画在一起，有点难看
%     figure(7)
%     ha = tight_subplot(1, 2, 0.05);
%     axes(ha(1));
%     plot(X,f,'r--',X,fv_prio,'b:');
%     legend('Exact','Reconstructed');
%     xlabel('x');
%     ylabel('f^{v}(x)');
%     title(['\delta =',num2str(delta),',','priori ',num2str(prio_relative_err)]);
%     
%     axes(ha(2));
%     plot(X,f,'r--',X,fv_poster,'b:');
%     legend('Exact','Reconstructed');
%     xlabel('x');
%     ylabel('f^{v}(x)');
%     title(['\delta =',num2str(delta),',','poster ',num2str(poster_relative_err)]);
%     filename = [type,',','delta =',num2str(delta),',','direct solve.png'];
%      print(7,'-dpng',filename);
else
    
    X1 = X(:,1);X2 = X(:,2);
    X1 = reshape(X1,2*N1,2*N1);X2 = reshape(X2,2*N1,2*N1);

    figure(5)
    mesh(X1,X2,fv_prio);
%     colorbar;
    xlabel('x');
    ylabel('y');
    zlabel('f^{v}(x)');
    title([type,',','\delta =',num2str(delta),',','priori relative err =',num2str(prio_relative_err)]);
    filename = [type,',','delta =',num2str(delta),',','priori solve relative err =',num2str(prio_relative_err),'.png'];
    print(5,'-dpng',filename);
   
    figure(6)
    mesh(X1,X2,fv_poster);
%     colorbar;
    xlabel('x');
    ylabel('y');
    zlabel('f^{v}(x)');
    title([type,',','\delta =',num2str(delta),',','poster relative err =',num2str(poster_relative_err)]);
    filename = [type,',','delta =',num2str(delta),',','poster solve relative err =',num2str(poster_relative_err),'.png'];
    print(6,'-dpng',filename);
    
    %% plot 平面图
    figure(7)
    subplot(2,2,1)
    imagesc(f);
     colorbar;
    title('f(x,y)');
    subplot(2,2,2)
    imagesc(g_delta);
     colorbar;
    title(['g^{\delta}(x,y),\delta=',num2str(delta)])
    
    subplot(2,2,3)
    imagesc(fv_prio);
     colorbar;
    title(['priori rule,',num2str(prio_relative_err)]);
    subplot(2,2,4)
    imagesc(fv_poster);
     colorbar;
    title(['poster rule,',num2str(poster_relative_err)]);
    finame = [type,',','delta =',num2str(delta),',','direct2D','.png'];
    print(7,'-dpng',finame)   % 保存为png格式的图片
    
    %% plot 3D图
    figure(9)
    subplot(2,2,1)
     mesh(X1,X2,f);
%      colorbar;
    title('f(x,y)');
    subplot(2,2,2)
    mesh(X1,X2,g_delta);
%     colorbar;
    title(['g^{\delta}(x,y),\delta=',num2str(delta)])
    
    subplot(2,2,3)
    mesh(X1,X2,fv_prio);
%     colorbar;
    title(['priori rule,',num2str(prio_relative_err)]);
    subplot(2,2,4)
    mesh(X1,X2,fv_poster);
%     colorbar;
    title(['poster rule,',num2str(poster_relative_err)]);
    finame = [type,',','delta =',num2str(delta),',','direct3D','.png'];
    print(9,'-dpng',finame)   % 保存为png格式的图片
end

end