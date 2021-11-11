function [lambda_opt,piecewise_lambda_opt]=test_TVBFD(delta_vec,piecewise_delta_vec,flag)
% purpose: 
% Using the TVBFD algorithm to solve the inverse source problem of
%    four different source function.
% Input : 
%  delta_vec :  a vector for diffrent noise level. 
%     For one-dimensional problems, this parameter is used to generate 
%     noise observations for smooth source terms; 
%     for two-dimensional problems, this parameter is used to generate 
%     noise observations for continues source terms.
%  piecewise_delta_vec: a vector for diffrent noise level.
%     For one-dimensional problems, this parameter is used to generate 
%     noise observations for piecewise linear source terms; 
%     for two-dimensional problems, this parameter is used to generate 
%     noise observations for piecewise constent source terms.
%  flag: 
%      Determine whether it's a one-dimensional problem or 
%      a two-dimensional problem.
%      flag = 1; one-dimensional problem.
%      flag = 2; two-dimensional problem.
% usage:
% delta_vec = [0.5*1e-2,1e-2,5*1e-2,0.1];
% piecewise_delta_vec = [0.5*1e-3,1e-3,0.5*1e-2,1e-2,5*1e-2,0.1];
% [lambda_opt,piecewise_lambda_opt]=test_TVBFD(delta_vec,piecewise_delta_vec,1)
% [lambda_opt,piecewise_lambda_opt]=test_TVBFD(delta_vec,piecewise_delta_vec,2)
%% 1D problem
if 1==flag
    % spatial [-15,15]
    N = 3000; 
    N1 = N/2;
    L = 15;
    h = 2*L/N;
    x = h*[-N1:N1-1]';

    % frequency:[-pi*N/(2L),pi*N/(2L)]
    ksi1 = pi/L*[0:N1-1 -N1:-1]'; % 这是为了避免做fft之后要用fftshift
    %% caculate intergral
    % J = Compute_Integral_term(ksi1,N);
    load('N3000_1D_integral.mat'); 

    %% define differential matrix
    e = ones(N,1);
    D = spdiags([-e,e],[0,1],N,N);
    D(N,N) = 0;
    %% smooth source function
    F = @(x) 5*exp(-(x+5).^2)+3*exp(-2*x.^2)+exp(-(x-4).^2);
    % Flinear_hat = @(x) 5.0/sqrt(2)*exp(-x.^2/4+5i*x)...
    %     +3.0/2.0*exp(-x.^2/8)+1.0/sqrt(2)*exp(-x.^2/4-4i*x);

    f = feval(F,x);

    delta_num = size(delta_vec,2);

    % TVBFD start
    beta=5;    
    % turing_lambda(f,J,x,ksi1,1e-2,beta,flag,D);  % 观察迭代过程时取消注释
    lambda_opt = zeros(1,delta_num);
    relative_err_opt = zeros(1,delta_num);
    broadcast = J;
    parfor k = 1:delta_num
    % lambda and beta are the htper-parameter waiting for turing
        [lambda_opt(k),relative_err_opt(k)] = turing_lambda(f,broadcast,x,ksi1,delta_vec(k),beta,0,D);
    end

    fid = fopen('TVBFD_result.txt','a'); % 读写方式打开文件
    fprintf(fid,'\n%s\n','smooth source function:');
    fprintf(fid,'%s\t','delta:');
    fprintf(fid,'%g\t%g\t%g\t%g\n',delta_vec);
    fprintf(fid,'%s\t','optimal lambda:');
    fprintf(fid,'%g\t%g\t%g\t%g\n',lambda_opt);
    fprintf(fid,'%s\t','optimal relative err:');
    fprintf(fid,'%g\t%g\t%g\t%g\n',relative_err_opt);

    % [relative_err,lambda_opt]=auto_turing_lambda(f,J,ksi1,delta_vec(2),beta,D);
    %% piecewise_linear_source function
    f_piecewise_linear = 2*piecewise_linear_source(x,-4,1)+...
        piecewise_linear_source(x,0,2)+3*piecewise_linear_source(x,4,1);

    
    piecewise_delta_num = size(piecewise_delta_vec,2);

    %% 调参数lambda
    piecewise_lambda_opt = zeros(1,piecewise_delta_num);
    piecewise_relative_err_opt = zeros(1,piecewise_delta_num);
    % lambda and beta are the htper-parameter waiting for turing

    parfor k = 1:piecewise_delta_num
        [piecewise_lambda_opt(k),piecewise_relative_err_opt(k)] = ...
            turing_lambda(f_piecewise_linear,broadcast,x,ksi1,piecewise_delta_vec(k),beta,1,D);
    end
    fid = fopen('TVBFD_result.txt','a'); % 读写方式打开文件
    fprintf(fid,'\n%s\n','piecewise linear source function:');
    fprintf(fid,'%s\t','delta:');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_delta_vec);
    fprintf(fid,'%s\t','optimal  lambda:');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_lambda_opt);
    fprintf(fid,'%s\t','optimal relative err:');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_relative_err_opt);

    % delta = 0.05;
    % g_delta_hat = generate_noise_measure(f_piecewise_linear,delta,J);
    % fm = oneD_TVBFD_apply(1.5*1e+4,beta,g_delta_hat,J,ksi1,delta,D);
    % b = norm(fm - f_piecewise_linear)/norm(f_piecewise_linear);
    % figure(3)
    % plot(x,fm,'b:',x,f_piecewise_linear,'r--','linewidth',1.5);
    % xlabel('x');
    % ylabel('f_{TV}(x)')
    % legend('f^m','f exact');
    % finalname=['Piecewise linear',',','\delta =',num2str(delta),',', 'relative err = ',num2str(b)];
    % title(finalname);
    % finame = ['Piecewise linear',',','delta =',num2str(delta),',', 'relative err = ',num2str(b),'.png'];
    % print(3,'-dpng',finame)   % 保存为png格式的图片

%%  2D problem
else
    %% discretization
    % spatial [-10,10]
    N = 500;  N1 = N/2;
    L = 10;
    h = 2*L/N;
    x = h*[-N1:N1-1]';
    [X,Y] = meshgrid(x,x);

    % frequency:[-pi*N/(2L),pi*N/(2L)]
    ksi1 = pi/L*[0:N1-1 -N1:-1]'; % 这是为了避免做fft之后要用fftshift
    [Ksi1,Ksi2] = meshgrid(ksi1,ksi1);
    Ksinorm2 = Ksi1.^2 + Ksi2.^2;

    %% caculate intergral
    load('N500_2D_integral.mat');
    % Ksinorm = sqrt(Ksinorm2);
    % % 将Ksinorm拉直为列向量，方便计算积分值
    % Ksi = Ksinorm(:);
    % J = Compute_Integral_term(Ksi,N);
    %% define differential matrix
    % differial matrix such that D2*u = grad_x u, u*D1' = grad_y u;
    e = ones(N,1);
    D1 = spdiags([-e,e],[0,1],N,N);
    D1(N,N) = 0;

    % D2 = spdiags([-e,e],[0,1],N,N); % 默认输入是方阵，所以可注释掉
    % D2(N,N) = 0;
    %% continuous source function
    F = @(x,y) 2*exp(-(x+5).^2-y.^2)+5*exp(-(x-4).^2-y.^2);
    % F_hat = @(x,y) 2.0/sqrt(2)*exp(-x.^2/4+5i*x).*exp(-y.^2/4)...
    %     +5.0/sqrt(2)*exp(-x.^2/4-4i*x).*exp(-y.^2/4);
    f = feval(F,X,Y);

    % delta_vec = [0.5*1e-2,1e-2,5*1e-2,0.1,0.2];
    delta_num = size(delta_vec,2);

    %% use TVBFD
    beta = 20;  % 二维β用20效果比较好
    lambda_opt = zeros(1,delta_num);
    relative_err_opt = zeros(1,delta_num);
    broadcast = J;
    parfor k = 1:delta_num
    % lambda and beta are the htper-parameter waiting for turing
        [lambda_opt(k),relative_err_opt(k)] = turing_lambda2D(f,broadcast,X,Y,...
            Ksi1,Ksi2,Ksinorm2,delta_vec(k),beta,0,D1);
    end
    fid = fopen('TVBFD_result.txt','a'); % 读写方式打开文件
    fprintf(fid,'\n%s\n','continuous source function:');
    fprintf(fid,'%s\t','delta:');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\n',delta_vec);
    fprintf(fid,'%s\t','optimal lambda:');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\n',lambda_opt);
    fprintf(fid,'%s\t','optimal relative err:');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\n',relative_err_opt);

    %% challenging example
    X1 = X(:);Y1 = Y(:);
    temp = piecewise_constant_source(X1,Y1);
    f_piecewise_constant = reshape(temp,N,N);

    % piecewise_delta_vec = [0.5*1e-3,1e-3,0.5*1e-2,1e-2,5*1e-2,0.1,0.2];
    piecewise_delta_num = size(piecewise_delta_vec,2);

    piecewise_lambda_opt = zeros(1,piecewise_delta_num);
    piecewise_relative_err_opt = zeros(1,piecewise_delta_num);
    parfor k = 1:piecewise_delta_num
    % lambda and beta are the htper-parameter waiting for turing
        [piecewise_lambda_opt(k),piecewise_relative_err_opt(k)] = ...
            turing_lambda2D(f_piecewise_constant,broadcast,X,Y,Ksi1,Ksi2,Ksinorm2,...
            piecewise_delta_vec(k),beta,1,D1);
    end
    fid = fopen('TVBFD_result.txt','a'); % 读写方式打开文件
    fprintf(fid,'%s\n','piecewise constant source function:');
    fprintf(fid,'%s\t','delta:');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_delta_vec);
    fprintf(fid,'%s\t','optimal lambda:');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_lambda_opt);
    fprintf(fid,'%s\t','optimal relative err:');
    fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_relative_err_opt);
end
end