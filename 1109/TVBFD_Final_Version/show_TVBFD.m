function show_TVBFD
%(delta_vec,lambda_opt,piecewise_delta_vec,piecewise_lambda_opt,flag)
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
% usage:
%     Just click the Run button
% 
%% 1D problem
% % spatial [-15,15]
% N = 3000; 
% N1 = N/2;
% L = 15;
% h = 2*L/N;
% x = h*[-N1:N1-1]';
% 
% % frequency:[-pi*N/(2L),pi*N/(2L)]
% ksi1 = pi/L*[0:N1-1 -N1:-1]'; % 这是为了避免做fft之后要用fftshift
% %% caculate intergral
% % J = Compute_Integral_term(ksi1,N);
% load('N3000_1D_integral.mat'); 
% 
% %% define differential matrix
% e = ones(N,1);
% D = spdiags([-e,e],[0,1],N,N);
% D(N,N) = 0;
% %% smooth source function
% F = @(x) 5*exp(-(x+5).^2)+3*exp(-2*x.^2)+exp(-(x-4).^2);
% % Flinear_hat = @(x) 5.0/sqrt(2)*exp(-x.^2/4+5i*x)...
% %     +3.0/2.0*exp(-x.^2/8)+1.0/sqrt(2)*exp(-x.^2/4-4i*x);
% 
% f = feval(F,x);
% delta_vec = [0.5*1e-2,1e-2,5*1e-2,0.1];
% delta_num = size(delta_vec,2);
% lambda_opt = [35000,	25000,	2500,	15000];
% % TVBFD start
% beta=5; 
% broadcast = J;
% parfor k = 1:delta_num
% TVBFD_figure1D(f,broadcast,x,ksi1,delta_vec(k),lambda_opt(k),beta,D,0);
% end
% %% piecewise_linear_source function
% f_piecewise_linear = 2*piecewise_linear_source(x,-4,1)+...
%     piecewise_linear_source(x,0,2)+3*piecewise_linear_source(x,4,1);
% piecewise_delta_vec = [0.5*1e-3,1e-3,0.5*1e-2,1e-2,5*1e-2,0.1];
% piecewise_delta_num = size(piecewise_delta_vec,2);
% piecewise_lambda_opt =  [2.74102e+06,	685608,	52500,	52500,	25000	,15000];
% parfor k = 1:piecewise_delta_num
%     TVBFD_figure1D(f_piecewise_linear,broadcast,x,ksi1,...
%         piecewise_delta_vec(k),piecewise_lambda_opt(k),beta,D,1);
% end
% %%  2D problem
% clear;
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
% % 将Ksinorm拉直为列向量,方便计算积分值
% Ksi = Ksinorm(:);
% J = Compute_Integral_term(Ksi,N);
%% define differential matrix
% differial matrix such that D2*u = grad_x u, u*D1' = grad_y u;
e = ones(N,1);
D1 = spdiags([-e,e],[0,1],N,N);
D1(N,N) = 0;

% D2 = spdiags([-e,e],[0,1],N,N); % 默认输入是方阵,所以可注释掉
% D2(N,N) = 0;
%% continuous source function
F = @(x,y) 2*exp(-(x+5).^2-y.^2)+5*exp(-(x-4).^2-y.^2);
% F_hat = @(x,y) 2.0/sqrt(2)*exp(-x.^2/4+5i*x).*exp(-y.^2/4)...
%     +5.0/sqrt(2)*exp(-x.^2/4-4i*x).*exp(-y.^2/4);
f = feval(F,X,Y);

delta_vec = [0.5*1e-2,1e-2,5*1e-2,0.1,0.2];
delta_num = size(delta_vec,2);
lambda_opt = [15000,	15000	,15000	,5000,	5000];
%% use TVBFD
beta = 20;  % 二维β用20效果比较好
broadcast = J;
parfor k = 1:delta_num
    TVBFD_figure2D(f,broadcast,X,Y,Ksi1,Ksi2,Ksinorm2,delta_vec(k),...
        lambda_opt(k),beta,D1,0);
end
  %% challenging example
X1 = X(:);Y1 = Y(:);
temp = piecewise_constant_source(X1,Y1);
f_piecewise_constant = reshape(temp,N,N);

piecewise_delta_vec = [0.5*1e-3,1e-3,0.5*1e-2,1e-2,5*1e-2,0.1,0.2];
piecewise_delta_num = size(piecewise_delta_vec,2);
piecewise_lambda_opt =  [2.74102e+06,	685608	,27537.6,	15000	,15000	,15000,	15000];
parfor k = 1:piecewise_delta_num
    TVBFD_figure2D(f_piecewise_constant,broadcast,X,Y,Ksi1,Ksi2,...
        Ksinorm2,piecewise_delta_vec(k),piecewise_lambda_opt(k),beta,D1,1);
end
end