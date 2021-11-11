function test_direct_method
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
% %% smooth source function
% F = @(x) 5*exp(-(x+5).^2)+3*exp(-2*x.^2)+exp(-(x-4).^2);
% % Flinear_hat = @(x) 5.0/sqrt(2)*exp(-x.^2/4+5i*x)...
% %     +3.0/2.0*exp(-x.^2/8)+1.0/sqrt(2)*exp(-x.^2/4-4i*x);
% 
% f = feval(F,x);
% 
% delta_vec = [0.5*1e-2,1e-2,5*1e-2,0.1];
% delta_num = size(delta_vec,2);
% %% direct solve method
% % use a priori parameter choice rule
% priori_relative_err = zeros(1,delta_num);
% % use a posteriori parameter choice rule
% poster_relative_err = zeros(1,delta_num);
% broadcast = J;
% parfor k = 1:delta_num
%     [priori_relative_err(k),poster_relative_err(k)] = direct_solve_method(delta_vec(k),f,broadcast,...
%         ksi1,N1,x,0);
% end
% fid = fopen('result.txt','a+'); % 读写方式打开文件
% fprintf(fid,'\n%s\n','smooth source function:');
% fprintf(fid,'%s\t','delta:');
% fprintf(fid,'%g\t%g\t%g\t%g\n',delta_vec);
% fprintf(fid,'%s\t','priori relative err:');
% fprintf(fid,'%g\t%g\t%g\t%g\n',priori_relative_err);
% fprintf(fid,'%s\t','poster relative err:');
% fprintf(fid,'%g\t%g\t%g\t%g\n',poster_relative_err);
% 
% %% piecewise_linear_source function
% f_piecewise_linear = 2*piecewise_linear_source(x,-4,1)+...
%     piecewise_linear_source(x,0,2)+3*piecewise_linear_source(x,4,1);
% 
% piecewise_delta_vec = [0.5*1e-3,1e-3,0.5*1e-2,1e-2,5*1e-2,0.1];
% piecewise_delta_num = size(piecewise_delta_vec,2);
% 
% %% direct solve method
% % use a priori parameter choice rule
% piecewise_priori_relative_err = zeros(1,piecewise_delta_num);
% % use a posteriori parameter choice rule
% piecewise_poster_relative_err = zeros(1,piecewise_delta_num);
% parfor k = 1:piecewise_delta_num
%     [piecewise_priori_relative_err(k),piecewise_poster_relative_err(k)]...
%         = direct_solve_method(piecewise_delta_vec(k),f_piecewise_linear,...
%        broadcast,ksi1,N1,x,1);
% end
% fid = fopen('result.txt','a'); % 写方式打开文件,追加数据到文件末尾
% fprintf(fid,'\n %s \n','piecewise linear source function:');
% fprintf(fid,'%s\t','delta:');
% fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_delta_vec);
% fprintf(fid,'%s\t','priori relative err:');
% fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_priori_relative_err);
% fprintf(fid,'%s\t','poster relative err:');
% fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_poster_relative_err);
% clear;
%%  2D problem
% E = 120;  
% p=3;
% gamma=1;
% delta_vec = [0.5*1e-3,1e-3,0.5*1e-2,1e-2,5*1e-2,0.1,0.2];
% v = (E./delta_vec).^(1/(p+gamma));
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

%% caculate intergral
load('N500_2D_integral');
% Ksinorm2 = Ksi1.^2 + Ksi2.^2;
% Ksinorm = sqrt(Ksinorm2);
% % 将Ksinorm拉直为列向量，方便计算积分值
% J = Compute_Integral_term(Ksinorm(:),N);
%% continuous source function
F = @(x,y) 2*exp(-(x+5).^2-y.^2)+5*exp(-(x-4).^2-y.^2);
% F_hat = @(x,y) 2.0/sqrt(2)*exp(-x.^2/4+5i*x).*exp(-y.^2/4)...
%     +5.0/sqrt(2)*exp(-x.^2/4-4i*x).*exp(-y.^2/4);
f = feval(F,X,Y);

X1 = X(:);Y1 = Y(:);
Ksi = [Ksi1(:) Ksi2(:)];
S = [X1 Y1];

delta_vec = [0.5*1e-2,1e-2,5*1e-2,0.1,0.2];
delta_num = size(delta_vec,2);

%% direct solve method
% use a priori parameter choice rule
priori_relative_err = zeros(1,delta_num);
% use a posteriori parameter choice rule
poster_relative_err = zeros(1,delta_num);

broadcast = J;
parfor k = 1:delta_num
    [priori_relative_err(k),poster_relative_err(k)] = ...
        direct_solve_method(delta_vec(k),f,broadcast,...
        Ksi,N1,S,2);
end
fid = fopen('result.txt','a'); % 读写方式打开文件
fprintf(fid,'\n %s \n',' continuous source function:');
fprintf(fid,'%s\t','delta:');
fprintf(fid,'%g\t%g\t%g\t%g\t%g\n',delta_vec);
fprintf(fid,'%s\t','priori relative err:');
fprintf(fid,'%g\t%g\t%g\t%g\t%g\n',priori_relative_err);
fprintf(fid,'%s\t','poster relative err:');
fprintf(fid,'%g\t%g\t%g\t%g\t%g\n',poster_relative_err);

%% challenging example
temp = piecewise_constant_source(X1,Y1);
f_piecewise_constant = reshape(temp,N,N);

piecewise_delta_vec = [0.5*1e-3,1e-3,0.5*1e-2,1e-2,5*1e-2,0.1,0.2];
piecewise_delta_num = size(piecewise_delta_vec,2);

%% direct solve method
% use a priori parameter choice rule
piecewise_priori_relative_err = zeros(1,piecewise_delta_num);
% use a posteriori parameter choice rule
piecewise_poster_relative_err = zeros(1,piecewise_delta_num);
parfor k = 1:piecewise_delta_num
    [piecewise_priori_relative_err(k),piecewise_poster_relative_err(k)] ...
        = direct_solve_method(piecewise_delta_vec(k),f_piecewise_constant,...
        broadcast,Ksi,N1,S,3);
end
fid = fopen('result.txt','a'); % 写方式打开文件,追加数据到文件末尾
fprintf(fid,'\n %s\n','piecewise constant source function:');
fprintf(fid,'%s\t','delta:');
fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_delta_vec);
fprintf(fid,'%s\t','priori relative err:');
fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_priori_relative_err);
fprintf(fid,'%s\t','poster relative err:');
fprintf(fid,'%g\t%g\t%g\t%g\t%g\t%g\t%g\n',piecewise_poster_relative_err);
end