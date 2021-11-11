function [J] = Compute_Integral_term(Ksi,N)
% denote:
% $w(s) = E_{\alpha, \alpha}\left(-|\xi|^{\gamma}(T-s)^{\alpha}\right) h(s) $
% $$
% w(s) \approx \frac{1}{\Delta t}\left[w\left(t_{i}\right)\left(t_{i+1}-s\right)+w\left(t_{i+1}\right)\left(s-t_{i}\right)\right], \quad t_{i} \leq s \leq t_{i+1}\label{gs12}
% $$
% the code is used as caculate the intergral:
% $$
% \begin{aligned}
% &\int_{0}^{T}(T-s)^{\alpha-1} w(s) d s=\sum_{i=0}^{k-1} \int_{t_{i}}^{t_{i+1}}(T-s)^{\alpha-1} w(s) d s \\
% &\approx \frac{1}{\Delta t} \sum_{i=0}^{k-1}\left[w\left(t_{i}\right) \int_{t_{i}}^{t_{i+1}}\left(t_{i+1}-s\right)(T-s)^{\alpha-1} d s+w\left(t_{i+1}\right) \int_{t_{i}}^{t_{i+1}}\left(s-t_{i}\right)(T-s)^{\alpha-1} d s\right]
% \end{aligned}
% $$

%   Input: 
% Ksi:  must be a vector.
% for two Dimensional case, suppose the discrete mesh grid in frequency is
% [Ksi1,Ksi2] = meshgrid(ksi1,ksi1);
% Ksinorm2 = Ksi1.^2 + Ksi2.^2;
% Ksinorm = sqrt(Ksinorm2);
% then Ksi =  Ksinorm(:);
% N: Spatial discrete scale size
gamma=1;
alpha=0.5;

M = 128;
dt = 1/M;
t = [0:dt:1];
vecLen = size(Ksi,1);

J = zeros(vecLen,1);
%% 1D problem
 % w 每一列是时间在变，每一行是ksi在变
KT = -(abs(Ksi).^gamma)*((1-t).^alpha);
E = mlf(alpha,alpha,KT,6);
w = E.*exp(-t); 

parfor i = 1:M
    a = t(i);
    b = t(i+1);
    fun1 = @(s)((b-s).*(1-s).^(alpha-1));
    fun2 = @(s)((s-a).*(1-s).^(alpha-1));
    J = J + w(:,i)*quadgk(fun1,a,b)...
       + w(:,i+1)*quadgk(fun2,a,b);
end
J = M*J;
if N == vecLen
    save('N3000_1D_Integral','J');
else
    J = reshape(J,N,N);
    save('N500_2D_Integral.mat','J');
end
% w = zeros(vecLen,size(t,2));
% for k = 1:vecLen
% w(k,:) = (mlf(alpha,alpha, -(abs(Ksi(k))^gamma)*(1-t).^alpha ,6).*exp(-t)); 
% end
% for  k = 1:vecLen
%     for i = 1:M
%         a = t(i);
%         b = t(i+1);
%         fun1 = @(s)((b-s).*(1-s).^(alpha-1));
%         fun2 = @(s)((s-a).*(1-s).^(alpha-1));
%        J(k) = J(k) + w(k,i)*quadgk(fun1,a,b)...
%            + w(k,i+1)*quadgk(fun2,a,b);
%     end
% end
% J = M*J;
% J = reshape(J,N,N);
end