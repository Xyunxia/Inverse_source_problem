function TVBFD_figure1D(f,J,x,ksi1,delta,lambda_opt,beta,D,flag)
g_delta_hat = generate_noise_measure(f,delta,J);
fm = oneD_TVBFD_apply(lambda_opt,beta,g_delta_hat,J,ksi1,delta,D);
if flag == 1
    source_type = 'Piecewise linear';
else
    source_type = 'smooth';
end
b = norm(f-fm)/norm(f);
figure(3)
plot(x,fm,'b:',x,f,'r--','linewidth',1.5);
xlabel('x');
ylabel('f(x)');
legend('f^m','f exact');
finalname=[source_type,',','\delta =',num2str(delta),',', 'relative err = ',num2str(b)];
title(finalname);
finame = [source_type,',','delta =',num2str(delta),',', 'relative err =',num2str(b),'.png'];
print(3,'-dpng',finame)   % 保存为png格式的图片
end