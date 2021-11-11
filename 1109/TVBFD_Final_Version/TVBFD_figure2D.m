function TVBFD_figure2D(f,J,X,Y,Ksi1,Ksi2,Ksinorm2,delta,lambda,beta,D1,flag)
%% generate the noise mesurement of the u_hat(ksi,T)
g_delta_hat = generate_noise_measure(f,delta,J);
fm = TVBFD(g_delta_hat,J,Ksi1,Ksi2,Ksinorm2,delta,...
        lambda,beta,D1);
if flag == 1
    source_type = 'Piecewise constant';
else
    source_type = 'continuous';
end
b = norm(f-fm)/norm(f);
figure(5)
mesh(X,Y,fm);
xlabel('x');
ylabel('y');
zlabel('f_{TV}^m(x,y)');
title(['delta = ',num2str(delta),'relative err = ',num2str(b),'.png'])
finame = ([source_type,',','delta =',num2str(delta),',', 'relative err =',num2str(b),'.png']);
print(5,'-dpng',finame)
end