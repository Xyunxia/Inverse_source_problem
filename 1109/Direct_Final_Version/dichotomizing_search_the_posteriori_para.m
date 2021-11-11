function [Sv_opt_hat,v_opt] = dichotomizing_search_the_posteriori_para(g_delta_hat,delta,Ksi,N1,flag,tau) 
% purpore: use the dichotomizing search the posteriori parameter v
% such that G(v_{\delta}) = ||S_{v_{\delta}}-g^{\delta}||_{L_2} = \tau \delta
wequal = tau*delta;
supe = wequal+5*1e-3;
infe = wequal-5*1e-3;
m=1; 
[r,w] = size(g_delta_hat);

Cf = sqrt(r*w);
% ksi1 = pi/L*[0:N1-1 -N1:-1]';  % 这是为了避免做fft之后要用fftshift
if flag ==0 || flag == 1 
    %% 1D
    a = Ksi(2)+ 0.01053;  b = max(Ksi)+1.0234;
    temp1 = find(abs(Ksi)<a);
    Sva_hat = zeros(r,w);
    Sva_hat(temp1) = g_delta_hat(temp1);
    temp2 = find(abs(Ksi)<b);
    Svb_hat = zeros(r,w);
    Svb_hat(temp2) = g_delta_hat(temp2);
    difa = norm(Sva_hat-g_delta_hat)/Cf;
    difb = norm(Svb_hat-g_delta_hat)/Cf;
 
    while(1)
       mid = (a+b)/2.0;
       xvalid = find(abs(Ksi)<mid);
        Sv_hat = zeros(r,w);
        Sv_hat(xvalid) = g_delta_hat(xvalid);
        dif = norm(Sv_hat-g_delta_hat)/Cf;
        if difa>infe && difa<supe
            v_opt = a;
            Sv_opt_hat = Sva_hat;
            break;
        elseif difb>infe && difb<supe
            v_opt = b;
            Sv_opt_hat = Svb_hat;break;
        elseif dif > infe && dif<supe
            v_opt = mid;
            Sv_opt_hat = Sv_hat;break;
        end
        if dif>supe 
            if difa>supe && difb<infe
                a = mid;
                difa = dif;
            elseif difa>supe && difb>supe
                b = b+1/N1;
                temp2 = find(abs(Ksi)<b);
                Svb_hat = zeros(r,w);
                Svb_hat(temp2) = g_delta_hat(temp2);
                difb = norm(Svb_hat-g_delta_hat)/Cf;
            end
        else
             if difa>supe && difb<infe
                 b = mid;
                 difb = dif;
             elseif difa<infe && difb<infe
                 a = a-1/N1;
                 temp1 = find(abs(Ksi)<a);
                Sva_hat = zeros(r,w);
                Sva_hat(temp1) = g_delta_hat(temp1);
                difa = norm(Sva_hat-g_delta_hat)/Cf;
             end
        end
        if m>5*1e+3 
            fprintf('surpass the max iteration!\n');
            v_opt = 0;
            Sv_opt_hat = 0;
            break;
        end
        if abs(b-a)<5*1e-4
            fprintf('the cell is too small!\n');
            v_opt = mid;
            Sv_opt_hat = Sv_hat;
            break;
        end
        m = m+1;
    end
else
    %% 2D
%     [Ksi,Ksi2] = meshgrid(Ksi,Ksi);
    temp1 = Ksi(:,1);temp2 = Ksi(:,2);
    a = 1.010534;  b = max(temp1)+1.02346;
    
    ind1 = find(abs(temp1)<a);
    ind2 = find(abs(temp2)<a);
    index = intersect(ind1,ind2'); % 找两个数组共同包含的元素
    vector_g_delta_hat = g_delta_hat(:);
    Valid_temp = zeros(r*w,1);
    Valid_temp(index) = vector_g_delta_hat(index);
    Sva_hat = reshape(Valid_temp,r,w);
    
    ind1 = find(abs(temp1)<b);
    ind2 = find(abs(temp2)<b);
    index = intersect(ind1,ind2'); % 找两个数组共同包含的元素
    Valid_temp = zeros(r*w,1);
    Valid_temp(index) = vector_g_delta_hat(index);
    Svb_hat = reshape(Valid_temp,r,w);
    
    difa = norm(Sva_hat-g_delta_hat)/Cf;
    difb = norm(Svb_hat-g_delta_hat)/Cf;
    while(1)
        mid = (a+b)/2.0;
        ind1 = find(abs(temp1)<mid);
        ind2 = find(abs(temp2)<mid);
        index = intersect(ind1,ind2'); % 找两个数组共同包含的元素
        Valid_temp = zeros(r*w,1);
        Valid_temp(index) = vector_g_delta_hat(index);
        Sv_hat = reshape(Valid_temp,r,w);
        dif = norm(Sv_hat-g_delta_hat)/Cf;
        if difa>infe && difa<supe
            v_opt = a;
            Sv_opt_hat = Sva_hat;
            break;
        elseif difb>infe && difb<supe
            v_opt = b;
            Sv_opt_hat = Svb_hat;break;
        elseif dif > infe && dif<supe
            v_opt = mid;
            Sv_opt_hat = Sv_hat;break;
        end
        if dif>supe 
            if difa>supe && difb<infe
                a = mid;
                difa = dif;
            elseif difa>supe && difb>supe
                b = b+1/N1;
                ind1 = find(temp1<b);
                ind2 = find(temp2<b);
                index = intersect(ind1,ind2'); % 找两个数组共同包含的元素
                Valid_temp = zeros(r*w,1);
                Valid_temp(index) = vector_g_delta_hat(index);
                Svb_hat = reshape(Valid_temp,r,w);
                difb = norm(Svb_hat-g_delta_hat)/Cf;
            end
        else
             if difa>supe && difb<infe
                 b = mid;
                 difb = dif;
             elseif difa<infe && difb<infe
                a = a-1/N1;
                ind1 = find(abs(temp1)<a);
                ind2 = find(abs(temp2)<a);
                index = intersect(ind1,ind2'); % 找两个数组共同包含的元素
                Valid_temp = zeros(r*w,1);
                Valid_temp(index) = vector_g_delta_hat(index);
                Sva_hat = reshape(Valid_temp,r,w);
                difa = norm(Sva_hat-g_delta_hat)/Cf;
             end
        end
        if m>5*1e+3 
            fprintf('surpass the max iteration!\n');
            v_opt = 0;
            Sv_opt_hat = 0;
            break;
        end
        if abs(b-a)<5*1e-4
            fprintf('the cell is too small!\n');
            v_opt = mid;
            Sv_opt_hat = Sv_hat;
            break;
        end
        m = m+1;
    end

end