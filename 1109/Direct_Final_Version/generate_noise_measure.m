function [g_delta_hat,g,g_delta] = generate_noise_measure(f,delta,J)
[r,w] = size(f);
if w==1
    f_hat = fft(f);
    % solve the time-space fractional partial eqation
    g_hat = f_hat.*J;
    g = abs(ifft(g_hat));
    %% generate the noise mesurement of the u_hat(ksi,T)
    g_delta = g + delta*randn(r,w);
    g_delta_hat = fft(g_delta);
else
    f_hat = fft2(f);
    % solve the time-space fractional partial eqation
    g_hat = f_hat.*J;
    g = abs(ifft2(g_hat));
    %% generate the noise mesurement of the u_hat(ksi,T)
    g_delta = g + delta*randn(r,w);
    g_delta_hat = fft2(g_delta);
end