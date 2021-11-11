function  [d_tidle_x,d_tidle_y,dx,dy]=solveDsubproblem(d_tidle_x_old,d_tidle_y_old,dx_old,dy_old,u,gama,D1)
% declaration: a vectorial shrinkage is defined as:
% shrink(x,gama) = x/(abs(x))*max(abs(x) - gama,0);
% the function has calculate the result of d- subproblem
% in the isotropic TV Denoising
% dx^(k+1) = max(s^k -1/gama,0)(dx^(k)+\partial_x u-d_tidle_x)/s^k and
% dy^(k+1) = max(s^k -1/gama,0)(dy^(k)+\partial_y u-d_tidle_y)/s^k
% s^k = sqrt(|dx^(k)+\partial_x u-d_tidle_x|^2+|dy^(k)+\partial_y u-d_tidle_y|^2)
% simultaneous update d_tidle_x and d_tidle_y

% differial matrix such that D2*u = grad_x u, u*D1' = grad_y u;
%    e = ones(max(width,height),1);
%    D1 = spdiags([-e,e],[0,1],height,height);
%    D1(height,height) = 0;
%    D2 = spdiags([-e,e],[0,1],width,width);
%    D2(width,width) = 0;
   grad_y = u*D1';
   grad_x = D1*u;
   thresh = 1.0/gama;

   dx_new = dx_old + grad_x - d_tidle_x_old ;
   dy_new = dy_old + grad_y - d_tidle_y_old ;
   grad_u_norm = sqrt(dx_new.^2+dy_new.^2);
   temp = max(grad_u_norm-thresh,0);
   dx = temp.*dx_new./grad_u_norm;
   dy = temp.*dy_new./grad_u_norm;
   d_tidle_x = 2*dx - dx_new;
   d_tidle_y = 2*dy - dy_new;
   clear temp;
   clear grad_u_norm;
   
   dx(isnan(dx)) = 0;
   dy(isnan(dy)) = 0;
   d_tidle_x(isnan(d_tidle_x)) = 0;
   d_tidle_y(isnan(d_tidle_y)) = 0;
   dx(isinf(dx)) = 0;
   dy(isinf(dy)) = 0;
   d_tidle_x(isinf(d_tidle_x)) = 0;
   d_tidle_y(isinf(d_tidle_y)) = 0;
 
end