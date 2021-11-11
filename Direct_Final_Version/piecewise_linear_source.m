function f = piecewise_linear_source(x,a,b)
temp1 = double( (x>=a-b) & (x<=a) );
temp2 = double( (x>a) & (x<=a+b) );
f = (x-a+b).*temp1 + (a+b - x).*temp2;
f = f/b;
end
