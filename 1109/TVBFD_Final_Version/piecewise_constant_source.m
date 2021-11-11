function f = piecewise_constant_source(x,y)
% 用&而不用&&是因为 ：
% &既能判断标量值又能判断向量值，当输入是数组时返回值也是数组。
% 而 && 只能处理标量值输入，当输入数组返回的也是数值。
temp1 = double((x<-1&x>-3)&(y>1&y<5));  
temp2 = double(x>2&x<6)&(y>-6&y<-2);
f = 2.*temp1 + 4.*temp2;
end