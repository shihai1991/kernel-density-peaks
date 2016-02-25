function [s] = kernal( i,j,a,b)

%% 估计想表达径向基核函数(RBF)
% delta = 16;
% s =exp(-((i-a)*(i-a)+(j-b)*(j-b))/(2*delta^2)); % 效果好

%% 傅里叶核函数
% delta = 5;
% s = (1-delta^2)/(2*(1-2*delta*cos(sqrt(((i-a)*(i-a)+(j-b)*(j-b))))+delta^2));
 
%% 二次有理核
%  delta = 24;
%  s = 1-(i-a)*(i-a)+(j-b)*(j-b)/((i-a)*(i-a)+(j-b)*(j-b)+delta);
 
 %% 多元二次核
%  delta = 2;
%  s = sqrt((i-a)*(i-a)+(j-b)*(j-b)+delta^2);
 
 %% 逆多元二次核
 delta = 2.5;
 s = 1/sqrt((i-a)*(i-a)+(j-b)*(j-b)+delta^2); %效果好
end