function [ y ] = gmath( x,a )
% test function
if size(x,2)~=size(a,2)
    error('x and a must have the same size');
end
m=size(x,2);

y=1;

for i=1:m
    y=y*(abs(4*x(i)-2)+a(i))/(1+a(i));
end




end

