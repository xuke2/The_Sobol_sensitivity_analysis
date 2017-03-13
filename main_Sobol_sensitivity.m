% ---------------------------------------------------------------------
% Using Monte Carlo method to analysis the global sensitivities of
% each variables
% Sampling method : sobol sequence
% Test function         : gmath function
% ---------------------------------------------------------------------
% Reference:
% 1.Archer, G. E. B., A. Saltelli, and Sobol I. M.. "Sensitivity measures,
%   ANOVA-like techniques and the use of bootstrap." Journal of Statistical
%   Computation and Simulation 58(2),1997: 99-120 .
% 2.Sobol, I. M. "Global sensitivity indices for nonlinear mathematical
%   models and their Monte Carlo estimates." Mathematics and computers in
%   simulation 55(1), 2001: 271-280.
% 3.I.M. Sobol¡¯, S. Tarantola, D. Gatelli, S.S. Kucherenko, W. Mauntz,
%   Estimating the approximation error when fixing unessential factors in
%   global sensitivity analysis, Reliability Engineering & System Safety,
%   92(7) 2007: 957-960.
% ---------------------------------------------------------------------
clear all;
clc;

% ---------------------------------------------------------------------
% settings of gmath function
% the dimension of gmath function     
n=8;                
% the patemeters of gmath function which  controls the sensitivity of each variable   
a=[0,1,4.5,9,99,99,1,99];   
% the number of sample points for  Monte Carlo simulation
num_sample=1e4;   

% ---------------------------------------------------------------------
% initilize some parameters
Di=zeros(n,1);
Dii=zeros(n,1);
S=zeros(n,1);
TSI=zeros(n,1);

% ---------------------------------------------------------------------
% generate sample points using Sobol sequence
q=qrandstream('sobol',2*n);
sample=qrand(q,num_sample);


% ---------------------------------------------------------------------
% calculate F0 and D
f=zeros(1,num_sample);
ff=zeros(1,num_sample);
for i=1:num_sample
    f(i)=gmath(sample(i,1:n),a);
    ff(i)=gmath(sample(i,n+1:2*n),a);
end
F0=(sum(f.*ff)/num_sample)^0.5;
D=sum(f.^2)/num_sample-F0^2;

% ---------------------------------------------------------------------
% calculate Di and Dii-
for k=1:n
    for i=1:num_sample
        kesi1=sample(i,1:n);
        kesi2=sample(i,n+1:2*n);
       
        kesi3=kesi2;
        kesi3(:,k)=kesi1(:,k);
        
        kesi4=kesi1;
        kesi4(:,k)=kesi2(:,k);
        
        Di(k,:)=Di(k,:)+f(i)*(gmath(kesi3,a)-ff(i));
        Dii(k,:)=Dii(k,:)+f(i)*(gmath(kesi4,a)-ff(i));
    end
    
    Di(k,:)=Di(k,:)/num_sample;
    Dii(k,:)=Dii(k,:)/num_sample;
    S(k)=Di(k,:)/D;
    TSI(k)=Dii(k,:)/D;
end

% ---------------------------------------------------------------------
% print the results on screen
disp('1th sensitivity  global sensitivity')
disp([S   TSI]);
bar([S   TSI]);





