% To optimize the function using Bisection method
function [xm]=bisection_method(a,b,e,f)
%a=-5
%b=20
syms x
xm=(a+b)/2;
f(x)=10+(x -4).^2+10*exp(-x/4-2)-3*sin(0.2*x*pi);
f1(x)=diff(f);
i=0;
while abs(f1(xm))>0.00001       % convergence criterion 
      i=i+1;                    % number of iterations 
      xm=(a+b)/2;
   if f1(xm)>0 
      b=xm;                     % dropping [b-xm]
  else 
      a=xm;                     % dropping [a-xm] 
   end
  double(xm)
  
t=linspace(a,20,100); 
plot(t,f(t),'k','LineWidth',2);
plot(t,f1(t),'k','LineWidth',2);
title('Bisection method'); 
xlabel('x'); 
ylabel('f(x)');
grid on 
hold on 
plot(xm,f(xm),'rx');            % function plot
hold on
plot(xm,f1(xm),'bx');           % derivative plot
xlabel('x');
ylabel('f(x)');
legend('function','function values','derivative');
grid on
end
disp(['The optimum value of x is ',num2str(xm)]);
disp(['Number of iterations = ',num2str(i)]);