% To optimize the function using interval halving method
function[x1,x2,xm,L]=intervalhalving_method(a,b,e,f)
clear all
close all
syms x
a=-5; %a=input('Enter the lower bound ');  
b=20; %b=input('Enter the upper bound '); 
L=abs(b-a);
x1=a+L/4;
x2=b-L/4;
xm=(a+b)/2;                       % The midpoint of a-b
%f(x)=2*x^2+(16/x);
f(x)=10+(x -4).^2+10*exp(-x/4-2)-3*sin(0.2*x*pi);
fx1=f(x1);
fx2=f(x2);
fxm=f(xm);
i=0;
while L>=0.00001                  % Termination criterion 
      i=i+1;                      % counter to calculate number of iterations
      L=abs(b-a);
      x1=a+L/4;
      x2=b-L/4;
      fx1=f(x1);
      fx2=f(x2);
      fxm=f(xm);
      if fx1<fxm                  % dropping [xm-b]
         b=xm;
         xm=x1;
         fxm=fx1;
  elseif fx2<fxm                  % dropping [a-xm]
         a=xm;
         xm=x2;
         fxm=fx2;
    else fx2>=fxm;                % dropping [a-x1] and [x2-b]
         a=x1;
         b=x2; 
         double(xm);
      end
     
t=linspace(a,10,100); 
plot(t,f(t),'k','LineWidth',2);
title('Interval halving method'); 
xlabel('x'); 
ylabel('f(x)');
grid on 
hold on -
    plot(x1,f(x1), '.r', 'MarkerSize', 20);
    hold on
    plot(x2,f(x2),'*b','MarkerSize',10);
    hold on
    legend('function','x1','x2');
    hold on
      
end
disp(['The optimum value of x = ',num2str(xm)]);
disp(['Number of iterations = ',num2str(i)]);
