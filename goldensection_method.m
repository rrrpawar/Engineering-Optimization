% To optimization using Golden Section Search Method
function [L,x1,x2]=goldensection_method(a,b,e,f)
syms x;
%a=-5; %a= input('Enter the lower bound ');
%b=20; %b= input('Enter the upper bound '); 
GR=0.618;
L= abs(b-a);
x1=a+(1-GR)*(b-a);       
x2=a+GR*(b-a);            
f(x)=10+(x-4).^2 +10*exp(-x/4-2)-3*sin(0.2*x*pi);
i=0;                         % iteration count
while (L>0.00001)            % termination criterion
    i=i+1;
    L= abs(b-a);
    fx1=subs(f(x),x1);
    fx2=subs(f(x),x2);
    if fx1<fx2
        b=x2;                % eliminating area (x2-b)
        x2=x1;
        x1=a+(1-GR)*(b-a);
    else 
        a=x1;                % eliminating area (a-x1)
        x1=x2;
        x2=a+GR*(b-a);
    end
    double(x1);

t=linspace(x1,10,100); 
plot(t,f(t),'k','LineWidth',2);
hold on
plot(x1,fx1,'ro');
hold on
plot(x2,fx2,'bo');
title('Golden Section method'); 
xlabel('x'); 
ylabel('f(x)');
legend('function','x1','x2');
grid on 
end
 disp(['The optimum value of x is ' ,num2str(x1)]);
 disp(['The number of iteration is ',num2str(i)]);