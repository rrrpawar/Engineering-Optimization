% Powell's method to optimize 
function[x1,x2,dx]=powell_method(a,b,e,f)
clear all
close all
x1=-5; 
dx=2;
syms x
f(x)=10+(x -4).^2+10*exp(-x/4-2)-3*sin(0.2*x*pi);
x2=x1+dx;
    if f(x1)>f(x2)                  % to find x3
        x3=x1+2*dx;
       else 
        x3=x1-dx;
    end
    
    if f(x1)<f(x2) && f(x1)<f(x3)   % finding xmin by finding the minimum value function
        xmin=x1;
    elseif f(x2)<f(x3) && f(x2)<f(x1)
        xmin=x2;
    elseif f(x3)<f(x1) && f(x3)<f(x2)
        xmin=x3;
    end 
   fmin= f(xmin);
   
   if fmin == f(x1)                 % to find xmin from fmin
     xmin = x1;
     
 elseif fmin == f(x2)
     xmin = x2;
     
 elseif fmin == f(x3)
     xmin = x3;
   end 
   fnew=fmin;
   
a2=(f(x2)-f(x1))/(x2-x1);
a3=(((f(x3)-f(x1))/(x3-x1))-((f(x2)-f(x1))/(x2-x1)))/(x3-x2);
xn=(((x2+x1)/2)-a2/(2*a3));          % finding xnew
fnew=f(xn); 
   
i=0; j=1; l=1;

while j>=0.00001 || l>=0.00001       % convergence criteria
    
    i=i+1;
     
     if f(x1)<f(x2) && f(x1)<f(x3)   % finding xmin by finding the minimum value function
        xmin=x1;
 elseif f(x2)<f(x3) && f(x2)<f(x1)
        xmin=x2;
 elseif f(x3)<f(x1) && f(x3)<f(x2)
        xmin=x3;
    end 
   fmin= f(xmin);
   
     if fmin == f(x1)
        xmin = x1;  
 elseif fmin == f(x2)
        xmin = x2;
 elseif fmin == f(x3)
        xmin = x3;
    end 
   
a1=f(x1);                         % a1,a2,a3 for the quadratic curve
a2=(f(x2)-f(x1))/(x2-x1);
a3=(((f(x3)-f(x1))/(x3-x1))-((f(x2)-f(x1))/(x2-x1)))/(x3-x2);
xn=(((x2+x1)/2)-a2/(2*a3));       % finding xnew
fnew=f(xn); 
j=abs(fmin-fnew);
l=abs(xmin-xn);

z=[x1 xmin xn];                   % to sort and find minimum
z1=sort(z);
x1=z1(1,1);
x2=z1(1,2);
x3=z1(1,3);         

double(xn)
t=linspace(x1,10,100); 
plot(t,f(t),'k','LineWidth',2);
title('Powell estimation method'); 
xlabel('x'); 
ylabel('f(x)');
grid on 
hold on 
plot(x1,f(x1), '.r', 'MarkerSize', 20);
hold on
plot(x2,f(x2),'*k','MarkerSize',10);
hold on
plot(x3,f(x3),'ob','MarkerSize',10);
legend('function','x1','x2','x3');
hold on
end
disp(['The number of iterations: ',num2str(i)]);
fprintf('The value of x is %d',double(xn));