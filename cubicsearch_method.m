% To optimize the function using cubic search method
function[x1]= cubicsearch_method(a,b,e,f)
x0=-5; %input('Enter the initial value: ');
dx=2;  %input('Enter the step size: ');
syms x
f(x)= 10+(x -4).^2+10*exp(-x/4-2)-3*sin(0.2*x*pi);
x1=x0+dx;
f1(x)=diff(f,x);
f1x0=f1(x0);
f1x1=f1(x1);
i=0;
k=1;
while (f1x0*f1x1)>0        % to calculate x2
    if  f1x0<0
        x2=x1+(2^k)*dx;
    else
        x2=x1-(2^k)*dx;
    end
    k=k+1;
    x0=x1;
    x1=x2;
    f1x0=f1(x0);
    f1x1=f1(x1);
end

fx0=f(x0);
fx1=f(x1);
f1x0=f1(x0);
f1x1=f1(x1);

z=3*((fx0-fx1)/(x1-x0))+f1x1+f1x0;

    if x0<x1
       w=sqrt(z^2-(f1x0*f1x1)); % to find the value of w and mu
    else
       w=-sqrt(z^2-(f1x0*f1x1));
    end

mu=(f1x1+w-z)/(f1x1-f1x0+2*w);  % to calculate mu

    if mu<0                     
       xn=x1;
elseif mu>= 0 && mu<=1
       xn=x1-mu*(x1-x0);
    else
       xn=x0;
    end
fxn=f(xn);
f1xn=f1(xn);
while abs(f1xn)>=0.00001 && abs(((xn-x0)/x0))>=0.00001 % convergence criterion
        while fxn>fx0
            xn=xn+(xn-x0)/2;
        end
    if (f1xn*f1x0)<0
        x1=x0;
        x0=xn;
        fx0=f(x0);
        fx1=f(x1);
        f1x0=f1(x0);
        f1x1=f1(x1);
    elseif (f1xn*f1x1)<0
        x0=xn;
    end

z=3*((fx0-fx1)/(x1-x0))+f1x0+f1x1; % to calculate z

    if x0<x1
       w=sqrt(z^2-(f1x0*f1x1));   % to calculate w
    else
       w=-sqrt(z^2-(f1x0*f1x1));
    end
mu=(f1x1+w-z)/(f1x1-f1x0+2*w); % to find mu
    if mu<0
       xn=x1;
       %sprintf('x_min=%f', xn)
       plot(xn,fxn,'ro')
       hold on
       plot (xn,f1xn,'bo')
elseif mu>=0 && mu<= 1
       xn=x1-mu*(x1-x0);
       %sprintf('x_min=%f', xn)
       xlabel('x'); 
       ylabel('f(x)');
       grid on 
       plot(xn,fxn,'bo')
       hold on
       plot (xn,f1xn,'ro')       
    else
       xn=x0;
       %sprintf('x_min=%f', xn)
       title('Cubic search');
       xlabel('x'); 
       ylabel('f(x)');
       grid on 
       plot(xn,fxn,'go')
       hold on
       plot (xn,f1xn,'ro')
    end
fxn=f(xn);
f1xn=f1(xn);
i=i+1;
%eval(xn);
end
disp('The minimum value of the function is: ');
disp(eval(xn));
disp('The minimum of iterations is: ');
disp(i);
 %disp(['The optimum value of x is ' ,num2str(xn)]);
 %-disp(['The number of iteration is ',num2str(i)]);