%% optimization using simplex method
function simplex_method
syms x1 x2 x3 x4
x=[x1;x2;x3;x4];
iter=0;
f(x1,x2,x3,x4) = (-4*x1+8*x2-2*x3+2*x4)^2+(x1-3*x2+x3-3*x4)^2+(-x1+x2+4*x3)^2+2*x1*x3+4*x2;
p=1/(4*sqrt(2))*(sqrt(5)+4-1);
q=1/(4*sqrt(2))*(sqrt(5)-1);
%x01=input('Enter the first element of input vector');
%x02=input('Enter the second element of input vector');
%x03=input('Enter the third element of input vector');
%x04=input('Enter the fourth element of input vector');
x01=1;x02=1;x03=1;x04=1;
x0=[x01;x02;x03;x04];
fx00=vpa(subs(f,(x),(x0)));
s1=[1;0;0;0]; s2=[0;1;0;0];
s3=[0;0;1;0]; s4=[0;0;0;1];
x01=x0+p*s1+q*(s2+s3+s4);  % find equidistant points from the starting point
fx01=vpa(subs(f,(x),(x01)));
x02=x0+p*s2+q*(s1+s3+s4);
fx02=vpa(subs(f,(x),(x02)));
x03=x0+p*s3+q*(s1+s2+s4);
fx03=vpa(subs(f,(x),(x03)));
x04=x0+p*s4+q*(s1+s3+s2);
fx04=vpa(subs(f,(x),(x04)));
if fx00>fx01 && fx00>fx02 && fx00>fx03 && fx00>fx04 % to find xh
    xh=x0;
    xs1=x01; xs2=x02;
    xs3=x03; xs4=x04;
elseif fx01>fx00 && fx01>fx02 && fx01>fx03 && fx01>fx04
    xh=x01;
    xs1=x02; xs2=x03;
    xs3=x04; xs4=x0;
elseif fx02>fx00 && fx02>fx01 && fx02>fx03 && fx02>fx04
    xh=x02;
    xs1=x01; xs2=x03;
    xs3=x04; xs4=x0;
elseif fx03>fx00 && fx03>fx01 && fx03>fx02 && fx03>fx04
    xh=x03;
    xs1=x02; xs2=x01;
    xs3=x04; xs4=x0; 
else
    xh=x04;
    xs1=x02; xs2=x03;
    xs3=x01; xs4=x0;
end
fx_h=vpa(subs(f,(x),(xh)));
centroid=1/(4)*(xs1+xs2+xs3+xs4); % to find centroid 
fc=vpa(subs(f,(x),(centroid)));
xr=centroid+ 1*(centroid-xh);
fx_r=vpa(subs(f,(x),(xr)));
cond=sqrt(1/5*((fx00-fc)^2+(fx01-fc)^2+(fx02-fc)^2+(fx03-fc)^2+(fx04-fc)^2));
fx_x=[fx00;fx01;fx02;fx03;fx04];
f1=min(fx_x);
func_eval=9;
while cond>0.0001                % convergence criteria 
if fx_r>f1&&fx_r<fx_h            % to replace xh by xr
    xh=xr;
elseif fx_r<f1 
    xe=centroid+2*(xr-centroid); % to replace xh by xe 
    fe=vpa(subs(f,(x),(xe)));
    if fe<fx_r
        xh=xe;
    else
        xh=xr;
    end
elseif fx_r>=fx_h
    xc=centroid+0.5*(xh-centroid);
    fx_c=vpa(subs(f,(x),(xc)));
 % to find the best point xl with minimun function value
 if fx00<fx01 && fx00<fx02 && fx00<fx03 && fx00<fx04   
    x1=x0;
    x11=x01; x12=x02;
    x13=x03; x14=x04;
elseif fx01<fx00 && fx01<fx02 && fx01<fx03 && fx01<fx04
    x1=x01;
    x11=x02; x12=x03;
    x13=x04; x14=x0;   
elseif fx02<fx00 && fx02<fx01 && fx02<fx03 && fx02<fx04
    x1=x02;
    x11=x01; x12=x03;
    x13=x04; x14=x0;  
elseif fx03<fx00 && fx03<fx01 && fx03<fx02 && fx03<fx04
    x1=x03;
    x11=x02; x12=x01;
    x13=x04; x14=x0;   
else
    x1=x04;
    x11=x02; x12=x03;
    x13=x01; x14=x0;   
end
    if fx_c>fx_h
        x0=x0+1/2*(x1-x0);
        x01=x01+1/2*(x1-x01);
        x02=x02+1/2*(x1-x02);
        x03=x03+1/2*(x1-x03);
        x04=x04+1/2*(x1-x04);
    else 
        xh=xc; % to replace xh by xc 
    end
end
    x0=xs1; x01=xs2;
    x02=xs3; x03=xs4; 
    x04=xh;
    fx00=vpa(subs(f,(x),(x0)));
    fx01=vpa(subs(f,(x),(x01)));
    fx02=vpa(subs(f,(x),(x02)));
    fx03=vpa(subs(f,(x),(x03)));
    fx04=vpa(subs(f,(x),(x04)));
% to find the highest function value (worst point)
   if fx00>fx01&&fx00>fx02&&fx00>fx03&&fx00>fx04
    xh=x0;
    xs1=x01; xs2=x02;
    xs3=x03; xs4=x04;
elseif fx01>fx00 && fx01>fx02 && fx01>fx03 && fx01>fx04
    xh=x01;
    xs1=x02; xs2=x03;
    xs3=x04; xs4=x0;
elseif fx02>fx00 && fx02>fx01 && fx02>fx03 && fx02>fx04
    xh=x02;
    xs1=x01; xs2=x03;
    xs3=x04; xs4=x0;
elseif fx03>fx00 && fx03>fx01 && fx03>fx02 && fx03>fx04
    xh=x03;
    xs1=x02; xs2=x01;
    xs3=x04; xs4=x0; 
else
    xh=x04;
    xs1=x02; xs2=x03;
    xs3=x01; xs4=x0;
end
fx_h=vpa(subs(f,(x),(xh)));
centroid=1/(4)*(xs1+xs2+xs3+xs4);
fc=vpa(subs(f,(x),(centroid)));
xr=centroid+ 1*(centroid-xh);
fx_r=vpa(subs(f,(x),(xr)));
cond=sqrt(1/5*((fx00-fc)^2+(fx01-fc)^2+(fx02-fc)^2+(fx03-fc)^2+(fx04-fc)^2));
fx_x=[fx00;fx01;fx02;fx03;fx04];
f1=min(fx_x);
func_eval=func_eval+11;
iter=iter+1;
end
disp(['Number of iterations = ',num2str(iter)]);
disp('The optimum value of x from simplex method is ');
double(xh)
disp(['Number of function evaluation = ',num2str(func_eval)]);