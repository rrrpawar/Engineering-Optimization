
function steepestdescent
syms x1 x2 x3 x4 a0
%x01=input('Enter the first element of input vector ');
%x02=input('Enter the second element of input vector ');
%x03=input('Enter the third element of input vector ');
%x04=input('Enter the fourth element of input vector ');
x01=1;x02=1;x03=1;x04=1;
x0=[x01;x02;x03;x04];
x=[x1;x2;x3;x4];
fun=@(x1,x2,x3,x4) (-4*x1+8*x2-2*x3+2*x4)^2+(x1-3*x2+x3-3*x4)^2+(-x1+x2+4*x3)^2+2*x1*x3+4*x2;
delfx=gradient(fun,[x1,x2,x3,x4]);       % finding the gradient
d=vpa(subs(delfx,{x1,x2,x3,x4},{x0(1),x0(2),x0(3),x0(4)}));
i=0;
normf=norm(d);
fun_eval=1;
delfx_eval=4;
while normf > 0.0001 
      x1 = x0 - a0.*d;
      fx1=fun(x1(1),x1(2),x1(3),x1(4));
      fun_eval=fun_eval+1;
      dfa0=diff(fx1,a0);
      delfx_eval=delfx_eval+5;
      a0=solve(dfa0==0,a0);
      x1 = x0 - a0.*d;
      clear a0;
      syms a0;
      x0=x1;
      d=vpa(subs(delfx,(x),(x0)));
      normf=norm(d);
      i=i+1;
end
disp(['Number of iterations = ',num2str(i)]);
disp('The optimum value of x is ');
double(x0)
disp(['Number of function evaluation = ',num2str(fun_eval)]);
disp(['Number of derivative evaluation = ',num2str(delfx_eval)]);
