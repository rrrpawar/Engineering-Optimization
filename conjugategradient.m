%% Optimize using Conjugate gradient Method
function conjugategradient
syms x1 x2 x3 x4 x00 x01 a
f = (-4*x1+8*x2-2*x3+2*x4)^2+(x1-3*x2+x3-3*x4)^2+(-x1+x2+4*x3)^2+2*x1*x3+4*x2;
x=[x1;x2;x3;x4];
i=0;                                 % iteration counter
df = gradient(f,[x1,x2,x3,x4]);      % to find the gradient 
%X1 = input('Enter the first element of input vector');
%X2 = input('Enter the second element of input vector');
%X3 = input('Enter the third element of input vector');
%X4 = input('Enter the fourth element of input vector');
X1=1;X2=1;X3=1;X4=1;
x00 = [X1;X2;X3;X4];
f0 = vpa(subs(f,(x),(x00)));        
grad = -vpa(subs(df,(x),(x00)));            
s0 = grad;
x01 = x00+a*s0;
f01 = vpa(subs(f,(x),(x01)));
a = solve(diff(f01)==0,a);              
x01 = x00+a*s0;
f01 = vpa(subs(f,(x),(x01)));
grad_1 = -vpa(subs(df,(x),(x01)));
func_eval=3;
grad_eval=2;
while ((norm(grad)>0.0001)||(abs(f01-f0)>0.0001)||(norm(grad_1)>0.0001))  % convergence criterion
    clear a
    syms a 
    grad = -vpa(subs(df,(x),(x00)));
    grad_1 = -vpa(subs(df,(x),(x01)));
    b = (norm(-grad_1)^2)/(norm(-grad)^2);      
    s1 = grad_1+b*s0;
    x00 = x01;
    x01 = x00+a*s1;
    f01 = vpa(subs(f,(x),(x01)));
    a = solve(diff(f01)==0,a);                
    x01 = x00+a*s1;
    grad_1 = -vpa(subs(df,(x),(x01)));
    f0 = vpa(subs(f,(x),(x00)));
    f01 = vpa(subs(f,(x),(x01)));
    i = i+1;
    func_eval=func_eval+3;                  
    grad_eval=grad_eval+3;                                  
    s0 = s1;                                  
end
disp(['Number of iterations = ',num2str(i)]);
disp('The optimum value of x is ');
double(x00)
disp(['Number of function evaluation = ',num2str(func_eval)]);
disp(['Number of gradient evaluation = ',num2str(grad_eval)]);

   
   