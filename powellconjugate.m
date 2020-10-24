%% optimization using powell's conjugate directions 
function powellconjugate
syms x1 x2 x3 x4 a a1 a2 a3 a4;
x=[x1;x2;x3;x4];
i=0;                        % iteration counter 
f = (-4*x1+8*x2-2*x3+2*x4)^2+(x1-3*x2+x3-3*x4)^2+(-x1+x2+4*x3)^2+2*x1*x3+4*x2;
%X1=input('Enter the first element of input vector ');
%X2=input('Enter the second element of input vector ');
%X3=input('Enter the third element of input vector ');
%X4=input('Enter the fourth element of input vector ');
X1=1;X2=1;X3=1;X4=1;
x01=[X1;X2;X3;X4];
fx01=subs(f,(x),(x01));
s11=[1;0;0;0];              % choosing the direction 
x11=x01+a*s11;
fx11=vpa(subs(f,(x),(x11)));
a=solve(diff(fx11)==0,a);   % minimizing alpha 
x11=x01+a*s11;
fx11=vpa(subs(f,(x),(x11)));
s21=[0;1;0;0];
x21=x11+a1*s21;
fx21=vpa(subs(f,(x),(x21)));
a1=solve(diff(fx21)==0,a1);  % minimizing alpha1 
x21=x11+a1*s21;
fx21=vpa(subs(f,(x),(x21)));
s31=[0;0;1;0];
x31=x21+a2*s31;
fx31=vpa(subs(f,(x),(x31)));
a2=solve(diff(fx31)==0,a2);  % minimizing alpha2
x31=x21+a2*s31;
fx31=vpa(subs(f,(x),(x31)));
s41=[0;0;0;1];
x41=x31+a3*s41;
fx41=vpa(subs(f,(x),(x41)));
a3=solve(diff(fx41)==0,a3);  % minimizing alpha3
x41=x31+a3*s41;
fx41=vpa(subs(f,(x),(x41)));
sp1=x41-x01;
x02=x01+a4*sp1;
fx02=vpa(subs(f,(x),(x02)));
simplify(fx02);
a4=solve(diff(fx02)==0,a4);  % minimizing alpha4
x02=x01+a4*sp1;
fx02=vpa(subs(f,(x),(x02)));
fd1=fx01-fx11;
fd2=fx11-fx21;
fd3=fx21-fx31;
fd4=fx31-fx41;
% finding the maximum difference among the consecutive function values 
    if fd1>fd2&&fd1>fd3&&fd1>fd4
       fdmax=fd1;
elseif fd2>fd1&&fd2>fd3&&fd2>fd4
       fdmax=fd2;
elseif fd3>fd1&&fd3>fd2&&fd3>fd4
       fdmax=fd3;
  else
       fdmax=fd4;
   end
condition=sqrt((fx01-fx02)/(fdmax));
n=double(norm(x02-x01));
func_eval=11;
deri_eval=5;

   if abs(a4)<abs(condition)            % condition for finding the direction 
      while norm(x02-x01)>0.0001        % convergence criterion 
            clear a a1 a2 a3 a4 
            syms a a1 a2 a3 a4
             x01=x02;
             fx01=vpa(subs(f,(x),(x01)));
             s11=[1;0;0;0];
             x11=x01+a*s11;
             fx11=vpa(subs(f,(x),(x11)));
             a=solve(diff(fx11)==0,a);
             x11=x01+a*s11;
             fx11=vpa(subs(f,(x),(x11)));
             s21=[0;1;0;0];
             x21=x11+a1*s21;
             fx21=vpa(subs(f,(x),(x21)));
             a1=solve(diff(fx21)==0,a1);
             x21=x11+a1*s21;
             fx21=vpa(subs(f,(x),(x21)));
             s31=[0;0;1;0];
             x31=x21+a2*s31;
             fx31=vpa(subs(f,(x),(x31)));
             a2=solve(diff(fx31)==0,a2);
             x31=x21+a2*s31;
             fx31=vpa(subs(f,(x),(x31)));
             s41=[0;0;0;1];
             x41=x31+a3*s41;
             fx41=vpa(subs(f,(x),(x41)));
             a3=solve(diff(fx41)==0,a3);
             x41=x31+a3*s41;
             fx41=vpa(subs(f,(x),(x41)));
             sp1=x41-x01;
             x02=x01+a4*sp1;
             fx02=vpa(subs(f,(x),(x02)));
             simplify(fx02);
             a4=solve(diff(fx02)==0,a4);
             x02=x01+a4*sp1;
             fx02=vpa(subs(f,(x),(x02)));
             fd1=fx01-fx11;
             fd2=fx11-fx21;
             fd3=fx21-fx31;
             fd4=fx31-fx41;
          if fd1>fd2&&fd1>fd3&&fd1>fd4
             fdmax=fd1;
      elseif fd2>fd1&&fd2>fd3&&fd2>fd4
             fdmax=fd2;
      elseif fd3>fd1&&fd3>fd2&&fd3>fd4
             fdmax=fd3;
        else
             fdmax=fd4;
         end
             condition=sqrt((fx01-fx02)/(fdmax));
             i=i+1;
             func_eval=func_eval+11;
             deri_eval=deri_eval+5;
      end
  else
      if fdmax==fd1
         s11=sp1;
  elseif fdmax==fd1
         s21=sp1;
  elseif fdmax==fd2
         s31=sp1;
    else 
         s41=sp1;
      end
     while norm(x02-x01)<0.0001
         if s11==sp1
            s21=[0;1;0;0];
            s31=[0;0;1;0];
            s41=[0;0;0;1];
     elseif s21==sp1
            s11=[1;0;0;0];
            s31=[0;0;1;0];
            s41=[0;0;0;1];
     elseif s31==sp1
            s21=[0;1;0;0];
            s11=[1;0;0;0];
            s41=[0;0;0;1];
       else 
            s21=[0;1;0;0];
            s31=[0;0;1;0];
            s11=[1;0;0;0];
        end
x01=x02;
fx01=vpa(subs(f,(x),(x01)));
x11=x01+a*s11;
fx11=vpa(subs(f,(x),(x11)));
a=solve(diff(fx11)==0,a);
x11=x01+a*s11;
fx11=vpa(subs(f,(x),(x11)));
x21=x11+a1*s21;
fx21=vpa(subs(f,(x),(x21)));
a1=solve(diff(fx21)==0,a1);
x21=x11+a1*s21
fx21=vpa(subs(f,(x),(x21)));
x31=x21+a2*s31;
fx31=vpa(subs(f,(x),(x31)));
a2=solve(diff(fx31)==0,a2);
x31=x21+a2*s31;
fx31=vpa(subs(f,(x),(x31)));
x41=x31+a3*s41;
fx41=vpa(subs(f,(x),(x41)));
a3=solve(diff(fx41)==0,a3);
x41=x31+a3*s41;
fx41=vpa(subs(f,(x),(x41)));
sp1=x41-x01;
x02=x01+a4*sp1;
fx02=vpa(subs(f,(x),(x02)));
simplify(fx02);
a4=solve(diff(fx02)==0,a4);
x02=x01+a4*sp1;
fx02=vpa(subs(f,(x),(x02)));
fd1=fx01-fx11;
fd2=fx11-fx21;
fd3=fx21-fx31;
fd4=fx31-fx41;
% finding maximum difference among the consecutive function values 
        if fd1>fd2&&fd1>fd3&&fd1>fd4     
           fdmax=fd1;
    elseif fd2>fd1&&fd2>fd3&&fd2>fd4
           fdmax=fd2;
    elseif fd3>fd1&&fd3>fd2&&fd3>fd4
           fdmax=fd3;
      else
           fdmax=fd4;
        end
        condition=sqrt((fx01-fx02)/(fdmax));
        i=i+1;
     end
  end
disp(['Number of iterations = ',num2str(i)]);
disp('The optimum value of x is ');
double(x02)
disp(['Number of function evaluation = ',num2str(func_eval)]);
disp(['Number of derivative evaluation = ',num2str(deri_eval)]);