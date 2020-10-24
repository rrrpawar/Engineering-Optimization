function master=mainprogram()
syms x;
list = {'I want to enter a new funtion ','f(x)=10+(x-4)^2+10*exp(-x/4-2)-3*sin(0.2*x*pi)'};
[index,~] = listdlg('ListString',list);
if (index==1)
    prompt={'Enter the function'};
    title='Input function';
    solution=inputdlg(prompt,title);
    
elseif (index==2)
    solution=10+(x-4)^2+10*exp(-x/4-2)-3*sin(0.2*x*pi);
end 

list = {'I want to enter new values','x1=-5,x2=20,dx=2,a=-5,b=20,convergence criteria=0.00001'};
[index2,~] = listdlg('ListString',list);
if (index2==1)
    prompt={'Enter x1'};
    title='Input lower bound';
    x1=inputdlg(prompt,title);
    prompt={'Enter x2'};
    title='Input upper bound';
    x2=inputdlg(prompt,title);
    prompt={'Enter stepsize'};
    title='Input step size';
    dx=inputdlg(prompt,title);
    prompt={'Enter epsilon'};
    title='Input convergence criteria';
    e=inputdlg(prompt,title);

elseif (index2==2)
    x1=-5;x2=20;dx=2;e=0.00001;a=-5;b=20; 
end 

list = {'Select the 1-D minimization method','exit'};
[index3,~] = listdlg('Liststring',list);
if (index3==1)
list = {'Interval Halving','Golden section search','Bisection','Powell','Cubic'};
elseif (index3==2)
       solution=inputdlg(prompt,title);
end

[index4,~] = listdlg('ListString',list);
if (index4==1)
    a= intervalhalving_method(x1,x2,e,solution);
    
elseif (index4==2)
    a= goldensection_method(x1,x2,e,solution);
    
elseif (index4==3)
    a= bisection_method(x1,x2,e,solution);
    
elseif(index4==4)
    a= powell_method(x1,x2,e,solution);
    
elseif(index4==5)
    a=cubicsearch_method(x1,x2,e,solution);
end    

end