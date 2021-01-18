function Main_grad_descent_2
% using dorpri8 for grad descent
clear variables
close all
clc

a = Class_bumpy_funs;
X1 = [];

x = linspace(-3,3,1e2)';
y = linspace(-3,3,1e2)';


[xx,yy] = meshgrid(x,y);


x0 = [0;-1];

% FAILS: (ROLLS DOWN THE SLOPE)
% x0 = [-1.840525771547726
%   -1.782448121543857];


    options1 = struct('outputfcn',@outFunMAP);
 
figure
[T,X] = ode8d(@(t,x) -a.grad_Peaks(x(1),x(2)),[0,10],x0,options1);
X = X.';
    

l = [0;cumsum(sqrt(sum(diff(X,1,1).^2,2)))];
    ff = peaks(X(:,1),X(:,2));

figure
    peaks
    hold on
    plot(X(:,1),X(:,2),'r','linewidth',2)
    plot(X(1,1),X(1,2),'or')
    plot(X(end,1),X(end,2),'xr')


figure
    plot(l,ff)
    yline(min(min(peaks)))

keyboard

    

     function status = outFun1(t,xv,flag)
        
        status = 0;
  
        if any(strcmp(flag,{'init','done'})), return; end
        
        X1 = [X1 xv];
        
        cla
        peaks;
        hold on
        plot(X1(1,:),X1(2,:),'.-r')
        
        
        
        title(['t= ' num2str(t)])
        drawnow
                 
    end
  
    
     function status = outFunMAP(t,y,flag)
        
        if ~isempty(y), y = y(1:2); end
        
        
        status = odeplot(t,y,flag);
        title(['t= ' num2str(t)])
%          keyboard
    end
    


end



