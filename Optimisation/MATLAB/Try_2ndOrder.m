function Try_2ndOrder
% maybe, I can't solve 2nd order odes?
clear variables
close all
clc

    options1 = odeset('outputfcn',@outFunMAP);

    
    
x0 = 0;
v0 = -1;
    
figure
[T,X] = ode15s(@rhs,[0,10],[x0;v0],options1);

figure
plot(T,X(:,1),'.')
hold on
plot(T,T.^2+exp(-T))

keyboard

    function out = rhs(t,xv)
        N = length(xv)/2;
        
        x = xv(1:N);
        v = xv(N+1:end);
        
        out = [v;
               2+exp(-x)];
        
        
    end

    function status = outFunMAP(t,y,flag)
        
        if ~isempty(y), y = y(1:2); end
        
        
        status = odeplot(t,y,flag);
        title(['t= ' num2str(t)])
         
    end
    


end









