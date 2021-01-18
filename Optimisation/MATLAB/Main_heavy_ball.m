function Main_heavy_ball
% using ode15s for grad descent
clear variables
close all
clc

a = Class_bumpy_funs;
b = 1;
XV = [];

x0 = [-2;2];
x0 = [0;-1];
x0 = [-1;0];
x0 = [0;0];
v0 = [0;0];rand(2,1);a.grad_Peaks(x0(1),x0(2));

    options1 = odeset('outputfcn',@outFunMAP);
    options1 = odeset('outputfcn',@outFun1);

    
figure
[T,X] = ode15s(@rhs,[0,10],[x0;v0],options1);


l = [0;cumsum(sqrt(sum(diff(X,1,1).^2,2)))];

figure
    peaks
    hold on
    ff = peaks(X(:,1),X(:,2));
    plot3(X(:,1),X(:,2),ff,'r','linewidth',2)
    plot3(X(1,1),X(1,2),ff(1),'or')
    plot3(X(end,1),X(end,2),ff(end),'xr')


figure
    plot(l,ff,'.-')
    yline(min(min(peaks)))

    
figure(3)
hold on
    plot(T,ff,'.-')
keyboard


    function status = outFun1(t,xv,flag)
        
        status = 0;
  
        if any(strcmp(flag,{'init','done'})), return; end
        
        XV = [XV xv];
        
%         ff = peaks(XV(1,:),XV(2,:));
 v = [xv(2),xv(3)]./norm([xv(2),xv(3)])/3;       
        cla
        peaks;
        hold on
        plot(XV(1,:),XV(2,:),'.-r')
        
        quiver(xv(1),xv(2),v(1),v(2),0,'linewidth',1,'color','k')
        
        
        title(['t= ' num2str(t)])
        drawnow
                 
    end


    function out = rhs(t,xv)
        N = length(xv)/2;
        
        
        
        x = xv(1:N);
        v = xv(N+1:end);
        
        out = [v;
               -b*v-a.grad_Peaks(x(1),x(2))];
        
        
    end

    
    
    
     function status = outFunMAP(t,y,flag)
        
        if ~isempty(y), y = y(1:2); end
        
        
        status = odeplot(t,y,flag);
        title(['t= ' num2str(t)])
         
    end
    


end



