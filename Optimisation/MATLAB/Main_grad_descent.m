function Main_grad_descent
% using ode15s for grad descent
clear variables
close all
clc

a = Class_bumpy_funs;
X1 = [];

x = linspace(-3,3,1e2)';
y = linspace(-3,3,1e2)';


[xx,yy] = meshgrid(x,y);


x0 = [0;-1];

    options1 = odeset('outputfcn',@outFunMAP);
    
   [T,X] = ode15s(@(t,x) -a.grad_Peaks(x(1),x(2)),[0,1000],x0,options1);

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

    

    function status = outFun1(t,x,flag)
        
        status = 0;
        
        if strcmp(flag,'init'), return; end
        
        X1 = [X1 x];
        
%         ff = peaks(XV(1,:),XV(2,:));
%  v = [xv(2),xv(3)]./norm([xv(2),xv(3)])/3;       
        cla
        peaks;
        hold on
        plot(X1(1,:),X1(2,:),'r')
%         quiver(xv(1),xv(2),v(1),v(2),0,'linewidth',1,'color','k')
        
%         plot3(XV(1,:),XV(2,:),ff,'r','linewidth',2)

%         quiver3(xv(1),xv(2),ff(end),xv(3),xv(4),0,'linewidth',1,'color','k')
%         hold off
        
        
        
        title(['t= ' num2str(t)])
        
        drawnow
%         keyboard
         
    end

    
    
     function status = outFunMAP(t,y,flag)
        
        if ~isempty(y), y = y(1:2); end
        
        
        status = odeplot(t,y,flag);
        title(['t= ' num2str(t)])
         
    end
    


end



