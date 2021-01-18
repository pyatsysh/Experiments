function Main_Stoch_grad_descent
% works

clear variables
close all
clc

a = Class_bumpy_funs;
X1 = [];

x = linspace(-3,3,1e2)';
y = linspace(-3,3,1e2)';


[xx,yy] = meshgrid(x,y);


x0 = (rand(2,1)-.5)*2*3;

    options_e = struct('step',.01,...
                       'outputfcn',@outFun1,...
                       'noisestrength',1);

% optimal values                   
Fo = a.Peaks(x0(1),x0(2));                   
Xo = x0;

   [T,X] = ode1eS(@(t,x) a.grad_Peaks(x(1),x(2)),[0,10],x0,options_e);

   
   figure
   plot(Fo,'.-')
   
   
keyboard
l = [0;cumsum(sqrt(sum(diff(X,1,1).^2,2)))];
    ff = peaks(X(:,1),X(:,2));

figure
    plot(l,ff)
    
    
figure
    peaks
    hold on
    plot(X(:,1),X(:,2),'r','linewidth',2)
    plot(X(1,1),X(1,2),'or')
    plot(X(end,1),X(end,2),'xr')


    yline(min(min(peaks)))

keyboard

    

    function status = outFun1(t,x,flag)
        
        
%         if strcmp(flag,'init'), return; end
        subplot(3,1,1)
            status = odeplot(t,x,flag);
 
        
        X1 = [X1 x];
        
        
        f1 = a.Peaks(x(1),x(2));
        ispl = false;
        if f1>Fo(end),
            Fo = [Fo;f1];
            Xo = [Xo x];
            ispl = true;
        end
        
        subplot(3,1,2)
        cla
            peaks;
            hold on
            plot(X1(1,:),X1(2,:),'r')
        
        subplot(3,1,3)
        hold on
            if ispl,
                plot(t,f1,'.')
            end
        
        subplot(3,1,3)    
        title(['t= ' num2str(t)])
 
%         title(['t= ' num2str(t)])
        
        
        drawnow
%         keyboard
    end

    
    
     function status = outFunMAP(t,y,flag)
        
        if ~isempty(y), y = y(1:2); end
        
        
        status = odeplot(t,y,flag);
        title(['t= ' num2str(t)])
         
    end
    


end



