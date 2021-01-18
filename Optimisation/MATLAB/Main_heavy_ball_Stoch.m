function Main_heavy_ball_Stoch
% doesn't work
clear variables
close all
clc

a = Class_bumpy_funs;
b = 1;
X1 = [];

    options_e = struct('step',.1,...
                       'outputfcn',@outFun1,...
                       'noisestrength',.1);

% optimal values 

x0 = [0;-1];
figure
Fo = a.Peaks(x0(1),x0(2));                   
Xo = x0;
v0 = [0;0];

   [T,X] = ode1eS(@rhs,[0,10],[x0;v0],options_e);

   
keyboard

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


  
    function status = outFun1(t,x,flag)
        
        status = 0;
        
        if strcmp(flag,'init'), return; end
        
        X1 = [X1 x];
        
        
        f1 = a.Peaks(x(1),x(2));
        ispl = false;
        if f1>Fo(end),
            Fo = [Fo;f1];
            Xo = [Xo x(1:2)];
            ispl = true;
        end
        
        subplot(2,1,1)
        cla
            peaks;
            hold on
            plot(X1(1,:),X1(2,:),'r')
        
        subplot(2,1,2)
        hold on
            if ispl,
                plot(t,f1,'.')
            end
        
        title(['t= ' num2str(t)])
        
        
        drawnow
%         keyboard
    end


    function out = rhs(t,xv)
        N = length(xv)/2;
        
        x = xv(1:N);
        v = xv(N+1:end);
        
        out = [v;
               -b*v+a.grad_Peaks(x(1),x(2))];
        
        
    end

    
    
    
     function status = outFunMAP(t,y,flag)
        
        if ~isempty(y), y = y(1:2); end
        
        
        status = odeplot(t,y,flag);
        title(['t= ' num2str(t)])
         
    end
    


end



