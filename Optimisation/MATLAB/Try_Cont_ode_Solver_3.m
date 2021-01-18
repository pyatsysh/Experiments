function Try_Cont_ode_Solver_3
% REFORMULATED WITH t AS A PARAMETER
% PROBLEM: INCRESING STIFFNESS K > 5 BREAKS IT -- TRY HIGHER ORDER
% RUNGE-KUTTA? 

% 
% TO SOLVE STIFF
%   * Higer-order Runge-Kutta, implicit
%   * Lower  tolerance
%   * Analytic jacobian into continuation

clear variables
close all
clc

K = 5; % stiffness (try 1, 1e1, 1e2, 1e3): K=1e3 is benchmark


y_prime = 0;

y0 = [2;0];

y_old = y0;

dt0 = 0.1;
t0 = 0;

T = t0;
DT = dt0;
Y = y0;

[t_exact,y_exact] = ode15s(@rhs_ode,[0,10],y0);
plot(t_exact,y_exact(:,1),'-or')
%   
keyboard

    hFig = figure;
    contHDFbeta2_1(@rhs_cont,y_old,t0,inf,...
                    'monitorfcn',@MonitorFun,...
                    'direction',1,...
                    'step',0.1,...
                    'adaptive',[],...
                    'tolfun',1e-8,...
                    'jacobian','on',...
                    'stopbtn',hFig);


   keyboard 

   
       contHDFbeta2_1(@rhs_cont,y_old,T(end),inf,...
                    'monitorfcn',@MonitorFun,...
                    'direction',1,...
                    'step',1e-1,...
                    'adaptive',[],...
                    'tolfun',1e-5,...
                    'jacobian',[],...
                    'stopbtn',hFig);

                
                
       contHDFbeta2_1(@rhs_cont,y_old,T(end),inf,...
                    'monitorfcn',@MonitorFun,...
                    'direction',1,...
                    'step',0.001,...
                    'adaptive',[.1,.1,.1],...
                    'tolfun',1e-8,...
                    'jacobian',[],...
                    'stopbtn',hFig);
                
 
   function out = rhs_ode1(t,y)
        out = -y.^2;
    end

    function dydt = rhs_ode(t,y)
        % ode15s benchmark problem
        dydt = [y(2); 
               K*(1-y(1)^2)*y(2)-y(1)];
    end

    function J = J_rhs_ode(t,y)
        
        J = [0,                     1;
             -K*y(2)*2*y(1)-1,   K*(1-y(1)^2)];
        
        
    end
    



    function [out,J] = rhs_cont(y,t)
        % Runge-Kutta order 1 (euler method) for derivatives:
        % forward  euler (explicit): y_prime = rhs_ode(t,y_old)
        % backward euler (implicit): y_prime = rhs_ode(t,y)
        
        
        
        dt = t-T(end);
        
        y_prime = rhs_ode(t,y);
        out = y-dt*y_prime-y_old; 
        
        
        % jac w.r.t. y
        if nargout == 2, 
            
            J = eye(length(y))-dt*J_rhs_ode(t,y);
            
        end
    end

    
   
    

    function [MonVal,status] = MonitorFun(y,t)
        
        
        
        MonVal = [];
        status = [];
        
        if t>10, status = NaN; end
        
        dt = t-T(end);
        
        y_old = y;
        
T  = [T;t];
DT = [DT;dt];
Y  = [Y y_old];
        
       cla
        plot(T,Y(1,:),'x-')
        hold on
        plot(t_exact,y_exact(:,1),'.-')
        title(['dt= ' num2str(dt)])
        
        drawnow
% keyboard        
    end




end



