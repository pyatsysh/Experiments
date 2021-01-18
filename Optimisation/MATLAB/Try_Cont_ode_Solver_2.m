function Try_Cont_ode_Solver_2
% working version
% build a stiff solver by using continuation to dertermine dt
% t is parameter
% 
% Motivation: ode15s is too magical. Can I compute delta t using
% continuation??
% 
% y'=-y^2, y(0)=1
% 
% Notes:
%       * there should not be any principled difference between implicit
%       and explicit schemes
%       * in continuation view, y and t are the unknown and dt is parameter



%   * ALGORITHM NEVER REDUCES t -- MAYBE, MAKES SENSE, BECAUSE IT TREATS IT AS A
%   PARAMETER. LET'S REFORMULATE WITH t AS A PARAMETER


% 
% TO SOLVE STIFF
%   * Higer-order Runge-Kutta, implicit
%   * Lower  tolerance
%   * Analytic jacobian into continuation

clear variables
% close all
clc

y_prime = 0;

y0 = [2;0];

y_old = y0;

dt0 = 0.1;

T = 0;
DT = dt0;
Y = y0;

    
    hFig = figure;
    contHDFbeta2_1(@rhs_cont,y_old,dt0,inf,...
                    'monitorfcn',@MonitorFun,...
                    'direction',1,...
                    'step',0.1,...
                    'adaptive',[.1,.1,.1],...
                    'tolfun',1e-8,...
                    'jacobian',[],...
                    'stopbtn',hFig);


   keyboard 

   
       contHDFbeta2_1(@rhs_cont,y_old,dt0,inf,...
                    'monitorfcn',@MonitorFun,...
                    'step',0.1,...
                    'adaptive',[],...
                    'jacobian',[],...
                    'stopbtn',hFig);


    function dydt = rhs_ode(t,y)
        % ode15s benchmark problem
        dydt = [y(2); 1000*(1-y(1)^2)*y(2)-y(1)];
    end

    
    function out = rhs_cont(y,dt)
        % euler method for derivatives:
        % forward  euler (explicit): y_prime = rhs_ode(t,y_old)
        % backward euler (implicit): y_prime = rhs_ode(t,y_new)
        
        y_new = y;
        t = T(end)+dt;
        
        y_prime = rhs_ode(t,y_new);
        out = y_new-dt*y_prime-y_old; 
        
    end


    function [MonVal,status] = MonitorFun(y,dt)
        
        MonVal = [];
        status = [];
        
        y_new = y;
        t = T(end)+dt;
        
        y_old = y_new;
        
T = [T;t];
DT = [DT;dt];
Y = [Y y_old];
        
       cla
        plot(T,Y(1,:),'x-')
        title(['dt= ' num2str(dt)])
        
        drawnow
% keyboard        
    end




end



