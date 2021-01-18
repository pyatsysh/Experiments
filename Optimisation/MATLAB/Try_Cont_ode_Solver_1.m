function Try_Cont_ode_Solver_1
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
%       * implementing simplest approximation for derivatives -- Euler 
%         (implicit: forward and explit: backward)
%       * there should not be any principled difference between implicit
%       and explicit schemes
%       * in continuation view, y and t are the unknown and dt is parameter

clear variables
close all
clc

y0 = 1;

y_old = y0;

dt0 = 0.1;

T = 0;
Y = y0;

    
    contHDFbeta2_1(@rhs_cont,y_old,dt0,100,...
                    'monitorfcn',@MonitorFun,...
                    'step',0.1,...
                    'adaptive',[],...
                    'jacobian',[],...
                    'stopbtn',figure);


    function out = rhs_ode(t,y)
        out = -y.^2;
    end
    
    function out = rhs_cont(y,dt)
        % euler method for derivatives:
        % forward  euler (explicit): rhs_ode(t,y_old)
        % backward euler (implicit): rhs_ode(t,y_new)
        
        
        y_new = y;
        t = T(end)+dt;
        
        out = y_new-dt*rhs_ode(t,y_old)-y_old; % for explicit scheme by forward euler rhs_ode(t,y_old)
        
    end


    function [MonVal,status] = MonitorFun(y,dt)
        
        MonVal = [];
        status = [];
        
        y_new = y;
        t = T(end)+dt;
        
        y_old = y_new;
        
T = [T;t];
Y = [Y y_old];
        
       cla
        plot(T,Y,'x-')
        hold on
        plot(T,1./(T+1))
        title(['dt= ' num2str(dt)])
        
keyboard        
    end




end



