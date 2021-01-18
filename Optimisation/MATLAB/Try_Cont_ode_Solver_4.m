function Try_Cont_ode_Solver_4
% t IS PARAMETER
% INCREASING ORDER OF SOLVER TO RADAU IIA
% see this stack exchange answer for algorithm:
% https://math.stackexchange.com/questions/2881557/how-to-write-the-radau-2nd-order-methods-butcher-s-table
% 
% 
% TO SOLVE STIFF
%   * If it gets fucked, take a step back, reduce stepsize and re-do, don't
%   just quit 
%   * Change adaptivity of continuation!
% 
%   * Higer-order Runge-Kutta, implicit -- done
%   * Lower  tolerance -- done
%   * Analytic jacobian into continuation


clear variables
close all
clc

K = 100; % stiffness (try 1, 1e1, 1e2, 1e3): K=1e3 is benchmark
too = 400;

y_prime = 0;

y0 = [2;0];

y_old = y0;

dt0 = 0.1;
t0 = 0;

T = t0;
DT = dt0;
Y = y0;
U1 = [y0];

[t_exact,y_exact] = ode15s(@rhs_ode,[0,too],y0);
% plot(t_exact,y_exact(:,1),'-or')
% %   
% keyboard

N = length(y0);
u1u2 = [y_old;y_old];


    hFig = figure;
    Cont_alpha_1(@rhs_cont,u1u2,t0,inf,...
                    'display','off',...
                    'monitorfcn',@MonitorFun,...
                    'direction',1,...
                    'step',0.1,...
                    'maxstep',inf,...
                    'adaptive',[.1,.1,.1],...
                    'tolfun',1e-8,...
                    'jacobian',[],...
                    'stopbtn',hFig);


   keyboard 


   ind = 1;
   
 u1u2 = [U1(:,end-ind); Y(:,end-ind)];
 t0 = T(end-ind);
 y_old = Y(:,end-ind);
 
 U1 = U1(:,1:end-ind);
 Y = Y(:,1:end-ind);
 T = T(1:end-ind); 
 DT = DT(1:end-ind);
 
      contHDFbeta2_1(@rhs_cont,u1u2,t0,inf,...
                    'monitorfcn',@MonitorFun,...
                    'direction',1,...
                    'step',.1,...
                    'maxstep',100,...
                    'adaptive',[1,.1,.1],...
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
    



    function [out,J] = rhs_cont(u1u2,t)       
        
        
        h = t-T(end);
        
        

        u1 = u1u2(1:N);
        u2 = u1u2(N+1:end);
        
        out1 = -u1+y_old+h*(5/12*rhs_ode(t+1/3*h,u1)...
                           -1/12*rhs_ode(t+h,u2));
                       
        out2 = -u2+y_old+h*(3/4*rhs_ode(t+1/3*h,u1)...
                           +1/4*rhs_ode(t+h,u2));

        
        out = [out1;out2];
        
        
        % jac w.r.t. y
        if nargout == 2, 
            
            J = eye(length(y))-dt*J_rhs_ode(t,y);
            
        end
    end

    
   
    

    function [MonVal,status] = MonitorFun(u1u2,t)
        
        u1 = u1u2(1:N);
        u2 = u1u2(N+1:end);        
        
        MonVal = [];
        status = [];
        
        if t>too, status = NaN; end
        
        h = t-T(end);
        
        y_old = u2;
        
T  = [T;t];
DT = [DT;h];
Y  = [Y y_old];
U1 = [U1 u1];
        
       cla
        plot(T,Y(1,:),'x-')
        hold on
        plot(t_exact,y_exact(:,1),'.-')
        title(['dt= ' num2str(h)])
        
        drawnow
% keyboard        
    end




end



