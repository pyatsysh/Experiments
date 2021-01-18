function Test_Gaussian_RBF_1
% testing analytic expressions for the derivatives of 1D Gaussian RBF w.r.t. its params
% 2D0: 
%       * vectorise over cwm

close all
clear variables
clc


N = 100;
L = 8;
[x,dx,Dx,s] = collocati(N,'saus',-L/2+.5, L/2-.5,9);


c = 0;
w = 1;
m = 2;


    % CHECK DERVIATIVE W.R.T. m
    [dm,f] = NumJac1(@(m)Grbf1([c,w,m],x), m);

    % CHECK DERVIATIVE W.R.T. c
    [dc,f] = NumJac1(@(c)Grbf1([c,w,m],x), c); 

    % CHECK DERVIATIVE W.R.T. w
    [dw,f] = NumJac1(@(w)Grbf1([c,w,m],x), w);
    
    % CHECK DERVIATIVE W.R.T. cwm
    [Jn,f] = NumJac1(@(cwm)Grbf1(cwm,x), [c,w,m]);

    J = Grbf1_jac([c,w,m],x);     
    
figure
    nexttile
        plot(x,Grbf1([c,w,m],x))
        title('fun')

        
    nexttile
        plot(x,Grbf1_dc([c,w,m],x))
        hold on
        plot(x,dc,'.')
        title('c')

    nexttile
        plot(x,Grbf1_dw([c,w,m],x))
        hold on
        plot(x,dw,'.')
        title('w')    
        
    nexttile
        plot(x,Grbf1_dm([c,w,m],x))
        hold on
        plot(x,dm,'.')
        title('m')
       
    nexttile
        plot(x,J(:,1),'.',x,Jn(:,1))
        hold on
        plot(x,J(:,2),'.',x,Jn(:,2))
        plot(x,J(:,3),'.',x,Jn(:,3))
        plot(x,J)
        hold on
        plot(x,Jn,'.')
        
        
norm(J-Jn)       


keyboard

   
figure
plot(x,J(:,1))
hold on

    function out = Grbf1(cwm,x)
        % input: center, width, magnitude
        % x is column, {a,b,c} are scalars
        
        c = cwm(1);
        w = cwm(2);
        m = cwm(3);
        
        
        out = m*exp(-(x-c).^2/2/w^2);
    end

    
    function out = Grbf1_dc(cwm,x)
        % input: center, width, magnitude
        % x is column, {a,b,c} are scalars
        
        c = cwm(1);
        w = cwm(2);
        m = cwm(3);

        
        out = (m.*(-c + x))./(exp((-c + x).^2./(2.*w.^2)).*w.^2);
    end

    function out = Grbf1_dm(cwm,x)
        % input: center, width, magnitude
        % x is column, {a,b,c} are scalars
        
        
        c = cwm(1);
        w = cwm(2);
        m = cwm(3);

        
        out = exp(-(x-c).^2/2/w^2);
    end


    
    function out = Grbf1_dw(cwm,x)
        % input: center, width, magnitude
        % x is column, {a,b,c} are scalars
        
        c = cwm(1);
        w = cwm(2);
        m = cwm(3);

        out = (m.*(-c + x).^2)./(exp((-c + x).^2./(2.*w.^2)).*w.^3);
    end
    

    function out = Grbf1_jac(cwm,x)
        % jacobian with respect to cwm
        
        c = cwm(1);
        w = cwm(2);
        m = cwm(3);
        
        out =  [(m.*(-c + x))./(exp((-c + x).^2./(2.*w.^2)).*w.^2), ...
                (m.*(-c + x).^2)./(exp((-c + x).^2./(2.*w.^2)).*w.^3),...
                exp(-(x-c).^2/2/w^2)];
        
    end
end