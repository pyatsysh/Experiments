classdef Class_bumpy_funs<handle
%  Collection of good examples to try optimisation

% 
% properties
%     x
%     y
%     xx
%     yy
% end

methods
    
    
    
    
    function z = Peaks(obj,x,y)
        % matlab's peaks function
        
        
        z =  3*(1-x).^2.*exp(-(x.^2) - (y+1).^2) ...
           - 10*(x/5 - x.^3 - y.^5).*exp(-x.^2-y.^2) ...
           - 1/3*exp(-(x+1).^2 - y.^2);

    end



   function out = grad_Peaks(obj,x,y)
        
        Dx = (-(2./3)).*exp(-2.*x - x.^2 - (1 + y).^2).*((-exp(2.*y)).*(1 + x) + ...
                9.*exp(2.*x).*(1 - 2.*x.^2 + x.^3) + exp(1 + 2.*x + 2.*y).*(3 - 51.*x.^2 + ...
                30.*x.^4 + 30.*x.*y.^5));
        
            
        Dy =  (1./3).*exp(-2.*x - x.^2 - (1 + y).^2).*(2.*exp(2.*y).*y - 18.*exp(2.*x).*(-1 + x).^2.*(1 + y) - ...
                6.*exp(1 + 2.*x + 2.*y).*y.*(-2.*x + 10.*x.^3 + 5.*y.^3.*(-5 + 2.*y.^2)));
            
        out = [Dx;Dy];
    
    end

    
    function out = Ros(obj,x1,x2)
        % Rosenbrock's fun
        out = 100*(x2-x1.^2).^2+(1-x1).^2;
    end

    function out = grad_Ros(obj,x1,x2)
        % grad of Rosenbrock's fun
        % f = 100*(x2-x1^2)^2+(1-x1)^2
        
        out = [-400*(x2-x1^2).*x1-2*(1-x1);
               200*(x2-x1^2)];
    end
    
end

    
end





