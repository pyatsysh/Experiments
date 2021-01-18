function Test_Gaussian_RBF_2
% Analytic jac of a sum of gaussian rbfs
% c - center
% w - width
% m - magnitude
% CWM = [c(:) w(:) m(:)] = [centers, widths, magnitudes] as columns
% Q - flattened CWM: [c(1);w(1);m(1);c(2);w(2);m(2),..]
% 2D0: 
%       * add normalisation analytically

close all
clear variables
clc


N = 100;
L = 8;
[x,dx,Dx,s] = collocati(N,'saus',-L/2+.5, L/2-.5,9);


c1 = -2;
w1 = 1;
m1 = 2;

c2 = 0;
w2 = 1.5;
m2 = 1;

c3 = 2.5;
w3 = .5;
m3 = 3;

CWM = [c1 w1 m1;
       c2 w2 m2;
       c3 w3 m3];
   
   
Q = CWM2Q(CWM);
norm(Q2CWM(Q)-CWM)
    
CWM1 = [c1 w1 m1;
       c2 w2 m2];


Q1 = CWM2Q(CWM1);

Q2CWM(Q1)

keyboard
    % CHECK DERVIATIVE W.R.T. cwm
    [Jn1,f1] = NumJac1(@(cwm)Grbf1(cwm,x), [c1,w1,m1]);
    [Jn2,f2] = NumJac1(@(cwm)Grbf1(cwm,x), [c2,w2,m2]);
    [Jn3,f3] = NumJac1(@(cwm)Grbf1(cwm,x), [c3,w3,m3]);

    
J = GrbfN_jac(Q2CWM(Q),x);
    Jn = [Jn1 Jn2 Jn3];
    
    JJn =  NumJac1(@(Q)GrbfN(Q2CWM(Q),x), Q);
    
    
figure
    nexttile
        plot(x,f1+f2+f3)
        title('fun')

        
    
       
    nexttile
        plot(x,J(:,1),'.',x,Jn(:,1))
        hold on
        plot(x,J(:,2),'.',x,Jn(:,2))
        plot(x,J(:,3),'.',x,Jn(:,3))
        plot(x,J)
        hold on
        plot(x,Jn,'.')
        plot(x,JJn,'o')
        
        
norm(J-Jn)  
norm(J-JJn)
norm(JJn)


keyboard

  J = GrbfN_jac(Q2CWM(Q1),x);
JJn =  NumJac1(@(Q)GrbfN(Q2CWM(Q),x), Q1);
norm(J-JJn)

   
figure
plot(x,J(:,1))
hold on


end    

    function out = Q2CWM(Q)
        % convert flattened parameter list to matrix
        
          n = length(Q)/3; % numbe of RBFs
          out = reshape(Q,[],n)';
    
    end

    function Q = CWM2Q(CWM);
        % convert CWM matrix to flattened parameter list
        % Q = [c1,w1,m1,c3,w2,m2,...]
        
        
        Q = CWM';
        Q = Q(:);
        
    end

    function out = GrbfN(CWM,x)
        % sum of N 1D rbfs
        % input: CWM = [c w m] = [centers, widths, magnitudes] as columns
        % x is column, 
        
        
        n = size(CWM,1);
        out = 0;
        for i = 1:n
            out = out + Grbf1(CWM(i,:),x);
        end
    end

    function out = GrbfN_jac(CWM,x)
        % Jacobian: sum of N 1D rbfs
        % input: CWM = [c w m] = [centers, widths, magnitudes] as columns
        % x is column, 
        
        
        n = 3*size(CWM,1); % size of jacobian
        
        out = zeros(length(x),n);
        
        for i = 1:3:n
            out(:,[i,i+1,i+2]) = Grbf1_jac(CWM((i-1)/3+1,:),x);
        end
    end





    function out = Grbf1(cwm,x)
        % one gaussian rbf
        % input: center, width, magnitude
        % x is column, cwm is 1x3 vector
        
        c = cwm(1);
        w = cwm(2);
        m = cwm(3);
        
        
        out = m*exp(-(x-c).^2/2/w^2);
    end

    
    

    function out = Grbf1_jac(cwm,x)
        % jacobian with respect to cwm
        % x is column, cwm is 1x3 vector
        
        
        c = cwm(1);
        w = cwm(2);
        m = cwm(3);
        
        out =  [(m.*(-c + x))./(exp((-c + x).^2./(2.*w.^2)).*w.^2), ...
                (m.*(-c + x).^2)./(exp((-c + x).^2./(2.*w.^2)).*w.^3),...
                exp(-(x-c).^2/2/w^2)];
        
    end
