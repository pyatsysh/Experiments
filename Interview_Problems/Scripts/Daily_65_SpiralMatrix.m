function Print_Matr_Matlab
%  Daily coding problem (?)
%  print matrix in a spiral clockwise, starting form element [1,1]
clc
clear all
close all


A = [1, 2, 3, 4;
     5, 6, 7, 8;
     9, 10, 11, 12];
 
Pop; 
    function Pop
        % print outer layer clockwise
        % popping layers as they are printed
        % call recursively
        
        
        fprintf('%d\n',A(1,1:end))
        A(1,:) = [];
        
        
        if size(A,1)>1, % we have rows
            fprintf('%d\n',A(1:end,end)); 
            A(:,end) = [];
        
        
        
            fprintf('%d\n',A(end,end:-1:1))
            A(end,:) = [];
            
            if size(A,2)>1, % we have columns
                fprintf('%d\n',A(end:-1:1,1))
                A(:,1)=[];
            else
                return
            end
        else
            return
        end
    
        Pop;
    end





end