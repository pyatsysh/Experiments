function FindSign_1
% Problem: you are given n integers X1..Xn and integer S. Put signes in
% front of each X1..n so that their sum equals S.

% Idea inspired by 
% http://algolist.manual.ru/ai/ga/dioph.php

% Let's call "a" the vector of signes and try to minimise the function
% Delta = {S-a'*x}. We will use a genetic approach, where the pairs of a's 
% to be coupled to produce the new a are sampled from the discrete probability 
% distro that favours those a's which minimise Delta.


% SUMMARY 
% parts: 
% i  - Distro fun -- generated using distance fun. Sample N pairs from it. 
%      Distance is our target fun to be minimized/maximised
% ii - Paring fun. This computes one sample from two by combining them. 
%      Used to get N samples. Note: must allow probability of non-paring
%      pass to next generation.
% ii - Mutating fun. 
%      If all N samples are identical, we MUST mutate, because we are
%      stuck. Note: mutate a bit on every generation, prevents getting
%      stuck and accelerates convergence!

% NOTE: 
% Theory is constructed by representing the evolving data structures as
% strings. The most important theorem is Schemas theorem, which proves that
% if we are truly random, and add mutations -- then convergence is
% exponential with size of population or something like that. 

clc
clear all 
close all

[a,X,S,n] = Tester(1500);
counter = 0;
stuck_counter = 0;

    % population size
    N = 5;

    del = 0;
    while nnz(del)<N % while protects against accidentally solving with random IG
        % 1. Generate a1: 1-st generation of N samples

        a1 = randi([0,1],n,N);
        a1(a1==0) = -1;

        % 2. Compute distances (values of target fun)
        del = (a1'*X-S).^2;
        fprintf('randomly solved\n')

    end
    
% iterate    
fprintf('mean distance: %f\n', mean(sqrt(del)))
while nnz(del)==N
  
counter = counter+1;    
    
    % 3. Compute distance-based distro -- create function to make optional
    % (higher values of distance correspond to lower probability)
    pn = 1./del/sum(1./del);

    % 4. Generate pairs of samples from distance-based distro. 
    %    bin -- contains indexa of pairs
    x = rand(2*N,1);
    % Y = discretize(x,[0;cumsum(pn)]);
    [freq, edges, bin] = histcounts(x,[0;cumsum(pn)]); 

    % 5. Pair 2*N samples to get new generation -- create function to make optional
    seps = randi(n,N,1); % use 1-st from 1:seps, 2-nd from seps+1:end

    a2 = a1(:,bin(1:2:end));
    for i=1:N
        % no cross-over case
        if seps(i)==n, continue; end 
        
        % cross-over between pair
        a2(seps(i)+1:end,i) = a1(seps(i)+1:end,bin(2*i));
        
        % mutation on every step (this speeds-up convergence. Try commenting out!)
        j = randi(n);
        a2(j,i)=-a2(j,i);
       
    end


    
    % MUTATE IF a1 HAS THE SAME COLUMS (I.E., WE ARE STUCK)
    if norm(diff(a2,1,2))==0,
        fprintf('Got Stuck. Mutating\n')
        stuck_counter = stuck_counter+1;
        i = randi(n);
        j = randi(N); % mutate j columns
        a2(i,1:j) = -a2(i,1:j);
        
        % this is not random! So may intro bias!
        % it also jumps up in distance greatly
%         a2(1:n+1:end)=-a2(1:n+1:end); 
    end
    
    
    % re-assign a1 and del for the loop
    a1=a2;
    del = (a1'*X-S).^2;
    
    fprintf('mean distance: %f\n', mean(sqrt(del))) 
    
    
    
    
    
    
end

fprintf('%d-dim prob in %d steps; %d times stuck \n',n, counter, stuck_counter)
% a1
% S
% [a1(:,find(del==0)) a X]
% [X'*a1(:,find(del==0)) S]
 keyboard
    function [a0,X,S,n] = Tester(N)
        % GENERATE TEST-CASE: a0, X, S and n<N

        a0 = randi([-1,1],N,1); 
        a0(a0==0) = [];

        X  = randi(10,size(a0));
        S  = a0'*X;
        n  = length(a0);
    end

end


