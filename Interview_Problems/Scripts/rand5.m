function rand5
% get rand5 from rand7

NBOX = 5; % NBOX is number of boxes. NBOX = 5 gives rand5 from rand7

a = zeros(NBOX,1);
i = 1;
% counter = 0;

K = 1e4;
C = 5+1; % if C is divisible by NBOX, the probability of the last 
% used box (in statistics-generating) being NBOX doubles: 
% from first condition + from 2-nd condition

k = zeros(K,1);

p = 0;

% STATISTICS
for n =1:K
    
    % THE ALGORITHM:
    for counter = 1:C
        x = randi(7);

        if x<=NBOX,
            a(x) = a(x)+1;
            p = x; % output 

        else
            a(i) = a(i)+1;
            p = i; % output

        end

        % counter used to re-distribute the extra 6 and 7
        i = i+1;
        if i>NBOX, i=1; end

    end
    k(n)=p;



end
histogram(k)

end