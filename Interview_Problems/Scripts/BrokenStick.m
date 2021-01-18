% clc

N = 3e7;

t = [];

l1 = [];
l2 = [];
l3 = [];


x = rand(2,N);
MIN = min(x); 
MAX = max(x);


L3 = abs(x(2,:)-x(1,:));

% maximal length
T = max([MIN;1-MAX;L3]);

mean(T)

figure
histogram(T,'normalization','pdf')


q = linspace(1/3,1/2);
q1 = linspace(1/2,1);


hold on
plot(q,18*(q-1/3),'linewidth',2)
plot(q1,6*(1-q1),'linewidth',2)

z = reshape(MIN,3,1e7);
Z = max(z);

figure
histogram(Z,'normalization','pdf')
