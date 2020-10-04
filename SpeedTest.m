N = 2000;
a = rand(N);
b = rand(N);
c = zeros(N);
tic
for i = 1:N
    for j = 1:N
        c(i,j) = a(i,j) + b(i,j);
    end
end
toc;