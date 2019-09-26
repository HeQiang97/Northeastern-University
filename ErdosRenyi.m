
%ER
function [ net,sp, A ] = ErdosRenyi( N,ErdosRenyiProbablity )

mat = rand(N,N);
index = find(mat < ErdosRenyiProbablity);
A = zeros(N,N);
A(index) = 1;
sp = sparse(A);
sp = max(sp,sp');
net = CreateMap(sp);

end

