function normal = normalization(xx,kind)

[m,n]  = size(xx);
normal = zeros(m,n);

if kind == 1  
    for i = 1:n
        ma=max(xx);
        mi=min(xx);
       normal(:,i) = ( xx(:,i)-mi )./( ma-mi );
    end
end
if kind == 2
    for i = 1:n
        mea = mean( xx(:,i) );
        va = var( xx(:,i) );
      normal(:,i) = ( xx(:,i)-mea )/va;
    end
end
