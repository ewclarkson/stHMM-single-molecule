function z = logsumexp2(x, y)
    % 
    % Define the function "logsum" defined by f(x,y)=log(exp(x)+exp(y))
    % 
    
    if x == -Inf && y == -Inf
       z = -Inf;
    else
       if x>y
         z = x+log(1+exp(y-x));
       else
         z = y+log(1+exp(x-y));
       end
    end
end