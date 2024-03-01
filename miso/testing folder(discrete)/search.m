function [newlambda] = search(c,alpha,A,P0,M)
lambda_max = 100;
lambda_min = 0;
lambda = lambda_max;
while(1)
    w = pinv(power(abs(c),2)*A+lambda*eye(M))*(alpha*conj(c));
    if(round( power( norm(w),2 ) , 5)<P0)                 %lambda needs smaller
        lambda_max = lambda;
        lambda = lambda_min+(lambda_max-lambda_min)/2;
    elseif (round( power( norm(w),2 ) , 5)>P0)            %lambda needs larger
        lambda_min = lambda;
        lambda = lambda_min+(lambda_max-lambda_min)/2;
    elseif(round( power( norm(w),2 ) , 5) == P0)
        break;
    end
end
    newlambda = lambda;
end