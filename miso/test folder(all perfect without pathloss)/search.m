function [newlambda] = search(c,alpha,A,P0,M)
lambda_max = 100;
lambda_min = 0;
lambda = lambda_max;
while(1)
    w = pinv(power(abs(c),2)*A+lambda*eye(M))*(alpha*conj(c));
    if(abs(power( norm(w),2 )-P0)<10^-8)
        break;
    end
    if((power( norm(w),2 )-P0)<0)                 %lambda needs smaller
        lambda_max = lambda;
        lambda = lambda_min+(lambda_max-lambda_min)/2;
    elseif ((power( norm(w),2 )-P0)>0)            %lambda needs larger
        lambda_min = lambda;
        lambda = lambda_min+(lambda_max-lambda_min)/2;
    end

end
    newlambda = lambda;
end