function [new_v] = opt_v(Q,q,v)
    epsilon_mm = power(10,-8);
    D = eig(Q);                 %Eigenvalue
    lambda_Q = max(D,[],"all"); %Find maximum of Eigenvalue
    while true
        u = (Q-lambda_Q*eye(length(Q)))*v-q;
        
        next_v = -exp(1i*angle(u));

        prev = real(v'*u);
        after = real(next_v'*u);
        
        if(abs(prev-after)<epsilon_mm)
            break;
        else
        end
        v = next_v;
    end   
    new_v = next_v;
end