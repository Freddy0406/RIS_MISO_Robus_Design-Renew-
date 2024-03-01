function [new_v] = discetet_phase(bits,v)
    partition = power(2,bits);
    phase_array = 0:(2*pi/partition):2*pi;
    new_v = zeros(size(v));    
    for i = 1:length(v)
        for j = 1:length(phase_array)-1
            if(j==1)
                distance = abs(wrapTo2Pi(angle(v(i)))-phase_array(j));
                decision = j;
            else
                if(abs(wrapTo2Pi(angle(v(i)))-phase_array(j))<distance)
                    distance = abs(wrapTo2Pi(angle(v(i)))-phase_array(j));
                    decision = j;
                end
            end
        end
        new_v(i) = exp(1i*phase_array(decision));
        if(wrapTo2Pi(angle(new_v(i)))==2*pi)
            new_v(i) = exp(1i*phase_array(length(phase_array)-1));
        end
    end
end