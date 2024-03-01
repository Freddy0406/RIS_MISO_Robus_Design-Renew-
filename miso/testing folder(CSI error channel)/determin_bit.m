function[b1,b2,b3,b4] = determin_bit(temp,d)
%To determin 4 bits

if temp(1,1) == -3*d && temp(1,2) == 3*d		  	   
    b1=1;
    b2=0;
    b3=1;
    b4=1;
    elseif temp(1,1) == -3*d && temp(1,2) == d		  	    
    b1=1;
    b2=0;
    b3=1;
    b4=0;            
    elseif temp(1,1) == -3*d && temp(1,2) == -1*d	 	  	
    b1=1;
    b2=1;
    b3=1;
    b4=0;
    elseif temp(1,1) == -3*d && temp(1,2) == -3*d	 	  	
    b1=1;
    b2=1;
    b3=1;
    b4=1;
    elseif temp(1,1) == -1*d && temp(1,2) == 3*d	 	  	
    b1=1;
    b2=0;
    b3=0;
    b4=1;
    elseif temp(1,1) == -1*d && temp(1,2) == 1*d	 	  	
    b1=1;
    b2=0;
    b3=0;
    b4=0;
    elseif temp(1,1) == -1*d && temp(1,2) == -1*d	 	  	
    b1=1;
    b2=1;
    b3=0;
    b4=0;
    elseif temp(1,1) == -1*d && temp(1,2) == -3*d	 	  	
    b1=1;
    b2=1;
    b3=0;
    b4=1;
    elseif temp(1,1) == 1*d && temp(1,2) == 3*d	 	  	
    b1=0;
    b2=0;
    b3=0;
    b4=1;
    elseif temp(1,1) == 1*d && temp(1,2) == 1*d	 	  	
    b1=0;
    b2=0;
    b3=0;
    b4=0; 
    elseif temp(1,1) == 1*d && temp(1,2) == -1*d	 	  	
    b1=0;
    b2=1;
    b3=0;
    b4=0;
    elseif temp(1,1) == 1*d && temp(1,2) == -3*d	 	  	
    b1=0;
    b2=1;
    b3=0;
    b4=1;
    elseif temp(1,1) == 3*d && temp(1,2) == 3*d	 	  	
    b1=0;
    b2=0;
    b3=1;
    b4=1;
    elseif temp(1,1) == 3*d && temp(1,2) == 1*d	 	  	
    b1=0;
    b2=0;
    b3=1;
    b4=0;
    elseif temp(1,1) == 3*d && temp(1,2) == -1*d	 	  	
    b1=0;
    b2=1;
    b3=1;
    b4=0;
    else			           	 
    b1=0;
    b2=1;
    b3=1;
    b4=1;             
end

end