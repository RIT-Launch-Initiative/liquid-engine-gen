function [Ae] = calcAe(Astar,Me,gam)
    
    term1 = Astar/Me;
    term2 = (2/(gam+1))*(1+((gam-1)/2)*Me^2);
    exp = (gam+1)/(2*(gam-1));

    Ae = term1*term2^exp;

end

