function [Astar] = calcAstar(mdot,Po,To,R,gam)

    term1 = mdot/Po;
    term2 = sqrt(To*R/gam);
    term3 = (gam+1)/2;
    exp = -(gam+1)/(2*(gam-1));

    Astar = term1*term2*1/(term3^exp);

end

