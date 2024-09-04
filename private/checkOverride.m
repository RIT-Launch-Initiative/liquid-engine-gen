function [Ac] = checkOverride(Rc_override,Ac)
    if(Rc_override ~= 0)
        if(convlength(Rc_override,'in','m') < sqrt(Ac/pi))
            warning(sprintf(['Ideal chamber radius (%.2f [in]) greater than ' ...
                'overridden radius (%.2f [in]), please correct for optimal ' ...
                'performance.\n'],convlength(sqrt(Ac/pi),'m','in'),Rc_override))
            Rc = convlength(Rc_override,'in','m');
            Ac = pi*Rc^2;
        else
            Rc = convlength(Rc_override,'in','m');
            Ac = pi*Rc^2;
        end
    end
end
