function [xi,yi] = returnXYIntersectionPoint(xt,yt,tht_to,xb,yb,tht_b)
    xi = (xt*tand(tht_to)-xb*tand(tht_b)+yb-yt)/(tand(tht_to)-tand(tht_b));
    yi = (tand(tht_to)*tand(tht_b)*(xt-xb)+tand(tht_to)*yb-tand(tht_b)*yt)/(tand(tht_to)-tand(tht_b));
end

