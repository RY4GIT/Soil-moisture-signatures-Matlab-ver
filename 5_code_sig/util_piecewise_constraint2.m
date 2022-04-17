function [c,ceq] = util_piecewise_constraint2(P2)
    c = [];
    ceq(1) = P2(1) + P2(2)*P2(3) - P2(6);
    ceq(2) = P2(1) + P2(2)* (P2(3)+P2(4)) - P2(5) ; 
end
