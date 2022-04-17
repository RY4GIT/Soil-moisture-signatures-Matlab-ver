function [c,ceq] = util_piecewise_constraint(P)
    c = [];
    ceq(1) = P(1) + P(2)*P(3) - P(5);
    ceq(2) = P(1) + P(2)* (P(3)+P(4)) - P(6) ; 
end
