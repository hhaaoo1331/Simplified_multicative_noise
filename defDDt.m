    function [D,Dt] = defDDt
        % defines finite difference operator D
        % and its transpose operator
        D = @(U) ForwardD(U);
        Dt = @(X,Y) Dive(X,Y);