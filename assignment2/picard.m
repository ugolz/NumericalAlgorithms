function [x,iter]=picard(g,x,tol,maxIter)
    err=tol+1;
    iter = 0;
    %error=zeros(maxIter, 1);
    while iter<=1 | err>tol & iter<=maxIter
        xs=g(x);
        fx=xs(1,:);
        err = abs(fx-x(1,:)); %test on the residual
        iter=iter+1;
        x=xs;
    end
    return
end