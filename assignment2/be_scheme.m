function [t_be,y_be0,e_be,r_be,t_be1,y_be1]=be_scheme(t_be,y_be0,eq,t0,tf,Ns,N,tol) 
    for n=Ns:N
        h=2^(-n);
        tall=t0:h:tf;
        if Ns==1
            ysol=fsol(tall,0.1);
        else
            options=odeset('AbsTol',1e-6,'RelTol',1e-3);
            [t_ode45,ysol]=ode45(eq,tall,y_be0,options);
        end
        iter =0;
        tic
        for j=2:length(tall)+1
            t_be(j)=t_be(j-1)+h;
            G=@(x) y_be0(:,j-1)+h*eq(t_be(j),x);
            [y_be0(:,j),iter_be] = picard(G,y_be0(:,j-1),tol,100);
            iter = iter+iter_be;
        end
        time = toc;
        iter = iter/length(tall);
        
        fprintf('Time = %i: TS %i: %i \n',time,length(tall),iter);
        
        y_err = y_be0(1,:);
        for i=1:length(tall)
            yerr(i)= ysol(i)-y_err(i);
        end
        e_be(n+1)= max(norm(yerr));
        if n>=2
            r_be(n+1) = e_be(n+1)/e_be(n);
        end
        if n==Ns
            t_be1=t_be(:);
            y_be1=y_be0(:,:);
        end
    end
    return
end