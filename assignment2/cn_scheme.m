function [t_cn,y_cn0,e_cn,r_cn,t_cn1,y_cn1]=cn_scheme(t_cn,y_cn0,eq,t0,tf,Ns,N,tol) 
    for n=Ns:N
        h=2^(-n);
        tall=t0:h:tf;
        if Ns==1
            ysol=fsol(tall,0.1);
        else
            options=odeset('AbsTol',1e-6,'RelTol',1e-3);
            [t_ode45,ysol]=ode45(eq,tall,y_cn0,options);
        end
        iter =0;
        tic
        for j=2:length(tall)+1
            t_cn(j)=t_cn(j-1)+h;
            G=@(x) y_cn0(:,j-1)+0.5*h*(eq(t_cn(j),x)+eq(t_cn(j-1),y_cn0(:,j-1)));
            [y_cn0(:,j),iter_be] = picard(G,y_cn0(:,j-1),tol,100);
            iter = iter+iter_be;
        end
        time = toc;
        iter = iter/length(tall);
        
        fprintf('Time = %i: TS %i: %i \n',time,length(tall),iter);
        
        y_err = y_cn0(1,:);
        for i=1:length(tall)
            yerr(i)= ysol(i)-y_err(i);
        end
        e_cn(n+1)= max(norm(yerr));
        if n>=2
            r_cn(n+1) = e_cn(n+1)/e_cn(n);
        end
        if n==Ns
            t_cn1=t_cn(:);
            y_cn1=y_cn0(:,:);
        end
    end
    return
end