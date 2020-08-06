function [t_fe,y_fe0,e_fe,r_fe,t_fe1,y_fe1]=fe_scheme(t_fe,y_fe0,eq,t0,tf,Ns,N)
    
    clear tall;
    for n=Ns:N
        h=2^(-n);
        tall=t0:h:tf;
        if Ns==1
            ysol=fsol(tall,0.1);
        else
            options=odeset('AbsTol',1e-6,'RelTol',1e-3);
            [t_ode45,ysol]=ode45(eq,tall,y_fe0,options);
        end
        tic
        for j=2:length(tall)+1
            t_fe(j)=t_fe(j-1)+h;
            f = eq(t_fe(j-1),y_fe0(:,j-1));
            y_fe0(:,j) = y_fe0(:,j-1)+h*f;
        end
        CPU_fe=toc;
        fprintf('Time = %i\n',CPU_fe);
        y_err = y_fe0(1,:);
        for i=1:length(tall)
            yerr(i)= ysol(i)-y_err(i);
        end
        e_fe(n+1)= max(norm(yerr));
        if n>=2
            r_fe(n+1) = e_fe(n+1)/e_fe(n);
        end
        if n==Ns
            t_fe1=t_fe(:);
            y_fe1=y_fe0(:,:);
        end
    end
    return
end
