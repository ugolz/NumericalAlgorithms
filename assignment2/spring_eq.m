function ydot = spring_eq(t,y)
    k=1;
    m=0.25;
    ydot=y*0;
    ydot(1)=y(2);
    ydot(2)=-(k/m)*y(1);
    return
end

