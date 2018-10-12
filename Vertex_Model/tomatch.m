function A=tomatch(X,Y)
    s=ismember(Y,X);
    s0=find(s==0);
    s1=find(s==1);
    A=0;
    if all(s0(2:end)-s0(1:end-1)==1) || all(s1(2:end)-s1(1:end-1)==1)
        A=1;
    end
end