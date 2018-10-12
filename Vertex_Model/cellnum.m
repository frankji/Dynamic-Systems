function A =cellnum(cells)
    A=zeros(1,length(cells));
    for i = 1:length(cells)
        A(i)=length(cells{i});
    end
end