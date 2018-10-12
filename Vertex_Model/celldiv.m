function [newcells, newlocs, state]=celldiv(cells,...
    locs,L, cellix,split1,split2)
    ixs=cells{cellix};
    N=length(ixs);
    p1=find(ixs==split1);
    p3=find(ixs==split2);
    %State is used to check crossover
    state=1;
    if(p1 > p3)
        temp=p1;
        p1=p3;
        p3=temp;
    end
    p2=mod(p1+1-1,N)+1;
    p4=mod(p3+1-1,N)+1;
    newcells=cells;
    for i=[1:cellix-1,(cellix+1):length(newcells)]
        cell=newcells{i};
        pp1=find(cell==ixs(p1));
        pp2=find(cell==ixs(p2));
        tail=length(cell)-1;
        if sum(ismember([ixs(p1),ixs(p2)],cell))==2 
            newcells{i}=[cell(1:pp1),length(locs)+1,cell(pp1+1:end)];
            if pp2-pp1==-1
                newcells{i}=[cell(1:pp2),length(locs)+1,cell(pp2+1:end)];
            end
            if pp2-pp1==tail
                newcells{i}=[cell(1:pp2),length(locs)+1,cell(pp2+1:end)];
            end
        end
        cell=newcells{i};
        pp3=find(cell==ixs(p3));
        pp4=find(cell==ixs(p4));
        if sum(ismember([ixs(p3),ixs(p4)],cell))==2 
            newcells{i}=[cell(1:pp3),length(locs)+2,cell(pp3+1:end)];
            if pp4-pp3==-1 
               newcells{i}=[cell(1:pp4),length(locs)+2,cell(pp4+1:end)];
            end
            if pp4-pp3==tail
               newcells{i}=[cell(1:pp4),length(locs)+2,cell(pp4+1:end)];
            end
        end
    end
    newcells=[newcells(1:cellix-1);newcells((cellix+1):end)];
    r12=mod(locs(ixs(p2),:) - locs(ixs(p1),:) + L/2, L) - L/2;
    m12=mod(r12./2+locs(ixs(p1),:)+L/2,L)-L/2;
    m12(m12<0)=m12(m12<0)+L(m12<0);
    r34=mod(locs(ixs(p4),:) - locs(ixs(p3),:) + L/2, L) - L/2;
    m34=mod(r34./2+locs(ixs(p3),:)+L/2,L)-L/2;
    m34(m34<0)=m34(m34<0)+L(m34<0);
    m_m=mod(m12 - m34+ L/2, L) - L/2;
    newlocs=[locs;m12;m34];
    N=length(newlocs);
    all_sign=zeros(1,length(ixs));
    if p3 < length(ixs)
        circle=[ixs(1:p1),N-1,ixs(p2:p3),N,ixs(p4:end)];
    end
    if p3 == length(ixs)
        circle=[ixs(1:p1),N-1,ixs(p2:p3),N];
    end
    for i = 1:length(circle)
        m=circle(i);
        temp1=mod(-m34+newlocs(m,:)+L/2,L)-L/2;
        temp=cross([temp1,0],[m_m,0]);
        all_sign(i)=temp(3);
    end
    newcells{end+1}=circle(all_sign<=0);
    newcells{end+1}=circle(all_sign>=0);
    %The sequence of new cell should be tandem in original circle.
    if tomatch(circle(all_sign<=0),...
            circle)*tomatch(circle(all_sign>=0),circle)==0
        newcells=cells;
        newlocs=locs;
        state=0;
    end
    connection=zeros(1,length(newlocs));
    for i = 1:length(newcells)
        cell=newcells{i};
        for ix = cell
            connection(ix)=connection(ix)+1;
        end
    end
    if all(connection==3)==0 
        newcells=cells;
        newlocs=locs;
        state=0;
    end
end