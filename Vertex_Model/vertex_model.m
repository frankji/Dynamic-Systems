%% Minimizing the Vertex Model
% Plot the location
locs=csvread('random_network.csv');
cells0=csvread('random_network_cells.csv');
cells=num2cell(cells0,2);
L0=[5.58,4.3];
drawpolys(locs,cells,L0);
title('Initial Network');
saveas(gcf,'plot1.eps','epsc');

% Calculate the energy
[E, E_elastic, E_tension, E_contractility]=vertexenergy(locs,cells,L0,1,1,1,1);
disp('Calculate the energy:');
disp('For all parameters equal 1, we only need to check the values');
disp(['E_elastic = ',num2str(E_elastic)]);
disp(['E_tension = ',num2str(E_tension)]);
disp(['E_contractility = ',num2str(E_contractility)]);
disp('These values are consistent with those in the reference');

% Calculate the forces on the network
epsilon=0.000001;
[E,~]=vertexenergy(locs,cells,L0,1,1,1,1);
for i = 1:1
    [Fs,~]=vertexforces(locs,cells,L0,1,1,1,1);
    locs1=locs+epsilon.*Fs;
    [E1,~]=vertexenergy(locs1,cells,L0,1,1,1,1);
    dE=(E-E1)./...
    (epsilon.*sum(sum(Fs.^2)));
end

disp(['dE = ', num2str(dE),' with epsilon = 0.000001']);


locs1=locs;
epsilon=0.05;
thres=100;
K=1;
A0=1;
i=0;
E_III=zeros(1,100000);
while thres > 10^-6
    i=i+1;
    [E, ~]=vertexenergy(locs1,cells,L0,1,-0.85,0.1,1);
    E_III(i)=E;
    [Fs,~]=vertexforces(locs1,cells,L0,1,-0.85,0.1,1);
    locs1=locs1+epsilon.*Fs;
    thres=sqrt(sum((sum(Fs.^2,2))));
%    [E1,~]=vertexenergy(locs1,cells,L0,1,-0.85,0.1,1);
%    dE=(E-E1)./...
%    (epsilon.*sum(sum(fs.^2)));
end
E_III=E_III(1:i);
figure;
drawpolys(locs1,cells,L0);
title('Case III');
saveas(gcf,'plot5_a.eps','epsc');
figure;
index=unique(round(logspace(1,log10(length(E_III)),100)));
loglog(index.*0.05,E_III(index)./(length(cells)*K*A0^2));
xlim([min(index.*0.05),max(index.*0.05)])
title('Case III');
ylabel('$$\frac{E(t)}{N_cKA_0^2}$$','Interpreter','latex')
xlabel('Time')
saveas(gcf,'plot5_b.eps','epsc');


locs1=locs;
epsilon=0.05;
thres=100;
K=1;
A0=1;
i=0;
E_I=zeros(1,100000);
while thres > 10^-6
    i=i+1;
    [E, ~]=vertexenergy(locs1,cells,L0,1,0.12,0.04,1);
    E_I(i)=E;
    [Fs,~]=vertexforces(locs1,cells,L0,1,0.12,0.04,1);
    locs1=locs1+epsilon.*Fs;
    thres=sqrt(sum((sum(Fs.^2,2))));
%    [E1,~]=vertexenergy(locs1,cells,L0,1,-0.85,0.1,1);
%    dE=(E-E1)./...
%    (epsilon.*sum(sum(fs.^2)));
end

[E, E_elastic, E_tension, E_contractility]=vertexenergy(locs1,...
    cells,L0,1,0.12,0.04,1);

disp('Calculate the energy:');
disp(['E_elastic/(Nc*K*A0^2) = ',...
    num2str(E_elastic/(length(cells)*K*A0^2))]);
disp(['E_tension/(Nc*K*A0^2) = ',...
    num2str(E_tension/(length(cells)*K*A0^2))]);
disp(['E_contractility/(Nc*K*A0^2) = ',...
    num2str(E_contractility/(length(cells)*K*A0^2))]);
disp(['E/(Nc*K*A0^2) = ',...
    num2str(E/(length(cells)*K*A0^2))]);
disp('These values are consistent with those in the reference');

E_I=E_I(1:i);
figure;
drawpolys(locs1,cells,L0);
title('Case I');
saveas(gcf,'plot6_a.eps','epsc');

figure;
index=unique(round(logspace(1,log10(length(E_I)),100)));
loglog(index.*0.05,E_I(index)./(length(cells)*K*A0^2));
xlim([min(index.*0.05),max(index.*0.05)])
title('Case I');
ylabel('$$\frac{E(t)}{N_cKA_0^2}$$','Interpreter','latex')
xlabel('Time')
saveas(gcf,'plot6_b.eps','epsc');



epsilon_L=0.00001;
epsilon=0.00001;
[E, E_elastic, E_tension, E_contractility, E_h]=scaleenergy(locs,...
    cells,L0,1,1,1,1,1);
[Fs, F_elastic, F_tension, F_contractility]=scaleforces(locs,...
    cells,L0,1,1,1,1,1);
locs1=(1-epsilon_L*E_h).*(locs+epsilon.*Fs);
L1=(1-epsilon_L*E_h).*L0;
[E1, ~]=scaleenergy(locs1,...
    cells,L1,1,1,1,1,1);
dE=(E-E1)./...
    (epsilon.*sum(sum(Fs.^2))+epsilon_L*E_h^2);
disp(['dE = ', num2str(dE)])

epsilon_L=0.001;
epsilon=0.05;
%Verify energy and forces with scaling coefficient
%Case I
locs1=locs;
thres=1;
L1=L0;
while thres>10^-4
    [E, E_elastic, E_tension, E_contractility, E_h]=scaleenergy(locs1,...
        cells,L1,1,1,0.12,0.04,1);
    [Fs, F_elastic, F_tension, F_contractility]=scaleforces(locs1,...
       cells,L1,1,1,0.12,0.04,1);
    locs1=(1-epsilon_L*E_h).*(locs1+epsilon.*Fs);
    L1=(1-epsilon_L*E_h).*L1;
    [E1, ~]=scaleenergy(locs1,...
      cells,L1,1,1,0.12,0.04,1);
     thres=(sum(sum(Fs.^2))+E_h.^2)./Nc;
%dE=(E-E1)./...
    %(epsilon.*sum(sum(Fs.^2))+epsilon_L*E_h^2);
end
figure;
drawpolys(locs1,cells,L1)
title('Case I with scaling coefficient')
saveas(gcf,'plot7_I.eps');

[E, ~]=vertexenergy(locs1,...
    cells,L1,1,0.12,0.04,1);

disp('Calculate the energy:');
disp(['E/(Nc*K*A0^2) = ',...
    num2str(E/(length(cells)*K*A0^2))]);
disp('This value is consistent with that for case I in the reference');

%Case II
locs1=locs;
thres=1;
L1=L0;
while thres>10^-4
    [E, E_elastic, E_tension, E_contractility, E_h]=scaleenergy(locs1,...
        cells,L1,1,1,0,0.1,1);
    [Fs, F_elastic, F_tension, F_contractility]=scaleforces(locs1,...
       cells,L1,1,1,0,0.1,1);
    locs1=(1-epsilon_L*E_h).*(locs1+epsilon.*Fs);
    L1=(1-epsilon_L*E_h).*L1;
    [E1, ~]=scaleenergy(locs1,...
      cells,L1,1,1,0,0.1,1);
     thres=(sum(sum(Fs.^2))+E_h.^2)./Nc;
%dE=(E-E1)./...
    %(epsilon.*sum(sum(Fs.^2))+epsilon_L*E_h^2);
end
figure;
drawpolys(locs1,cells,L1)
title('Case II with scaling coefficient')
saveas(gcf,'plot7_II.eps');

[E, ~]=vertexenergy(locs1,...
    cells,L1,1,0,0.1,1);

disp('Calculate the energy:');
disp(['E/(Nc*K*A0^2) = ',...
    num2str(E/(length(cells)*K*A0^2))]);
disp('This value is consistent with that for case II in the reference');

%Case III
locs1=locs;
thres=1;
L1=L0;
while thres>10^-4
    [E, E_elastic, E_tension, E_contractility, E_h]=scaleenergy(locs1,...
        cells,L1,1,1,-0.85,0.1,1);
    [Fs, F_elastic, F_tension, F_contractility]=scaleforces(locs1,...
       cells,L1,1,1,-0.85,0.1,1);
    locs1=(1-epsilon_L*E_h).*(locs1+epsilon.*Fs);
    L1=(1-epsilon_L*E_h).*L1;
    [E1, ~]=scaleenergy(locs1,...
      cells,L1,1,1,-0.85,0.1,1);
     thres=(sum(sum(Fs.^2))+E_h.^2)./Nc;
%dE=(E-E1)./...
    %(epsilon.*sum(sum(Fs.^2))+epsilon_L*E_h^2);
end
figure;
drawpolys(locs1,cells,L1)
title('Case III with scaling coefficient')
saveas(gcf,'plot7_III.eps');



%% Full operation
epsilon=0.05;
epsilon_L=0.0005;
d=0.05;

%% Case I
K=1;
lambda=0.12;
gamma=0.04;
newcells=cells;
newlocs=locs;
A0=1;
L1=L0;

i=0;
while i <= 360
    i=i+1;
    state=0;
%Randomly pick cell to split
    while state==0
        rand_cell=randsample(length(newcells),1);
        cell2pick=newcells{rand_cell};
        splits=randsample(length(cell2pick),2);
        split1=cell2pick(splits(1));
        split2=cell2pick(splits(2));

        [newcells, newlocs,state]=celldiv(newcells,...
            newlocs,L1, rand_cell,split1,split2);
    end
    thres=-1;
    tries=1;
    Nc=length(newcells);
    T1s=1;
    clf;
     drawpolys(newlocs,newcells,L1);
    drawnow;
    E_com=1000;
    while (thres < 0 || T1s+T2s > 0)
     [E, E_elastic, E_tension, E_contractility,...
         E_h]=scaleenergy(newlocs,newcells,L1,1,K,lambda,gamma,A0);
     [Fs, F_elastic, F_tension, F_contractility]=scaleforces(newlocs,...
      newcells,L1,1,K,lambda,gamma,A0);
       newlocs=(1-epsilon_L*E_h).*(newlocs+epsilon.*Fs);
       L1=(1-epsilon_L*E_h).*L1;
       [newcells, newlocs,T1s,T2s,tries] = transition(newcells, newlocs,...
          L1, [K,lambda,gamma], A0, d);
        E_cal=sum(sum(Fs.^2))./Nc+E_h.^2;
        thres=E_cal-E_com;
        E_com=E_cal;
    end

end

locs_I=newlocs;
cells_I=newcells;
L_I=L1;

figure;
drawpolys(locs_I,cells_I,L_I);
title('Case I');
saveas(gcf,'full_model_case_I.eps','epsc')

%% Case II
K=1;
lambda=0.0;
gamma=0.1;
newcells=cells;
newlocs=locs;
A0=1;
L1=L0;

i=0;
while i <= 360
    i=i+1;
    state=0;
%Randomly pick cell to split
    while state==0
        rand_cell=randsample(length(newcells),1);
        cell2pick=newcells{rand_cell};
        splits=randsample(length(cell2pick),2);
        split1=cell2pick(splits(1));
        split2=cell2pick(splits(2));

        [newcells, newlocs,state]=celldiv(newcells,...
            newlocs,L1, rand_cell,split1,split2);
    end
    thres=-1;
    tries=1;
    Nc=length(newcells);
    T1s=1;
    clf;
     drawpolys(newlocs,newcells,L1);
    drawnow;
    E_com=1000;
    while (thres < 0 || T1s+T2s > 0)
     [E, E_elastic, E_tension, E_contractility,...
         E_h]=scaleenergy(newlocs,newcells,L1,1,K,lambda,gamma,A0);
     [Fs, F_elastic, F_tension, F_contractility]=scaleforces(newlocs,...
      newcells,L1,1,K,lambda,gamma,A0);
       newlocs=(1-epsilon_L*E_h).*(newlocs+epsilon.*Fs);
       L1=(1-epsilon_L*E_h).*L1;
       [newcells, newlocs,T1s,T2s,tries] = transition(newcells, newlocs,...
          L1, [K,lambda,gamma], A0, d);
        E_cal=sum(sum(Fs.^2))./Nc+E_h.^2;
        thres=E_cal-E_com;
        E_com=E_cal;
    end

end

locs_II=newlocs;
cells_II=newcells;
L_II=L1;

figure;
drawpolys(locs_II,cells_II,L_II);
title('Case II');
saveas(gcf,'full_model_case_II.eps','epsc')

%% Case III
K=1;
lambda=-0.85;
gamma=0.1;
newcells=cells;
newlocs=locs;
A0=1;
L1=L0;
epsilon=0.05;

i=0;
while i <= 360
    i=i+1;
    state=0;
%Randomly pick cell to split
    while state==0
        rand_cell=randsample(length(newcells),1);
        cell2pick=newcells{rand_cell};
        splits=randsample(length(cell2pick),2);
        split1=cell2pick(splits(1));
        split2=cell2pick(splits(2));

        [newcells, newlocs,state]=celldiv(newcells,...
            newlocs,L1, rand_cell,split1,split2);
    end
    thres=-1;
    tries=1;
    Nc=length(newcells);
    T1s=1;
    E_com=1000;
    clf;
    drawpolys(newlocs,newcells,L1);
    drawnow;
    while (thres < 0|| T1s+T2s > 0)
     [E, E_elastic, E_tension, E_contractility,...
         E_h]=scaleenergy(newlocs,newcells,L1,1,K,lambda,gamma,A0);
     [Fs, F_elastic, F_tension, F_contractility]=scaleforces(newlocs,...
      newcells,L1,1,K,lambda,gamma,A0);
       newlocs=(1-epsilon_L*E_h).*(newlocs+epsilon.*Fs);
       L1=(1-epsilon_L*E_h).*L1;
       [newcells, newlocs,T1s,T2s,tries] = transition(newcells, newlocs,...
          L1, [K,lambda,gamma], A0, d);
        E_cal=sum(sum(Fs.^2))./Nc+E_h.^2;
        thres=E_cal-E_com;
        E_com=E_cal;
    end
end

locs_III=newlocs;
cells_III=newcells;
L_III=L1;

figure;
drawpolys(locs_III,cells_III,L_III);
title('Case III');
saveas(gcf,'full_model_case_III.eps','epsc')


%% Calculate the number of polygons
[N_III, edges_III]=histcounts(cellnum(cells_III),...
    3:max(cellnum(cells_III))+1);
N_III=N_III./length(cells_III);
figure;
exp_polys = [1, 6.78, 34.61, 38.28, 14.28, 2.17, 0.06] / 100;
exp_polys_err = [0.77, 4.176, 4.06, 6.29, 3.36, 1.76, 0.24] / 100;
plot_III=[exp_polys;N_III(1:7)]';
bar(plot_III);
hold on;
errorbar([1:7]-0.15, exp_polys, exp_polys_err, ...
'.', 'Color', [0.3,0.7,0.3], 'Linewidth', 3);
legend('Exp','III');
title('Case III');
xlabel('#sides');
set(gca,'XTickLabel',{'3','4','5','6','7','8','9'});
hold off;
saveas(gcf,'hist_case_III.eps','epsc');



[N_II, edges_II]=histcounts(cellnum(cells_II),...
    3:max(cellnum(cells_II))+1);
N_II=N_II./length(cells_II);
figure;
exp_polys = [1, 6.78, 34.61, 38.28, 14.28, 2.17, 0.06] / 100;
exp_polys_err = [0.77, 4.176, 4.06, 6.29, 3.36, 1.76, 0.24] / 100;
plot_II=[exp_polys;N_II(1:7)]';
bar(plot_II);
hold on;
errorbar([1:7]-0.15, exp_polys, exp_polys_err, ...
'.', 'Color', [0.3,0.7,0.3], 'Linewidth', 3);
legend('Exp','II');
title('Case II');
xlabel('#sides');
set(gca,'XTickLabel',{'3','4','5','6','7','8','9'});
hold off;
saveas(gcf,'hist_case_II.eps','epsc');


[N_I, edges_I]=histcounts(cellnum(cells_I),...
    3:max(cellnum(cells_I))+1);
N_I=N_I./length(cells_I);
figure;
exp_polys = [1, 6.78, 34.61, 38.28, 14.28, 2.17, 0.06] / 100;
exp_polys_err = [0.77, 4.176, 4.06, 6.29, 3.36, 1.76, 0.24] / 100;
plot_I=[exp_polys;N_I(1:7)]';
bar(plot_I);
hold on;
errorbar([1:7]-0.15, exp_polys, exp_polys_err, ...
'.', 'Color', [0.3,0.7,0.3], 'Linewidth', 3);
legend('Exp','I');
title('Case I');
xlabel('#sides');
set(gca,'XTickLabel',{'3','4','5','6','7','8','9'});
hold off;
saveas(gcf,'hist_case_I.eps','epsc');


%% Calculate the relative area of polygons
[An_I,A_I,C_I]=cellarea(cells_I,locs_I,L_I);
figure;
plot(C_I(C_I<9 & C_I>3), An_I(C_I<9 & C_I>3)./A_I,'o-', 'Color', ...
    [0.3,0,0.3], 'Linewidth', 2);
hold on;
exp_areas = [0.56, 0.82, 1.08, 1.36, 1.52];
exp_areas_err = [0.02, 0.01, 0.01, 0.02, 0.05];
errorbar(4:8, exp_areas, exp_areas_err, ...
'o-', 'Color', [0.3,0.7,0.3], 'Linewidth', 2);
title('Case I')
xlabel('Polygon');
ylabel('$$\frac{<A_n>}{<A>}$$','Interpreter','latex');
legend({'Case I','Exp'},'Location','northwest');
hold off;
saveas(gcf,'area_case_I.eps','epsc');


[An_II,A_II,C_II]=cellarea(cells_II,locs_II,L_II);
figure;
plot(C_II(C_II<9 & C_II>3), An_II(C_II<9 & C_II>3)./A_II,'o-', 'Color', ...
    [0.3,0,0.3], 'Linewidth', 2);
hold on;
exp_areas = [0.56, 0.82, 1.08, 1.36, 1.52];
exp_areas_err = [0.02, 0.01, 0.01, 0.02, 0.05];
errorbar(4:8, exp_areas, exp_areas_err, ...
'o-', 'Color', [0.3,0.7,0.3], 'Linewidth', 2);
title('Case II')
xlabel('Polygon');
ylabel('$$\frac{<A_n>}{<A>}$$','Interpreter','latex');
legend({'Case II','Exp'},'Location','northwest');
hold off;
saveas(gcf,'area_case_II.eps','epsc');


[An_III,A_III,C_III]=cellarea(cells_III,locs_III,L_III);
figure;
plot(C_III(C_III<9 & C_III>3), An_III(C_III<9 & C_III>3)./A_III,...
    'o-', 'Color', ...
    [0.3,0,0.3], 'Linewidth', 2);
hold on;
exp_areas = [0.56, 0.82, 1.08, 1.36, 1.52];
exp_areas_err = [0.02, 0.01, 0.01, 0.02, 0.05];
errorbar(4:8, exp_areas, exp_areas_err, ...
'o-', 'Color', [0.3,0.7,0.3], 'Linewidth', 2);
title('Case III')
xlabel('Polygon');
ylabel('$$\frac{<A_n>}{<A>}$$','Interpreter','latex');
legend({'Case III','Exp'},'Location','northwest');
hold off;
saveas(gcf,'area_case_III.eps','epsc');


%% Shape index
index_I=shapeindex(cells_I,locs_I,L_I);
index_II=shapeindex(cells_II,locs_II,L_II);
index_III=shapeindex(cells_III,locs_III,L_III);

disp(['p0 of Case I = ', num2str(index_I)]);
disp(['p0 of Case II = ', num2str(index_II)]);
disp(['p0 of Case III = ', num2str(index_III)]);

disp(['Based on the calculated value,',...
    'all these cases are more likely to be fluid, ',...
    'because they are greater than 3.81.'])
disp(['in these cases, Case I seems to be a little glassy, ',...
    'and Case II and Case III are more fluid-like.'])

disp(['Compared with the target shape index, ', ...
    'an increasing pattern of case shape index',...
    ' along with the increase of target shape index.']);
    