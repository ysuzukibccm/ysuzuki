t_ga_cpu=cputime;

ge = 1; % number of source points
num = 3*ge;
popsize = 30;

Ini=lhsdesign(popsize,num);

%range of source point coordinates and range of cylinder diameter
for i=1:popsize
    for j=1:3*ge
        if rem(j,3)~=0
            Ini(i,j)=6*Ini(i,j);
        elseif j == 3
            Ini(i,j)=Ini(i,j)*0.08;
        else
            Ini(i,j) = 2 * Ini(i,j) -1;
        end
    end
end

bpl=zeros(1,3*ge);
upl=zeros(1,3*ge);

for i=1:3*ge
    if rem(i,3)==1
        bpl(i)=0;
        upl(i)=6;
    elseif rem(i,3)==2
        bpl(i)=0;
        upl(i)=6;
    elseif i==3
        bpl(i)=0;
        upl(i)=0.08;
    else
        bpl(i)=-1;
        upl(i)=1;
    end
end

opts = gaoptimset('PlotFcns',{@gaplotbestf,@gaplotstopping}, 'Generations', 300 ,'PopulationSize', popsize, 'InitialPopulation', Ini,'mutationFcn',{@mutationadaptfeasible,0.3});
[S, val] = ga(@make_GAssample_n, num, [], [], [], [], bpl, upl,[], opts);

t_GA=cputime-t_ga_cpu;
fid = fopen( 'conclu_GA.txt', 'w+');
fprintf(fid,'CPU time consumed for conducting GA is %f. The optimization results are as follows. Tsai-Wu value = %f. The optimization paramters = (%f,%f,%f)',t_GA,val,S(1),S(2),S(3));
fclose(fid);
