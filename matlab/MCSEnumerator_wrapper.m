function [mcs] = MCSEnumerator_wrapper(cbmodel, max_dels, substrate, product, min_yield, min_growth, notknockable)
% Helper function to run MCSEnumerator
% Modified from the original code by von Kamp and Klamt (2014) 
% to facilitate setting up the cobra model directly.

cnap = CNAcobra2cna(cbmodel);
cnap.reacMax(cnap.reacMax >= 1000) = inf;
cnap.reacMin(cnap.reacMin <= -1000) = -inf;

substrate_idx = find(strcmp(cbmodel.rxns, substrate));
product_idx = find(strcmp(cbmodel.rxns, product));
n = length(cbmodel.rxns);
T = zeros(1, n);
T(1, product_idx) = 1; %#ok<*FNDSB>
T(1, substrate_idx) = -min_yield;
t = 0;

% desired modes setup
biomass = find(cbmodel.c);
D= zeros(1, n);
D(1, biomass)= -1;
d= -min_growth;


preprocess = 1;

st=initsmat(cnap.stoichMat,cnap.mue,cnap.macroComposition,cnap.macroDefault,cnap.specInternal);
[~,q]=size(st);

if(preprocess)
    % FVA to find blocked reactions and to check feasibility of target flux vectors and of desired flux vectors
    disp(' ');
    check_via_fva=1;
    
    if(check_via_fva)
        disp('Flux Variability Analysis ...');
        cgp= Cplex();
        cgp.Param.emphasis.numerical.Cur= 1;
        cgp.Param.simplex.tolerances.optimality.Cur= cgp.Param.simplex.tolerances.optimality.Min;
        cgp.Model.A= st;
        cgp.Model.ub= cnap.reacMax;
        %cgp.Model.ub= inf(cnap.numr,1);
        cgp.Model.lb= cnap.reacMin;
        %cgp.Model.lb= -inf(cnap.numr,1);
        %cgp.Model.lb(stirr) = 0;
        cgp.Model.lhs= zeros(size(cgp.Model.A, 1), 1);
        cgp.Model.rhs= zeros(size(cgp.Model.A, 1), 1);
        cgp.Model.obj=zeros(q, 1);
        cgp.Model.sense= 'maximize';
        %res= cgp.solve()
        cgp.DisplayFunc=[];
        
        fva= {cgp.Model.lb, cgp.Model.ub};
        sense= {'minimize', 'maximize'};
        for i= 1:q
            %for i= 1:0
            cgp.Model.obj(:)= 0;
            cgp.Model.obj(i)= 1;
            for j= 1:2
                cgp.Model.sense= sense{j};
                x= cgp.solve();
                if x.status ~= 1 && x.status~=2
                    fprintf('status %d %s for %s, %d\n', x.status,x.statusstring, sense{j}, i);
                end
                if x.status == 1
                    fva{j}(i)= x.objval;
                end
            end
        end
        
        fvalb= fva{1};
        fvaub= fva{2};
        fva_tol= cgp.Param.simplex.tolerances.optimality.Min;
        fvalb(abs(fvalb) < fva_tol)= 0;
        fvaub(abs(fvaub) < fva_tol)= 0;
        same= fvalb == fvaub;
        blocked_fva= same & (fvalb == 0);
        disp(['Found and removed ',num2str(sum(blocked_fva)),' blocked reactions via FVA!']);
        disp(' ');
        
        %Check feasibility of target flux vectors
        cgp.addRows(-inf(size(T,1),1),T,t);
        x= cgp.solve();
        disp(' ');
        if x.status == 3
            disp('Target flux vectors infeasible!');
            fprintf('status %d %s \n', x.status,x.statusstring, sense{j}, i);
            disp('Exit due to infeasibility of target flux vectors - no cuts required to block target flux vectors!');
            return;
        else
            disp('Target flux vectors feasible!');
        end
        
    end
    
    
    extraConst=zeros(0,q);
    extraConstrhs=zeros(0,1);
    zerorow=zeros(1,q);
    for i=1:q
        if(cnap.reacMin(i)~=0 && cnap.reacMin(i)~=-inf)
            extraConst(end+1,:)=zerorow;
            extraConst(end,i)=-1;
            extraConstrhs(end+1,1)=-cnap.reacMin(i);
        end
        if(cnap.reacMax(i)~=inf)
            extraConst(end+1,:)=zerorow;
            extraConst(end,i)=1;
            extraConstrhs(end+1,1)=cnap.reacMax(i);
        end
    end
    
    if(numel(extraConst>0))
        disp(['Added ',num2str(size(extraConst,1)),' min/max constraints different from 0 or +/- inf!']);
        T=[T;extraConst];
        t=[t;extraConstrhs];
        if(~isempty(D))
            D=[D;extraConst];
            d=[d;extraConstrhs];
        end
    end
    
    disp([num2str(numel(t)),' constraints for target flux vectors.']);
    disp([num2str(numel(d)),' constraints for desired flux vectors.']);
    
    %%% reduction
    Tinvolved=find(any(abs(T)>cnap.epsilon,1));
    Dinvolved=find(any(abs(D)>cnap.epsilon,1));
    TDinvolved=unique([Tinvolved,Dinvolved]);
    notknockable2=unique([notknockable,TDinvolved]);
    
    %%Old version of network compression
    %removal of conservation relations
    %[zw, bc]= rref(st',cnap.epsilon);
    %disp(['Removed ',num2str(size(st,1)-numel(bc)),' metabolites from conservation relations.']);
    %st= st(bc, :);
    %
    %[rd, sub, irrev_rd, sub_irr_viol]= subsets_reduction(st, stirr,blocked_fva, notknockable2);
    
    %removal of conservation relations
    %[zw, bc]= rref(rd',cnap.epsilon);
    %disp(['Removed ',num2str(size(rd,1)-numel(bc)),' metabolites from conservation relations.']);
    %rd= rd(bc, :);
    
    %new compression
    %cnap1.stoichMat=st;
    %cnap1.reacMin=-ones(size(st,2),1);
    %cnap1.reacMin(find(stirr))=0;
    %cnap1.reacID=cnap.reacID;
    %cnap1.epsilon=cnap.epsilon;
    %cnap1=CNAgenerateMFNetwork(cnap1,1);
    %[rd, irrev_rd, sub]= CNAcompressMFNetwork(cnap1,notknockable2,[],1,0,1,blocked_fva,0);
    [rd, irrev_rd, sub]= CNAcompressMFNetwork(cnap,notknockable2,[],1,0,1,find(blocked_fva),0);
    sub=sub';  %Transposed version required here (compatible with metatool)
    
    disp(' ');
    disp(['Final size of reduced system: ',num2str(size(rd))]);
    disp(' ');
    
    q=size(rd,2);
    
    todel=[];
    Tidx=[];
    for i=1:numel(Tinvolved)
        zw=find(sub(:,Tinvolved(i)));
        if(isempty(zw))
            %disp(['Warning: reaction ',deblank(cnap.reacID(Tinvolved(i),:)),' contained in definition of target flux vectors is blocked.']);
            %disp(' ');
            todel=[todel, i];
        else
            Tidx=[Tidx zw];
        end
    end
    Tinvolved(todel)=[];
    Tsub=zeros(size(T,1),q);
    Tsub(:,Tidx)=T(:,Tinvolved);
    
    zw=all(Tsub==0,2);
    todel=[];
    for j=1:numel(zw)
        if(zw(j))
            if(t(j)==0)
                todel=[todel j];
                %disp(['Removed always fulfilled constraint for target flux vectors.']);
                %disp(' ');
            else
                disp('One constraint for target flux vectors is infeasible. Empty target set! Exit!');
                disp(' ');
                return;
            end
        end
    end
    Tsub(todel,:)=[];
    t(todel)=[];
    
    if(~isempty(D))
        todel=[];
        Didx=[];
        for i=1:numel(Dinvolved)
            zw=find(sub(:,Dinvolved(i)));
            if(isempty(zw))
                %disp(['Warning: reaction ',deblank(cnap.reacID(Dinvolved(i),:)),' contained in definition of desired flux vectors is blocked.']);
                %disp(' ');
                todel=[todel, i];
            else
                Didx=[Didx zw];
            end
        end
        Dinvolved(todel)=[];
        Dsub=zeros(size(D,1),q);
        Dsub(:,Didx)=D(:,Dinvolved);
        
        zw=all(Dsub==0,2);
        todel=[];
        for j=1:numel(zw)
            if(zw(j))
                if(d(j)==0)
                    todel=[todel j];
                    %disp(['Removed always fulfilled constraint for desired flux vectors.']);
                    %disp(' ');
                else
                    disp('One constraint for desired flux vectors is infeasible. Exit!');
                    disp(' ');
                    return;
                end
            end
        end
        Dsub(todel,:)=[];
        d(todel)=[];
    else
        Dsub=[];
    end
    
    disp([num2str(numel(t)),' remaining constraints for target flux vectors.']);
    disp([num2str(numel(d)),' remaining constraints for desired flux vectors.']);
    disp(' ');
    
    notknockablesub=[];
    numk=0;
    for i=1:numel(notknockable)
        zw=find(sub(:,notknockable(i)));
        if(numel(find(sub(zw,:)))==1)
            numk=numk+1;
            notknockablesub(numk)=zw; %#ok<*AGROW>
        end
    end
    
    %Check feasibility of desired flux vectors in reduced system
    if(~isempty(Dsub))
        disp('Checking feasibility of desired flux vectors ...');
        disp(' ');
        cgp=Cplex();
        cgp.Param.emphasis.numerical.Cur= 1;
        cgp.Param.simplex.tolerances.optimality.Cur= cgp.Param.simplex.tolerances.optimality.Min;
        cgp.Model.A= rd;
        cgp.Model.ub= Inf(size(cgp.Model.A, 2), 1);
        cgp.Model.lb= -Inf(size(cgp.Model.A, 2), 1);
        cgp.Model.lb(irrev_rd ~= 0)= 0;
        cgp.Model.lhs= zeros(size(cgp.Model.A, 1), 1);
        cgp.Model.rhs= zeros(size(cgp.Model.A, 1), 1);
        cgp.Model.obj=zeros(q, 1);
        cgp.Model.sense= 'maximize';
        cgp.DisplayFunc=[];
        %		cgp.addRows(-inf(size(Tsub,1),1),-Tsub,-t);  % new!
        cgp.addRows(-inf(size(Dsub,1),1),Dsub,d);
        
        x= cgp.solve();
        disp(' ');
        if x.status == 3
            disp('Desired flux vectors infeasible!');
            fprintf('status %d %s \n', x.status,x.statusstring);
            disp('Exit due to infeasibility of desired flux vectors!');
            return;
        else
            disp('Desired flux vectors feasible!');
        end
        disp(' ');
        
        %FVA to find essential reactions for desired flux vectors
        ess_desired=[];
        %fp=fopen('essreacs.txt','w');
        sense= {'minimize', 'maximize'};
        
        for i= 1:size(rd,2)
            ess=0;
            cgp.Model.obj(:)= 0;
            cgp.Model.obj(i)= 1;
            val=[NaN NaN];
            for j= 1:2
                cgp.Model.sense= sense{j};
                x= cgp.solve();
                if x.status ~= 1 && x.status~=2 && x.status~=4
                    fprintf('Error while testing esential reactions for desired flux vectors.');
                    fprintf(' Optimization status %d %s for %s, %d\n', x.status,x.statusstring, sense{j}, i);
                    return;
                end
                if(x.status==4)
                    disp(['Reaction ',num2str(i),': Warning: solver status is 4. Assume unbounded solution.']);
                elseif(x.status==1 && ((j==1 && x.objval>fva_tol) || (j==2 && x.objval<-fva_tol)))
                    ess=j;
                    val(j)=x.objval;
                end
            end
            if(ess)
                %zw=find(sub(:,i));
                %fprintf(fp,[cnap.reacID(zw,:),' ',num2str(ess),' ',num2str(val(1)),' ',num2str(val(2)),'\n']);
                ess_desired=[ess_desired i];
            end
        end
        
        %fclose(fp);
        disp(['Found ',num2str(numel(unique(ess_desired))),' essential reactions in reduced system for desired flux vectors!']);
        %cnap.reacID(ess_desired,:)
        disp(' ');
        notknockablesub=unique([notknockablesub ess_desired]);
    end
    
    disp(['Final number of reactions that cannot be knocked-out: ',num2str(numel(notknockablesub))]);
    disp(['Final number of reactions that can be knocked-out: ',num2str(q-numel(notknockablesub))]);
    disp(' ');
    disp('Preprocessing finished!');
    disp(' ');
    
else %without preprocessing; the T, D etc. must be prepared accordingly; reacMin and reacMax are not considered separately
    Dsub=D;
    Tsub=T;
    rd=initsmat(cnap.stoichMat,cnap.mue,cnap.macroComposition,cnap.macroDefault,cnap.specInternal);
    irrev_rd=cnap.reacMin>=0;
    q=size(rd,2);
    notknockablesub=notknockable;
    disp(' ');
    disp(['Size of network: ',num2str(size(rd))]);
    disp(['Size of desired system: ',num2str(size(Dsub))]);
    disp(['Size of target system: ',num2str(size(Tsub))]);
    disp(' ');
    
end


%% MCS computation

if (q-numel(notknockablesub))==0
    mcs=[];
    disp('No reaction can be knocked-out - no cut set can exist!');
    return;
end


cuts= true(1, q);
cuts(notknockablesub)=false;

options.cuts= cuts;
options.method= 5;
options.workmem= 45056;

clearvars -except max_dels rd irrev_rd Tsub t options Dsub d q sub

data.max_dels = max_dels;
data.rd = rd;
data.irrev_rd = irrev_rd;
data.Tsub = Tsub;
data.t = t;
data.Dsub = Dsub;
data.d = d;
data.sub = sub;
data.options = options;


[mcsn, ~]= k_shortest_mcs(inf, max_dels, rd, irrev_rd, [], Tsub, t, [], 1, options);


%Preparing tests for constrained cut sets
if(~isempty(Dsub))
    cgp=Cplex();
    cgp.Param.emphasis.numerical.Cur= 1;
    cgp.Param.simplex.tolerances.optimality.Cur= cgp.Param.simplex.tolerances.optimality.Min;
    cgp.DisplayFunc=[];
    q= size(rd, 2);
    cgp.Model.A= rd;
    cgp.Model.ub= Inf(q, 1);
    cgp.Model.lb= -Inf(q, 1);
    cgp.Model.lb(irrev_rd ~= 0)= 0;
    cgp.Model.lhs= zeros(size(cgp.Model.A, 1), 1);
    cgp.Model.rhs= zeros(size(cgp.Model.A, 1), 1);
    
    cgp.Model.A(end+1:end+size(Dsub,1),:)=Dsub;
    cgp.Model.lhs(end+1:end+size(Dsub,1))=-Inf;
    cgp.Model.rhs(end+1:end+size(Dsub,1))=d;
    
    cgp.Model.obj=zeros(q, 1);
    cgp.Model.sense= 'maximize';
end



%[mcsn, ~]= k_shortest_mcs(inf, max_dels, rd, irrev_rd, [], Tsub, t, [], 1, options);

disp(' ');
disp([num2str(size(mcsn,2)),' (compressed) MCSs found.']);

%% LP check to filter cMCS from MCS

if(~isempty(mcsn))
    if(~isempty(Dsub))
        disp(' ');
        disp('Testing which MCSs are constrained MCSs ...');
        ct= cputime;
        [feasible, ~]= lp_check_cmcs(mcsn, cgp);
        disp(cputime - ct);
        
        disp(' ');
        disp([num2str(sum(feasible)),' (compressed) constrained MCSs found.']);
        
        mcs= expand_mcs(mcsn(:, feasible), sub)';
        disp(' ');
        disp([num2str(size(mcs,1)),' (uncompressed) constrained MCSs found.']);
    else
        mcs= expand_mcs(mcsn, sub)';
        disp([num2str(size(mcs,1)),' (uncompressed) MCSs found.']);
    end
else
    mcs=[];
end

