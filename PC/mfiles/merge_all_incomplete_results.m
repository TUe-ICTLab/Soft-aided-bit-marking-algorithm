function merge_all_incomplete_results(t)
% This function will search and merge all incomplete results

homedir=pwd;
BCH.t=t;
BCH.w=9;
BCH.iter=7;
cntr=1;
for mm=1:20
    BCH.m=mm;
    BCH.n=(2^BCH.m)-1;          % BCH Codeword length
    BCH.k=(2^BCH.m)-1-BCH.m*BCH.t;      % Number of inforamtion bit in BCH codeword
    R=(BCH.k-(BCH.n+1)/2)/((BCH.n+1)/2);
    if R>0.2 && R<0.98,mvec(cntr)=mm;cntr=cntr+1;end
end


pp=1;
for m=[1,2,3,4]
    fprintf('-------------\n');
    fprintf('m= %s \n',num2str(m));
    fprintf('-------------\n');
    for mm=mvec
        BCH.m=mm;
        BCH.n=(2^BCH.m)-1;          % BCH Codeword length
        BCH.k=(2^BCH.m)-1-BCH.m*BCH.t;      % Number of inforamtion bit in BCH codeword
        R=(BCH.k-(BCH.n+1)/2)/((BCH.n+1)/2);
        R_str=sprintf('%1.2g',R);R_str(find(R_str=='.'))='_';
        fprintf('R= %s \n',R_str);
        clear EsN0dBall
        list=dir(['./results/',R_str,'/',num2str(2^m),'PAM/incomplete/*.mat']);
        Nfiles=size(list,1);
        if Nfiles>0
            for i=1:Nfiles
                pnt1=findstr('PAM_',list(i).name);
                pnt2=findstr('_dB_',list(i).name);
                EsN0dBall(i)=str2num(list(i).name(pnt1+4:pnt2-1));
            end
            EsN0dB=unique(EsN0dBall);
            for gg=1:length(EsN0dB)
                
                PosFECBERC=merge_incomplete_results(R,R_str,m,BCH,EsN0dB(gg));
                
                if sum(PosFECBERC)<500
                    list_weak{pp}=['m=',num2str(m),', BCH.m=',num2str(BCH.m),', BCH.t=',num2str(BCH.t),', R=',R_str,', EsN0dB=',num2str(EsN0dB(gg)),', PosFECBERC= ',num2str(sum(PosFECBERC))];
                    torun{pp}=['run_simulations(',num2str(m),',',num2str(BCH.m),',',num2str(BCH.t),',',num2str(EsN0dB(gg)),');'];
                    pp=pp+1;
                end
            end
        end
    end
    
end

if pp>1
    fprintf('ooooooooooooooooooooooooooooooooooooooooooooo\n');
    fprintf('The following points need more simulations:\n');
    for ll=1:size(list_weak,2)
        fprintf('%s \n',list_weak{ll});
    end
    fprintf('ooooooooooooooooooooooooooooooooooooooooooooo\n');
    fprintf('Run the following code:\n');

    for ll=1:size(list_weak,2)
        fprintf('%s \n',torun{ll});
    end
    
    file_name=['./get_smoother_results/run_extra_sims_t_',num2str(BCH.t),'.m'];
    fid=fopen(file_name,'w');
    fprintf(fid,'cd ..\\ \n');
    for ll=1:size(list_weak,2)
        fprintf(fid,'%s \n',torun{ll});
    end
    fclose(fid);
else
    fprintf('Nothing more to do!! \n');
end

return

        