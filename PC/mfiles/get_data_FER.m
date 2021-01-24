function [PreFECBERav,PosFECBERav,MIav,GMIav,HDMIav,m,EsN0dB]= get_data_FER(R_str,m,BCH)


pref_save_dir   = strcat('results/',R_str,'/',num2str(2^m),'PAM/');
if ispc,pref_save_dir(pref_save_dir=='/')='\';end

list=dir([pref_save_dir,num2str(2^m),'PAM_*_dB_m_',num2str(BCH.m),'_t_',num2str(BCH.t),'_w_',num2str(BCH.w),'_iter_',num2str(BCH.iter),'*']);
LL=length(list());
fprintf('Found %i files \n',LL);
PreFECBER=[];
PreFECFER=[];
PosFECBER=[];
EsN0dB=[];
PreFECBERav=[];
PosFECBERav=[];
MI=[];
GMI=[];
MIav=[];
GMIav=[];
HDMIav=[];
FERav=[];
if LL==0
    return
end
for ll=1:LL
    file=load(strcat(pref_save_dir,list(ll).name));
    %PreFECBERtmp{ll}=file.PreFECBER;
    %PosFECBERtmp{ll}=file.PosFECBER;
    EsN0dB=[EsN0dB,file.EsN0dB];
    PreFECBERav=[PreFECBERav,file.PreFECBERav];
    PosFECBERav=[PosFECBERav,file.PosFECBERav];
    MIav=[MIav;file.MIav];
    GMIav=[GMIav;file.GMIav];
    HDMIav=[HDMIav;file.HDMIav];
    %FERav=[FERav;file.FERav];
end
[A,B]=sort(EsN0dB);
EsN0dB=EsN0dB(B);

for ll=1:LL
    %PreFECBER{ll}=PreFECBERtmp{B(ll)};
    %PosFECBER{ll}=PosFECBERtmp{B(ll)};
    clrrnd{ll}=rand(3,1);
end
PreFECBERav=PreFECBERav(B);
PosFECBERav=PosFECBERav(B);
MIav=MIav(B);
GMIav=GMIav(B);
HDMIav=HDMIav(B);
%FERav=FERav(B);

m=file.m;
return