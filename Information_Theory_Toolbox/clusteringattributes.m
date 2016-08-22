clear
clc
[YYYYMMDD,HHMMLST,Zenithdeg,Azimuthdeg,ETRWm2,ETRNWm2,GloModWm2,GloModUnc,GloModSource,DirModWm2,DirModUnc,DirModSource,DifModWm2,DifModUnc,DifModSource,MeasGloWm2,MeasGloFlg,MeasDirWm2,MeasDirFlg,MeasDifWm2,MeasDifFlg,TotCC10ths,PrecipWatcm,PrecipWatFlg,AODunitless,AODFlg]=importfile('E:\UCLA\Dataset\2005.csv',2,8761);
multi=[Zenithdeg,Azimuthdeg,ETRWm2,ETRNWm2,GloModWm2,GloModUnc,GloModSource,DirModWm2,DirModUnc,DirModSource,DifModWm2,DifModUnc,DifModSource,MeasGloWm2,MeasGloFlg,MeasDirWm2,MeasDirFlg,MeasDifWm2,MeasDifFlg,TotCC10ths,PrecipWatcm,PrecipWatFlg,AODunitless,AODFlg];
%multi=[Zenithdeg,ETRWm2,ETRNWm2,GloModWm2,GloModUnc,GloModSource,DirModWm2,DirModUnc,DirModSource,DifModWm2,DifModUnc,DifModSource,MeasGloFlg,MeasDirFlg,MeasDifFlg,PrecipWatcm,PrecipWatFlg,AODunitless,AODFlg];
[D,select]=size(multi);
multi_norm=multi;
for i=1:1:select
%     minmulti=min(multi(:,i));
%     maxmulti=max(multi(:,i));
%     meanmulti=mean(multi(:,i));
%     varmulti=var(multi(:,i));
%     if minmulti==maxmulti
%          multi_norm(:,i)=128;
%     else
%     multi_norm(:,i)=multi_norm(:,i)+abs(minmulti)+1;
%     norm=(multi(:,i)-meanmulti)/varmulti;
%     multi_norm(:,i)=round(norm)*256+abs(minmulti);
%     end
   tab=tabulate(multi(:,i));
   [tablen,tabwid]=size(tab);
   for j=1:1:D
       for k=1:1:tablen
           if multi(j,i)==tab(k,1) 
               multi_norm(j,i)=k;
           end
       end
   end
end
combine=combntns(1:select,2);
comblen=length(combine);
MInfo=zeros(comblen,1);
Entro=zeros(comblen,1);
Holo=zeros(comblen,1);
% CT=zeros(comblen,1);
% CSd=zeros(comblen,1);
for i=1:1:comblen
    a=multi_norm(:,combine(i,1));
    b=multi_norm(:,combine(i,2));
    MInfo(i)=D*mutInfo(a,b);
    Entro(i)=entropy(a)+entropy(b);
    Holo(i)=MInfo(i)+Entro(i);
%     CTa=log2(max(a)-min(a))+log2(log2(D))-logfr(a);
%     CTb=log2(max(b)-min(b))+log2(log2(D))-logfr(b);
%     CT(i)=CTa-CTb;
%     CSd(i)=MInfo(i)+CT(i);
end
MiniHolo=min(Holo);
% MiniHoloNumber=find(Holo==MiniHolo);
MiniHoloNumber=find(Holo<0.06);
MiniHoloLen=length(MiniHoloNumber);
UnionCombine=zeros(MiniHoloLen,2);
for i=1:1:MiniHoloLen
    UnionCombine(i,1)=combine(MiniHoloNumber(i),1);
    UnionCombine(i,2)=combine(MiniHoloNumber(i),2);
end
UnionAll=[]; %zeros(1,select);
for j=1:1:MiniHoloLen
    UnionAll=union(UnionAll,UnionCombine(j,:));
    UnionAllcombine=combntns(UnionAll,2);
    
end

    UnionAllcomblen=length(UnionAllcombine);
    UnionAllMInfo=zeros(UnionAllcomblen,1);
    UnionAllEntro=zeros(UnionAllcomblen,1);
    UnionAllHolo=zeros(UnionAllcomblen,1);

    for i=1:1:UnionAllcomblen
    a=multi_norm(:,UnionAllcombine(i,1));
    b=multi_norm(:,UnionAllcombine(i,2));
    UnionAllMInfo(i)=D*mutInfo(a,b);
    UnionAllEntro(i)=entropy(a)+entropy(b);
    UnionAllHolo(i)=UnionAllMInfo(i)+UnionAllEntro(i);
    end
   
   ClusterMiniHoloNumber=find(UnionAllHolo<0.06);
 ClusterMiniHoloLen=length( ClusterMiniHoloNumber);
 ClusterUnionCombine=zeros( ClusterMiniHoloLen,2);
for i=1:1: ClusterMiniHoloLen
     ClusterUnionCombine(i,1)=UnionAllcombine( ClusterMiniHoloNumber(i),1);
     ClusterUnionCombine(i,2)=UnionAllcombine( ClusterMiniHoloNumber(i),2);
end
 ClusterUnionAll=[]; %zeros(1,select);
for j=1:1: ClusterMiniHoloLen
     ClusterUnionAll=union( ClusterUnionAll, ClusterUnionCombine(j,:));
end
cluster=zeros(D,length(ClusterUnionAll));
for i=1:1:length(ClusterUnionAll)
    cluster(:,i)=multi(:,ClusterUnionAll(i,:));
end
left_number=setdiff((1:1:24)',ClusterUnionAll);
left_len=length(left_number);
left_cluster=zeros(D,left_len);
for i=1:1:left_len
    left_cluster(:,i)=multi(:,left_number(i,:));
end
left_Sparse=zeros(D,left_len);
for i=1:1:length(left_number)
    left_Sparse(:,i)=OF_cal(multi(:,left_number(i,:)));
end
left_apriori=[];
for i=1:1:D
    if mean(left_Sparse(i,:))>0
        left_apriori=[left_apriori;left_Sparse(i,:)];
    end
end
load('code.mat')
[Rules,FreqItemsets] = findRules(left_apriori, 0.001, 0.5, 1000, 1, code, 'rules_multi.txt');
FPtree={};
for j=1:1:length(left_apriori)
    FPtree{j,1}=find(left_apriori(j,:)==1);
end
