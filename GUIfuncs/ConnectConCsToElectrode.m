function E=ConnectConCsToElectrode(Econcs,conn,type)

%Econcs：构成一根电极的所有连通分量，是三维体数据
%conn：连通性
%type：连通分量的连接方式，type=1时E=logical(E+connect_neighb);type=2时E=logical(E+connect)

%E：所有连通分量连接起来后得到的电极，是三维体数据

if nargin<3
    type=1;
end

if nargin<2
    conn=26;
end

Econcs=logical(Econcs);
[Econcs_L,Econcs_n]=bwlabeln(Econcs,conn);
[nx,ny,nz]=size(Econcs);
VolumeSize=[nx,ny,nz];
succeed=0;

if Econcs_n==1
    E=Econcs;
end

if Econcs_n==2
    ConCi=Econcs_L;ConCi(Econcs_L~=1)=0;ConCi=logical(ConCi);
    ConCj=Econcs_L;ConCj(Econcs_L~=2)=0;ConCj=logical(ConCj);
    [~,~,~,~,E]=ConnectTwoConC2(ConCi,ConCj,type);
end

if Econcs_n>2
    
    [~,~,~,~,Econcs_explained,~]=PCACalculation4(Econcs);
    
    if Econcs_explained(1)>=95
        [~,~,~,~,~,Ciendpoints_sub]=ConcEndpsDetection2(Econcs);
        P=Ciendpoints_sub(1,:);
        Q=Ciendpoints_sub(2,:);
        [connect,~,connect_neighb,~]=MakeLine2(P,Q,VolumeSize);
        if type==2
            E=logical(Econcs+connect);
        else
            E=logical(Econcs+connect_neighb);
        end
        [~,E_n]=bwlabeln(E,conn);
        if E_n==1
            succeed=1;
        end
    end
    
    if succeed==0
        E=Econcs;
        conci_rec=cell(Econcs_n,1);
        conci_stats=zeros(Econcs_n,4);
        conci_stats(:,1)=(1:1:Econcs_n)';
        for i=1:Econcs_n
            conci=Econcs_L;conci(Econcs_L~=i)=0;conci=logical(conci);
            conci_rec{i}=conci;
            conci_Centroid=regionprops3(conci,'Centroid');
            conci_Centroid=conci_Centroid{1,1};
            conci_Centrsub=conci_Centroid(:,[2,1,3]);
            conci_stats(i,2:4)=conci_Centrsub;
        end
        R=conci_stats(:,2);
        C=conci_stats(:,3);
        S=conci_stats(:,4);
        Rdiff=max(R)-min(R);
        Cdiff=max(C)-min(C);
        Sdiff=max(S)-min(S);
        diff=[Rdiff,Cdiff,Sdiff];
        if max(diff)==Rdiff
            [~,RL]=sort(R,'ascend');
            conci_stats=conci_stats(RL,:);
        end
        if max(diff)==Cdiff
            [~,CL]=sort(C,'ascend');
            conci_stats=conci_stats(CL,:);
        end
        if max(diff)==Sdiff
            [~,SL]=sort(S,'ascend');
            conci_stats=conci_stats(SL,:);
        end
        for i=1:(Econcs_n-1)
            P=conci_stats(i,2:4);
            Q=conci_stats(i+1,2:4);
            [connect,~,connect_neighb,~]=MakeLine2(P,Q,VolumeSize);
            if type==2
                E=logical(E+connect);
            else
                E=logical(E+connect_neighb);
            end
        end
    end
    
end

end