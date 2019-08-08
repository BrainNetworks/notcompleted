function [MR2avg_T,ElecInfo_cell]=ContactSeg4(M_avg152T1,e17EiContacts_RASC,MR2avg_geomtform,ROI,AAL,BN_Atlas,f)
%%
waitbar(0.9,f,'Contact modeling...');

e17EiCs_avgRASC=e17EiContacts_RASC;
MR2avg_T=MR2avg_geomtform.T;
MR2avg_T=MR2avg_T'; %详见文档《affine3d》
L=length(e17EiContacts_RASC);
for i=1:L
    EiCrasc=(e17EiContacts_RASC{i})';
    C=size(EiCrasc,2);
    EiCrasc=EiCrasc([2,1,3],:); %将[R,A,S]变成[A,R,S]，因为【MRINormalization.m】和【regCTNormalization.m】
                                %中,Rfixed=imref3d(size(avg152T1),[Amin Amax],[Rmin Rmax],[Smin Smax]);
                                %其中[Amin Amax],[Rmin Rmax],[Smin Smax]的排列顺序是[A,R,S]。
    EiCrasc=[EiCrasc;ones(1,C)];
    EiCavgrasc=(MR2avg_T*EiCrasc)';
    EiCavgrasc=EiCavgrasc(:,1:3);
    EiCavgrasc=EiCavgrasc(:,[2,1,3]); %将[A;R;S]换回成[R;A;S]
    e17EiCs_avgRASC{i}=EiCavgrasc;
end
%%
[nx1,ny1,nz1]=size(AAL);
e17EiMNI=cell(L,1);
e17Ei_MNI=cell(L,1);
for i=1:L
    EiCavgrasc=(e17EiCs_avgRASC{i})';
    C=size(EiCavgrasc,2);
    EiCavgrasc=[EiCavgrasc;ones(1,C)];
    EiMNI=(M_avg152T1\EiCavgrasc)';
    EiMNI=EiMNI(:,1:3);
    EiMNI=EiMNI(:,[2,1,3]); %将坐标转换为索引
    e17EiMNI{i}=EiMNI; %不取整
    e17Ei_MNI{i}=round(EiMNI); %取整
end
%%
waitbar(0.95,f,'Determine the brain area of each contact...');

filename='AALnote.txt';
delimiterIn=' '; %空白字符等效于集合 {' ','\f','\n','\r','\t','\v'} 中的任何转义序列
AALnote=importdata(filename,delimiterIn);
AALnote_L=length(AALnote);
AAL_note=cell(AALnote_L,4);
for i=1:AALnote_L
    C=strsplit(AALnote{i});%分隔字符
    AAL_note{i,1}=C{1,1};
    AAL_note{i,2}=C{1,2};
    AAL_note{i,3}=C{1,3};
    AAL_note{i,4}=C{1,4};
end

ROI2=ROI;
for i=1:AALnote_L
    if mod(i,2)==1
        ROI2(i).Nom_L=char(strcat(ROI(i).Nom_L,' -',{32},'左脑',AAL_note{i,4})); %{32}能生成空格
    else
        ROI2(i).Nom_L=char(strcat(ROI(i).Nom_L,' -',{32},'右脑',AAL_note{i,4})); %{32}能生成空格
    end
end

ElectrodeID=[];
ContactNum=[];
Index=[];
for i=1:L
    EiMNI=e17Ei_MNI{i};
    C=size(EiMNI,1);
    EiID=i*ones(C,1);
    Cinum=(1:1:C)';
    ElectrodeID=[ElectrodeID;EiID];
    ContactNum=[ContactNum;Cinum];
    EiIndex=sub2ind(size(AAL),EiMNI(:,1),EiMNI(:,2),EiMNI(:,3));
    Index=[Index;EiIndex];
end
N=length(ElectrodeID);
Electrode_ID=cell(N,1);
Contact_Num=cell(N,1);
for i=1:N
    Electrode_ID{i}=ElectrodeID(i);
    Contact_Num{i}=ContactNum(i);
end

AALID=AAL(Index);
AAL_ID=cell(N,1);
for i=1:N
    if AALID(i)==0
        AAL_ID{i}='--';
    else
        AAL_ID{i}=AALID(i);
    end
end

[~,l]=size(ROI2);
IDroi=zeros(l,1);
for i=1:l
    IDroi(i)=ROI2(i).ID;
end
AAL_Nom_C=cell(N,1);
AAL_Nom_L=cell(N,1);
for i=1:N
    if isempty(find(IDroi==AALID(i))) 
        AAL_Nom_C{i}='--';
        AAL_Nom_L{i}='--';
    else
        AAL_Nom_C{i}=ROI2(find(IDroi==AALID(i))).Nom_C;
        AAL_Nom_L{i}=ROI2(find(IDroi==AALID(i))).Nom_L;
    end
end

field1='Electrode_ID';
value1=Electrode_ID;
field2='Contact_Num';
value2=Contact_Num;
field3='AAL_ID';
value3=AAL_ID;
field4='AAL_Nom_C';
value4=AAL_Nom_C;
field5='AAL_Nom_L';
value5=AAL_Nom_L;
Elec_AALinfo_struct=struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5);
Elec_AALinfo_cell=struct2cell(Elec_AALinfo_struct)';
%%
[nx4,ny4,nz4]=size(BN_Atlas);
% for i=1:nz4
%     figure(149);
%     imshow(squeeze(BN_Atlas(:,:,i)),[]);
% end
%将AAL插值至与BN_Atlas相同
numrows=nx4;
numcols=ny4;
numplanes=nz4;
AAL_scale=imresize3(AAL,[numrows,numcols,numplanes],'cubic'); %三维三次插值
%%
%BN_Atlas_274_with_cerebellum_without_255.nii和brodmann.nii是两种不同的分区，
%前者是中科院自动化所蒋田仔团队作的，后者是Brodmann（一个很早期的科学家）。
%BN_Atlas_274_with_cerebellum_without_255.nii与BNA_subregions.xlsx相对应。
e17Ei_BNsub=cell(L,1);
Index2=[];
for i=1:L
    EiMNI=e17EiMNI{i};
    EiBNsub=EiMNI;
    EiBNsub(:,1)=round(EiMNI(:,1)*nx4/nx1);
    EiBNsub(:,2)=round(EiMNI(:,2)*ny4/ny1);
    EiBNsub(:,3)=round(EiMNI(:,3)*nz4/nz1);
    e17Ei_BNsub{i}=EiBNsub;
    EiIndex2=sub2ind(size(BN_Atlas),EiBNsub(:,1),EiBNsub(:,2),EiBNsub(:,3));
    Index2=[Index2;EiIndex2];
end

BN_Atlas_ID=BN_Atlas(Index2);
BNA_ID=cell(N,1);
for i=1:N
    if BN_Atlas_ID(i)==0
        BNA_ID{i}='--';
    else
        BNA_ID{i}=BN_Atlas_ID(i);
    end
end

[num,~,raw]=xlsread('BNA_subregions_02.xlsx');
BNA_Lobe=cell(N,1);
BNA_Gyrus=cell(N,1);
BNA_Hemisphere=cell(N,1);
BNA_Descriptions=cell(N,1);
for i=1:N
    flag=mod(BN_Atlas_ID(i),2);
    switch flag
        case 0
            if isempty(find(num(:,2)==BN_Atlas_ID(i)))
                BNA_Lobe{i}='--';
                BNA_Gyrus{i}='--';
                BNA_Hemisphere{i}='--';
                BNA_Descriptions{i}='--';
            else
                %char用于转换为字符变量，以便写入excel表格
                BNA_Lobe{i}=char(raw(find(num(:,2)==BN_Atlas_ID(i))+1,1)); %'Lobe'
                BNA_Gyrus{i}=char(strcat('Right Hemi -',{32},raw(find(num(:,2)==BN_Atlas_ID(i))+1,2))); %'Gyrus',{32}能生成空格
                BNA_Hemisphere{i}=char(raw(find(num(:,2)==BN_Atlas_ID(i))+1,4)); %'Hemisphere'
                BNA_Descriptions{i}=char(raw(find(num(:,2)==BN_Atlas_ID(i))+1,7)); %'Descriptions'
            end
        case 1
            if isempty(find(num(:,1)==BN_Atlas_ID(i)))
                BNA_Lobe{i}='--';
                BNA_Gyrus{i}='--';
                BNA_Hemisphere{i}='--';
                BNA_Descriptions{i}='--';
            else
                %char用于转换为字符变量，以便写入excel表格
                BNA_Lobe{i}=char(raw(find(num(:,1)==BN_Atlas_ID(i))+1,1)); %'Lobe'
                BNA_Gyrus{i}=char(strcat('Left Hemi -',{32},raw(find(num(:,1)==BN_Atlas_ID(i))+1,2))); % 'Gyrus',{32}能生成空格
                BNA_Hemisphere{i}=char(raw(find(num(:,1)==BN_Atlas_ID(i))+1,3)); %'Hemisphere'
                BNA_Descriptions{i}=char(raw(find(num(:,1)==BN_Atlas_ID(i))+1,7)); %'Descriptions'
            end
    end
end

field1='Electrode_ID';
value1=Electrode_ID;
field2='Contact_Num';
value2=Contact_Num;
field6='BNA_ID';
value6=BNA_ID;
field7='BNA_Lobe';
value7=BNA_Lobe;
field8='BNA_Gyrus';
value8=BNA_Gyrus;
field9='BNA_Hemisphere';
value9=BNA_Hemisphere;
field10='BNA_Descriptions';
value10=BNA_Descriptions;
Elec_BNAinfo_struct=struct(field1,value1,field2,value2,field6,value6,...
    field7,value7,field8,value8,field9,value9,field10,value10);
Elec_BNAinfo_cell=struct2cell(Elec_BNAinfo_struct)';
%%
e17EiContactsR=[];
e17EiContactsA=[];
e17EiContactsS=[];
for i=1:L
    EiCrasc=e17EiContacts_RASC{i};
    e17EiContactsR=[e17EiContactsR;EiCrasc(:,1)];
    e17EiContactsA=[e17EiContactsA;EiCrasc(:,2)];
    e17EiContactsS=[e17EiContactsS;EiCrasc(:,3)];    
end
e17EiContacts_R=num2cell(e17EiContactsR);
e17EiContacts_A=num2cell(e17EiContactsA);
e17EiContacts_S=num2cell(e17EiContactsS);
field11='R';
value11=e17EiContacts_R;
field12='A';
value12=e17EiContacts_A;
field13='S';
value13=e17EiContacts_S;

MNI_Index=num2cell(Index);
field14='MNI_Index';
value14=MNI_Index;

ElecInfo_struct=struct(field1,value1,field2,value2,...
    field11,value11,field12,value12,field13,value13,...
    field14,value14,...
    field3,value3,field4,value4,field5,value5,...
    field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
ElecInfo_cell=struct2cell(ElecInfo_struct)';
% ElecInfo_table=struct2table(ElecInfo_struct);

end