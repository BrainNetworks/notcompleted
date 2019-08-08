function [VolumeOut,Stats_num,Stats_struct]=RemoveInterference(VolumeIn,Mask,conn)

%VolumeIn：输入的三维体数据
%Mask：掩模，也是三维数组，维度与VolumeIn相同，如brainmask,brainmask3等
%conn：连通性

%VolumeOut：去除干扰后得到的输出三维体数据
%Stats_num：是n*23的二维变量数组，是输入VolumeIn的所有连通分量的统计信息，n是连通分量总数，包括每一连通分量的像素数、主成分系数和贡献比以及是否属于电极等
%Stats_struct：是n*23的结构数组，是输入VolumeIn的所有连通分量的统计信息，n是连通分量总数，包括每一连通分量的像素数、主成分系数和贡献比以及是否属于电极等

if nargin<3
    conn=26;
end

VolumeIn=logical(VolumeIn);
Mask=logical(Mask);

[VolumeIn_L,VolumeIn_n]=bwlabeln(VolumeIn,conn);
VIL_intcracranial=VolumeIn_L.*Mask; %intcracranial：颅内的
VIL_cranial=VolumeIn_L.*(~Mask); %cranial：颅骨上的

VolumeInConC_ID=(1:VolumeIn_n)';

VolumeInConC_PN=zeros(VolumeIn_n,1);
VILintcConC_PN=zeros(VolumeIn_n,1); %VILintc：VIL_intcracranial
VILcConC_PN=zeros(VolumeIn_n,1); %VILc：VIL_cranial
IsElectrodes=ones(VolumeIn_n,1);
for i=1:VolumeIn_n
    VolumeInConC_PN(i)=sum(sum(sum(VolumeIn_L==i))); %ConCPN：number of pixels
    VILintcConC_PN(i)=sum(sum(sum(VIL_intcracranial==i)));
    VILcConC_PN(i)=sum(sum(sum(VIL_cranial==i)));
end
VILintcConC_PN_mean=mean(VILintcConC_PN);

[VolumeInConC_PCAcoeff,VolumeInConC_PCAlatentp]=PCACalculation3(VolumeIn_L,VolumeIn_n);
[VILintcConC_PCAcoeff,VILintcConC_PCAlatentp]=PCACalculation3(VIL_intcracranial,VolumeIn_n);
[VILcConC_PCAcoeff,VILcConC_PCAlatentp]=PCACalculation3(VIL_cranial,VolumeIn_n);

VolumeOut=VolumeIn_L;
for i=1:VolumeIn_n
    if VILintcConC_PN(i)<2
        VolumeOut(VolumeOut==i)=0;
        IsElectrodes(i)=0;
    else
        if VILintcConC_PN(i)<VILcConC_PN(i)
            if (VolumeInConC_PCAlatentp(i,1)<80)||(VILintcConC_PCAlatentp(i,1)<80)
                VolumeOut(VolumeOut==i)=0;
                IsElectrodes(i)=0;
            end
        else
            if VILintcConC_PN(i)/VILintcConC_PN_mean<0.4
                if (VolumeInConC_PCAlatentp(i,1)<80)||(VILintcConC_PCAlatentp(i,1)<80)
                    VolumeOut(VolumeOut==i)=0;
                    IsElectrodes(i)=0;
                end
            else
                if (VolumeInConC_PCAlatentp(i,1)<88)||(VILintcConC_PCAlatentp(i,1)<88)
                    VolumeOut(VolumeOut==i)=0;
                    IsElectrodes(i)=0;
                end
            end
        end
    end
end

VolumeOut=logical(VolumeOut);

Stats_num=[VolumeInConC_ID,VolumeInConC_PN,VILintcConC_PN,VILcConC_PN,...
           VolumeInConC_PCAcoeff,VILintcConC_PCAcoeff,VILcConC_PCAcoeff,...
           VolumeInConC_PCAlatentp,VILintcConC_PCAlatentp,VILcConC_PCAlatentp,...
           IsElectrodes];

VolumeInConC_ID_cell=num2cell(VolumeInConC_ID);
VolumeInConC_PN_cell=num2cell(VolumeInConC_PN);
VILintcConC_PN_cell=num2cell(VILintcConC_PN); %VILintc：VIL_intcracranial
VILcConC_PN_cell=num2cell(VILcConC_PN); %VILc：VIL_cranial
VolumeInConC_PCAcoeff_cell=num2cell(VolumeInConC_PCAcoeff);
VolumeInConC_PCAlatentp_cell=num2cell(VolumeInConC_PCAlatentp);
VILintcConC_PCAcoeff_cell=num2cell(VILintcConC_PCAcoeff);
VILintcConC_PCAlatentp_cell=num2cell(VILintcConC_PCAlatentp);
VILcConC_PCAcoeff_cell=num2cell(VILcConC_PCAcoeff);
VILcConC_PCAlatentp_cell=num2cell(VILcConC_PCAlatentp);
IsElectrodes_cell=num2cell(IsElectrodes);

field1='VolumeInConC_ID';
value1=VolumeInConC_ID_cell;
field2='VolumeInConC_PN';
value2=VolumeInConC_PN_cell;
field3='VILintcConC_PN';
value3=VILintcConC_PN_cell;
field4='VILcConC_PN';
value4=VILcConC_PN_cell;

field5='VolumeInConC_PCAcoeff1';
value5=VolumeInConC_PCAcoeff_cell(1:VolumeIn_n,1);
field6='VolumeInConC_PCAcoeff2';
value6=VolumeInConC_PCAcoeff_cell(1:VolumeIn_n,2);
field7='VolumeInConC_PCAcoeff3';
value7=VolumeInConC_PCAcoeff_cell(1:VolumeIn_n,3);

field8='VILintcConC_PCAcoeff1';
value8=VILintcConC_PCAcoeff_cell(1:VolumeIn_n,1);
field9='VILintcConC_PCAcoeff2';
value9=VILintcConC_PCAcoeff_cell(1:VolumeIn_n,2);
field10='VILintcConC_PCAcoeff3';
value10=VILintcConC_PCAcoeff_cell(1:VolumeIn_n,3);

field11='VILcConC_PCAcoeff1';
value11=VILcConC_PCAcoeff_cell(1:VolumeIn_n,1);
field12='VILcConC_PCAcoeff2';
value12=VILcConC_PCAcoeff_cell(1:VolumeIn_n,2);
field13='VILcConC_PCAcoeff3';
value13=VILcConC_PCAcoeff_cell(1:VolumeIn_n,3);

field14='VolumeInConC_PCAlatentp1';
value14=VolumeInConC_PCAlatentp_cell(1:VolumeIn_n,1);
field15='VolumeInConC_PCAlatentp2';
value15=VolumeInConC_PCAlatentp_cell(1:VolumeIn_n,2);
field16='VolumeInConC_PCAlatentp3';
value16=VolumeInConC_PCAlatentp_cell(1:VolumeIn_n,3);

field17='VILintcConC_PCAlatentp1';
value17=VILintcConC_PCAlatentp_cell(1:VolumeIn_n,1);
field18='VILintcConC_PCAlatentp2';
value18=VILintcConC_PCAlatentp_cell(1:VolumeIn_n,2);
field19='VILintcConC_PCAlatentp3';
value19=VILintcConC_PCAlatentp_cell(1:VolumeIn_n,3);

field20='VILcConC_PCAlatentp1';
value20=VILcConC_PCAlatentp_cell(1:VolumeIn_n,1);
field21='VILcConC_PCAlatentp2';
value21=VILcConC_PCAlatentp_cell(1:VolumeIn_n,2);
field22='VILcConC_PCAlatentp3';
value22=VILcConC_PCAlatentp_cell(1:VolumeIn_n,3);

field23='IsElectrodes';
value23=IsElectrodes_cell;

Stats_struct=struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,...
                    field6,value6,field7,value7,field8,value8,field9,value9,field10,value10,...
                    field11,value11,field12,value12,field13,value13,field14,value14,field15,value15,...
                    field16,value16,field17,value17,field18,value18,field19,value19,field20,value20,...
                    field21,value21,field22,value22,field23,value23);

end