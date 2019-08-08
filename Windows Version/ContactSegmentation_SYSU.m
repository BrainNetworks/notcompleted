function varargout = ContactSegmentation_SYSU(varargin)
% CONTACTSEGMENTATION_SYSU MATLAB code for ContactSegmentation_SYSU.fig
%      CONTACTSEGMENTATION_SYSU, by itself, creates a new CONTACTSEGMENTATION_SYSU or raises the existing
%      singleton*. 
% 
%      H = CONTACTSEGMENTATION_SYSU returns the handle to a new CONTACTSEGMENTATION_SYSU or the handle to
%      the existing singleton*. 
%
%      CONTACTSEGMENTATION_SYSU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTACTSEGMENTATION_SYSU.M with the given input arguments.
%  
%      CONTACTSEGMENTATION_SYSU('Property','Value',...) creates a new CONTACTSEGMENTATION_SYSU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ContactSegmentation_SYSU_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ContactSegmentation_SYSU_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ContactSegmentation_SYSU

% Last Modified by GUIDE v2.5 23-Jul-2019 22:25:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ContactSegmentation_SYSU_OpeningFcn, ...
                   'gui_OutputFcn',  @ContactSegmentation_SYSU_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ContactSegmentation_SYSU is made visible.
function ContactSegmentation_SYSU_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ContactSegmentation_SYSU (see VARARGIN)

% Choose default command line output for ContactSegmentation_SYSU
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ContactSegmentation_SYSU wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% set(gcf,'position',get(0,'screensize'));
% set(gcf,'outerposition',get(0,'screensize'));

global GUIfuncs_path
GUIfuncs_path=[pwd,'\GUIfuncs'];
addpath(GUIfuncs_path);

global spm12path
spm12path=[pwd,'\spm12'];
addpath(spm12path);

global StartImage ratio
StartImage=ones(700,700)*0.9;
ratio=2;

axes(handles.AxialPlane);
imshow(StartImage,[0 1]);
axes(handles.SagittalPlane);
imshow(StartImage,[0 1]);
axes(handles.CoronalPlane);
imshow(StartImage,[0 1]);
axes(handles.TDPlane);
imshow(StartImage,[0 1]);

set(handles.SaveDataType,'String',{'Please choose data to save：'});

global NonReconJudge ReconstructionJudge RASTransformJudge RegistrationJudge SegmentationJudge
NonReconJudge=0;ReconstructionJudge=0;RASTransformJudge=0;RegistrationJudge=0;SegmentationJudge=0;

global dicom_flag nii_flag mat_flag
dicom_flag=0;nii_flag=0;mat_flag=0;


% --- Outputs from this function are returned to the command line.
function varargout = ContactSegmentation_SYSU_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in ImportDataType.
function ImportDataType_Callback(hObject, eventdata, handles)
% hObject    handle to ImportDataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ImportDataType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ImportDataType

global AllType DataType DataType_Label
DataType_Label=get(hObject,'Value');
DataType=AllType{DataType_Label};


% --- Executes during object creation, after setting all properties.
function ImportDataType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImportDataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global AllType
AllType=get(hObject,'String');


% --- Executes on button press in DataImport.
function DataImport_Callback(hObject, eventdata, handles)
% hObject    handle to DataImport (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% [file1, path1, index1] = uigetfile( ...
% {'*.m;*.mlx;*.fig;*.mat;*.slx;*.mdl',...
%     'MATLAB Files (*.m,*.mlx,*.fig,*.mat,*.slx,*.mdl)';
%    '*.m;*.mlx','Code files (*.m,*.mlx)'; ...
%    '*.fig','Figures (*.fig)'; ...
%    '*.mat','MAT-files (*.mat)'; ...
%    '*.mdl;*.slx','Models (*.slx, *.mdl)'; ...
%    '*.*',  'All Files (*.*)'}, ...
%    'Select a File', ...
%    'MultiSelect', 'on');

global NonReconJudge ReconstructionJudge RASTransformJudge RegistrationJudge SegmentationJudge
NonReconJudge=1;ReconstructionJudge=0;RASTransformJudge=0;RegistrationJudge=0;SegmentationJudge=0;

global DataType Images Images2 Images3 Images_IPP1 Images_IOP Images_PixelSpacing Images_SliceThickness 
global FolderName selpath Fname1 Fname2_all Fname2 path1 path2 DataType_Label dicom_flag nii_flag mat_flag

if DataType_Label<2
    f=warndlg('Please choose data type first.','Warning');
else   
    
    if strcmp(DataType,'dicom')
        
        if isempty(selpath)
            selpath='D:\';
        end
        selpath_old=selpath;
        selpath = uigetdir(selpath,'Select a Folder');
        if selpath~=0
            dicom_flag=1;nii_flag=0;mat_flag=0;
            FolderName=[];Fname1=[];Fname2=[];Fname2_all=[];

            set(handles.ElecNum,'String',[]);
            set(handles.ChooseElectrode,'Value',1);
            set(handles.ChooseElectrode,'String','Please choose electrode：');
            set(handles.ChooseContact,'Value',1);
            set(handles.ChooseContact,'String','Please choose contact：');
            set(handles.RASCoordinate,'String',[]);
            set(handles.MNIIndex,'String',[]);
            set(handles.AALBrainArea,'String',[]);
            set(handles.BNBrainArea,'String',[]);

            set(handles.ChooseInterpolation,'Enable','on');
            set(handles.DataReconstruction,'Enable','on');
            set(handles.ChooseCoordinate,'Enable','on');
            set(handles.CoordinateTransformation,'Enable','on');
            
            for i=1:length(selpath)
                if strcmp(selpath(length(selpath)-i+1),'\')
                    break;
                else
                    FolderName=strcat(selpath(length(selpath)-i+1),FolderName);
                end
            end
            selpath=[selpath,'\'];

            [Images,Images_IPP1,Images_IOP,Images_PixelSpacing,Images_SliceThickness]=ReadDicomData2(selpath);

            global T1_path T1 T1_IPP1 T1_IOP T1_PixelSpacing T1_SliceThickness
            if strcmp(FolderName,'MR_T1')
                T1_path=selpath;
                T1=Images;T1_IPP1=Images_IPP1;T1_IOP=Images_IOP;T1_PixelSpacing=Images_PixelSpacing;T1_SliceThickness=Images_SliceThickness;
                SaveDataStr=get(handles.SaveDataType,'String');
                L=size(SaveDataStr,1);
                SaveDataStr2=cell(L+1,1);
                for i=1:L
                    SaveDataStr2{i,1}=SaveDataStr{i,1};
                end
                SaveDataStr2{L+1,1}='MR_T1';
                set(handles.SaveDataType,'String',SaveDataStr2);
            end

            global CT_path CT CT_IPP1 CT_IOP CT_PixelSpacing CT_SliceThickness
            if strcmp(FolderName,'CT')
                CT_path=selpath;
                CT=Images;CT_IPP1=Images_IPP1;CT_IOP=Images_IOP;CT_PixelSpacing=Images_PixelSpacing;CT_SliceThickness=Images_SliceThickness;
                SaveDataStr=get(handles.SaveDataType,'String');
                L=size(SaveDataStr,1);
                SaveDataStr2=cell(L+1,1);
                for i=1:L
                    SaveDataStr2{i,1}=SaveDataStr{i,1};
                end
                SaveDataStr2{L+1,1}='CT';
                set(handles.SaveDataType,'String',SaveDataStr2);
            end

            global fsmask_path fsmask fsmask_IPP1 fsmask_IOP fsmask_PixelSpacing fsmask_SliceThickness
            if strcmp(FolderName,'FsBrainmask')
                fsmask_path=selpath;
                fsmask=Images;fsmask_IPP1=Images_IPP1;fsmask_IOP=Images_IOP;fsmask_PixelSpacing=Images_PixelSpacing;fsmask_SliceThickness=Images_SliceThickness;
                SaveDataStr=get(handles.SaveDataType,'String');
                L=size(SaveDataStr,1);
                SaveDataStr2=cell(L+1,1);
                for i=1:L
                    SaveDataStr2{i,1}=SaveDataStr{i,1};
                end
                SaveDataStr2{L+1,1}='FsBrainmask';
                set(handles.SaveDataType,'String',SaveDataStr2);
            end
            selpath=selpath(1:(end-length(FolderName)-1));
        else
            dicom_flag=0;nii_flag=0;mat_flag=0;
            selpath=selpath_old;
        end
        
    end

    if strcmp(DataType,'nii')

        if isempty(path1)
            path1=[pwd,'\'];
        end
        Fname1_old=Fname1;
        path1_old=path1;
        [Fname1,path1] = uigetfile([path1,'.nii'],'Select an Nii File');
        if path1~=0
            dicom_flag=0;nii_flag=1;mat_flag=0;            
            Fname2_all=[];Fname2=[];FolderName=[];
            
            set(handles.ElecNum,'String',[]);
            set(handles.ChooseElectrode,'Value',1);
            set(handles.ChooseElectrode,'String','Please choose electrode：');
            set(handles.ChooseContact,'Value',1);
            set(handles.ChooseContact,'String','Please choose contact：');
            set(handles.RASCoordinate,'String',[]);
            set(handles.MNIIndex,'String',[]);
            set(handles.AALBrainArea,'String',[]);
            set(handles.BNBrainArea,'String',[]);

            set(handles.ChooseInterpolation,'Enable','on');
            set(handles.DataReconstruction,'Enable','on');
            set(handles.ChooseCoordinate,'Enable','on');
            set(handles.CoordinateTransformation,'Enable','on');
            
            f = waitbar(0,'Data importing. Please wait...');

            WholePath=fullfile(path1,Fname1);
            Images_struct=spm_vol(WholePath);waitbar(0.5,f);
            [Images,~]=spm_read_vols(Images_struct);waitbar(1,f);
            
            global avg152T1_path avg152T1 avg152T1_IPP1 avg152T1_IOP avg152T1_PixelSpacing avg152T1_SliceThickness
            if strcmp(Fname1,'avg152T1.nii')
                avg152T1_path=[WholePath,'\'];
                avg152T1=Images;
                avg152T1_PixelSpacing=[2;2];
                avg152T1_SliceThickness=2;
                avg152T1_IOP=[0;-1;0;1;0;0]; %见《MNI空间说明》文档第9页
                avg152T1_IPP1=[-90,126,-72]; %见《MNI空间说明》文档第9页
            end
            Images_IPP1=avg152T1_IPP1;Images_IOP=avg152T1_IOP;Images_PixelSpacing=avg152T1_PixelSpacing;Images_SliceThickness=avg152T1_SliceThickness;
            SaveDataStr=get(handles.SaveDataType,'String');
            L=size(SaveDataStr,1);
            SaveDataStr2=cell(L+1,1);
            for i=1:L
                SaveDataStr2{i,1}=SaveDataStr{i,1};
            end
            SaveDataStr2{L+1,1}='avg152T1';
            set(handles.SaveDataType,'String',SaveDataStr2);
            
            close(f);
        else
            dicom_flag=0;nii_flag=0;mat_flag=0;
            Fname1=Fname1_old;
            path1=path1_old;
        end
        
    end

    global nx1 ny1 nz1
    global IrAV IrSV IrCV
    global Axial_WinMin Axial_WinMax AxialColormap AxialColormap_Label AxialAlpha
    global Sagittal_WinMin Sagittal_WinMax SagittalColormap SagittalColormap_Label SagittalAlpha
    global Coronal_WinMin Coronal_WinMax CoronalColormap CoronalColormap_Label CoronalAlpha
    global TDColormap TDColormap_Label TDAlpha
    global Colormap AlphaValue Plane
    global StartImage

    if strcmp(DataType,'mat')

        if isempty(path2)
            path2=[pwd,'\'];
        end
        Fname2_all_old=Fname2_all;
        path2_old=path2;
        [Fname2_all,path2] = uigetfile([path2,'.mat'],'Select an MATLAB-mat File','MultiSelect', 'on');
        if path2~=0
            dicom_flag=0;nii_flag=0;mat_flag=1;
            Fname1=[];FolderName=[];
            
            f = waitbar(0,'Data importing. Please wait...');
            
            if ischar(Fname2_all)
                Fname2_all2=Fname2_all;
                Fname2_all=cell(1,1);
                Fname2_all{1,1}=Fname2_all2;
            end
            for i=1:length(Fname2_all)
                Fname2=Fname2_all{1,i};
                WholePath=fullfile(path2,Fname2);
                if (~strcmp(Fname2,'T12avg_geomtform.mat'))&&(~strcmp(Fname2,'regCT2avg_geomtform.mat'))...
                        &&(~strcmp(Fname2,'e17EiContacts_sub.mat'))&&(~strcmp(Fname2,'e17EiContacts_RASC.mat'))&&(~strcmp(Fname2,'ElecInfo_cell.mat'))

                    set(handles.ElecNum,'String',[]);
                    set(handles.ChooseElectrode,'Value',1);
                    set(handles.ChooseElectrode,'String','Please choose electrode：');
                    set(handles.ChooseContact,'Value',1);
                    set(handles.ChooseContact,'String','Please choose contact：');
                    set(handles.RASCoordinate,'String',[]);
                    set(handles.MNIIndex,'String',[]);
                    set(handles.AALBrainArea,'String',[]);
                    set(handles.BNBrainArea,'String',[]);

                    set(handles.ChooseInterpolation,'Enable','off');
                    set(handles.DataReconstruction,'Enable','off');
                    set(handles.ChooseCoordinate,'Enable','off');
                    set(handles.CoordinateTransformation,'Enable','off');

                    set(handles.ChoosePlane,'Enable','on');
                    set(handles.WindowMin,'Enable','on');
                    set(handles.WindowMax,'Enable','on');
                    set(handles.ChooseColormap,'Enable','on');
                    set(handles.Trans,'Enable','on');

                    Images=load(WholePath);
                    Images=struct2cell(Images);
                    Images=Images{1,1};
                    Images=single(Images);                    
                    [nx1,ny1,nz1]=size(Images);

                    set(handles.AxialSlider,'Max',nx1);
                    set(handles.SagittalSlider,'Max',nz1);
                    set(handles.CoronalSlider,'Max',ny1);

                    set(handles.AxialSlider,'Value',1);
                    set(handles.SagittalSlider,'Value',1);
                    set(handles.CoronalSlider,'Value',1);

                    set(handles.AxialSlider,'Min',1);
                    set(handles.SagittalSlider,'Min',1);
                    set(handles.CoronalSlider,'Min',1);

                    set(handles.AxialSlider,'SliderStep',[0.02 0.1]);
                    set(handles.SagittalSlider,'SliderStep',[0.02 0.1]);
                    set(handles.CoronalSlider,'SliderStep',[0.02 0.1]);

                    graymin=min(min(min(Images)));
                    graymax=max(max(max(Images)));

                    if strcmp(Fname2,'T12avg.mat')||strcmp(Fname2,'regCT2avg.mat')
                        Images3=Images;
                        Images=InvRotateAvg(Images);
                    end

                    SagCtr=Images(:,:,round(nz1/2));
                    CorCtr=fliplr(squeeze(Images(:,round(ny1/2),:)));
                    AxiCtr=fliplr(squeeze(Images(round(nx1/2),:,:)));

                    IrAV=AxiCtr;
                    IrSV=SagCtr;
                    IrCV=CorCtr;

                    axes(handles.SagittalPlane);
                    imshow(SagCtr,[graymin graymax]);
                    axes(handles.CoronalPlane);
                    imshow(CorCtr,[graymin graymax]);
                    axes(handles.AxialPlane);
                    imshow(AxiCtr,[graymin graymax]);

                    set(handles.AxialSlider,'Value',round(nx1/2));
                    set(handles.SagittalSlider,'Value',round(nz1/2));
                    set(handles.CoronalSlider,'Value',round(ny1/2));

                    set(handles.WindowMin,'Max',graymax);
                    set(handles.WindowMin,'Value',graymax);
                    set(handles.WindowMin,'Min',graymin);
                    set(handles.WindowMin,'Value',graymin);

                    set(handles.WindowMax,'Max',graymax);
                    set(handles.WindowMax,'Value',graymax);
                    set(handles.WindowMax,'Min',graymin);

                    Axial_WinMin=graymin;Axial_WinMax=graymax;AxialColormap=gray;AxialColormap_Label=1;AxialAlpha=1;
                    Sagittal_WinMin=graymin;Sagittal_WinMax=graymax;SagittalColormap=gray;SagittalColormap_Label=1;SagittalAlpha=1;
                    Coronal_WinMin=graymin;Coronal_WinMax=graymax;CoronalColormap=gray;CoronalColormap_Label=1;CoronalAlpha=1;
                    TDColormap=gray;TDColormap_Label=1;TDAlpha=1;

                    Plane='Axial Plane';Colormap=gray;AlphaValue=1;
                    set(handles.ChoosePlane,'Value',1);
                    set(handles.ChooseColormap,'Value',1);
                    set(handles.Trans,'String',[]);

                    axes(handles.TDPlane);
                    imshow(StartImage,[0 1]);
                    set(handles.TDEdit,'String',[]);
                    
                else
                    Images2=load(WholePath);
                    Images2=struct2cell(Images2);
                    Images2=Images2{1,1};            
                end

                if strcmp(Fname2,'regCT.mat')
                    global regCT
                    regCT=Images;
                    SaveDataStr=get(handles.SaveDataType,'String');
                    L=size(SaveDataStr,1);
                    SaveDataStr2=cell(L+1,1);
                    for i=1:L
                        SaveDataStr2{i,1}=SaveDataStr{i,1};
                    end
                    SaveDataStr2{L+1,1}='regCT';
                    set(handles.SaveDataType,'String',SaveDataStr2);
                end

                if strcmp(Fname2,'mask2regCT.mat')
                    global mask2regCT
                    mask2regCT=Images;
                    SaveDataStr=get(handles.SaveDataType,'String');
                    L=size(SaveDataStr,1);
                    SaveDataStr2=cell(L+1,1);
                    for i=1:L
                        SaveDataStr2{i,1}=SaveDataStr{i,1};
                    end
                    SaveDataStr2{L+1,1}='mask2regCT';
                    set(handles.SaveDataType,'String',SaveDataStr2);
                end

                if strcmp(Fname2,'T12avg.mat')
                    global T12avg
                    T12avg=Images3;
                    SaveDataStr=get(handles.SaveDataType,'String');
                    L=size(SaveDataStr,1);
                    SaveDataStr2=cell(L+1,1);
                    for i=1:L
                        SaveDataStr2{i,1}=SaveDataStr{i,1};
                    end
                    SaveDataStr2{L+1,1}='T12avg';
                    set(handles.SaveDataType,'String',SaveDataStr2);
                end

                if strcmp(Fname2,'regCT2avg.mat')
                    global regCT2avg
                    regCT2avg=Images3;
                    SaveDataStr=get(handles.SaveDataType,'String');
                    L=size(SaveDataStr,1);
                    SaveDataStr2=cell(L+1,1);
                    for i=1:L
                        SaveDataStr2{i,1}=SaveDataStr{i,1};
                    end
                    SaveDataStr2{L+1,1}='regCT2avg';
                    set(handles.SaveDataType,'String',SaveDataStr2);
                end

                if strcmp(Fname2,'T12avg_geomtform.mat')
                    global T12avg_geomtform
                    T12avg_geomtform=Images2;
                    SaveDataStr=get(handles.SaveDataType,'String');
                    L=size(SaveDataStr,1);
                    SaveDataStr2=cell(L+1,1);
                    for i=1:L
                        SaveDataStr2{i,1}=SaveDataStr{i,1};
                    end
                    SaveDataStr2{L+1,1}='T12avg_geomtform';
                    set(handles.SaveDataType,'String',SaveDataStr2);
                end

                if strcmp(Fname2,'regCT2avg_geomtform.mat')
                    global regCT2avg_geomtform
                    regCT2avg_geomtform=Images2;
                    SaveDataStr=get(handles.SaveDataType,'String');
                    L=size(SaveDataStr,1);
                    SaveDataStr2=cell(L+1,1);
                    for i=1:L
                        SaveDataStr2{i,1}=SaveDataStr{i,1};
                    end
                    SaveDataStr2{L+1,1}='regCT2avg_geomtform';
                    set(handles.SaveDataType,'String',SaveDataStr2);
                end

                if strcmp(Fname2,'e17EiContacts_sub.mat')
                    global e17EiContacts_sub
                    e17EiContacts_sub=Images2;
                    SaveDataStr=get(handles.SaveDataType,'String');
                    L=size(SaveDataStr,1);
                    SaveDataStr2=cell(L+1,1);
                    for i=1:L
                        SaveDataStr2{i,1}=SaveDataStr{i,1};
                    end
                    SaveDataStr2{L+1,1}='e17EiContacts_sub';
                    set(handles.SaveDataType,'String',SaveDataStr2);
                end

                if strcmp(Fname2,'e17EiContacts_RASC.mat')
                    global e17EiContacts_RASC
                    e17EiContacts_RASC=Images2;
                    SaveDataStr=get(handles.SaveDataType,'String');
                    L=size(SaveDataStr,1);
                    SaveDataStr2=cell(L+1,1);
                    for i=1:L
                        SaveDataStr2{i,1}=SaveDataStr{i,1};
                    end
                    SaveDataStr2{L+1,1}='e17EiContacts_RASC';
                    set(handles.SaveDataType,'String',SaveDataStr2);
                end

                if strcmp(Fname2,'ElecInfo_cell.mat')
                    global ElecInfo_cell
                    ElecInfo_cell=Images2;
                    SaveDataStr=get(handles.SaveDataType,'String');
                    L=size(SaveDataStr,1);
                    SaveDataStr2=cell(L+1,1);
                    for i=1:L
                        SaveDataStr2{i,1}=SaveDataStr{i,1};
                    end
                    SaveDataStr2{L+1,1}='ElecInfo_cell';
                    set(handles.SaveDataType,'String',SaveDataStr2);
                end
            end
            waitbar(1,f);
            close(f);
        else
            dicom_flag=0;nii_flag=0;mat_flag=0;
            Fname2_all=Fname2_all_old;
            path2=path2_old;
        end
        
    end

    if (dicom_flag==1)||(nii_flag==1)

        set(handles.ChoosePlane,'Enable','on');
        set(handles.WindowMin,'Enable','on');
        set(handles.WindowMax,'Enable','on');
        set(handles.ChooseColormap,'Enable','on');
        set(handles.Trans,'Enable','on');

        [nx1,ny1,nz1]=size(Images);

        set(handles.AxialSlider,'Max',nx1);
        set(handles.SagittalSlider,'Max',nz1);
        set(handles.CoronalSlider,'Max',ny1);

        set(handles.AxialSlider,'Value',1);
        set(handles.SagittalSlider,'Value',1);
        set(handles.CoronalSlider,'Value',1);

        set(handles.AxialSlider,'Min',1);
        set(handles.SagittalSlider,'Min',1);
        set(handles.CoronalSlider,'Min',1);

        set(handles.AxialSlider,'SliderStep',[0.02 0.1]);
        set(handles.SagittalSlider,'SliderStep',[0.02 0.1]);
        set(handles.CoronalSlider,'SliderStep',[0.02 0.1]);

        graymin=min(min(min(Images)));
        graymax=max(max(max(Images)));

        if strcmp(Fname1,'avg152T1.nii')
            Images=InvRotateAvg(Images);
        end

        SagCtr=Images(:,:,round(nz1/2));
        CorCtr=fliplr(squeeze(Images(:,round(ny1/2),:)));
        AxiCtr=fliplr(squeeze(Images(round(nx1/2),:,:)));

        IrAV=AxiCtr;
        IrSV=SagCtr;
        IrCV=CorCtr;

        axes(handles.SagittalPlane);
        imshow(SagCtr,[graymin graymax]);
        axes(handles.CoronalPlane);
        imshow(CorCtr,[graymin graymax]);
        axes(handles.AxialPlane);
        imshow(AxiCtr,[graymin graymax]);

        set(handles.AxialSlider,'Value',round(nx1/2));
        set(handles.SagittalSlider,'Value',round(nz1/2));
        set(handles.CoronalSlider,'Value',round(ny1/2));

        set(handles.WindowMin,'Max',graymax);
        set(handles.WindowMin,'Value',graymax);
        set(handles.WindowMin,'Min',graymin);
        set(handles.WindowMin,'Value',graymin);

        set(handles.WindowMax,'Max',graymax);
        set(handles.WindowMax,'Value',graymax);
        set(handles.WindowMax,'Min',graymin);

        Axial_WinMin=graymin;Axial_WinMax=graymax;AxialColormap=gray;AxialColormap_Label=1;AxialAlpha=1;
        Sagittal_WinMin=graymin;Sagittal_WinMax=graymax;SagittalColormap=gray;SagittalColormap_Label=1;SagittalAlpha=1;
        Coronal_WinMin=graymin;Coronal_WinMax=graymax;CoronalColormap=gray;CoronalColormap_Label=1;CoronalAlpha=1;
        TDColormap=gray;TDColormap_Label=1;TDAlpha=1;

        Plane='Axial Plane';Colormap=gray;AlphaValue=1;
        set(handles.ChoosePlane,'Value',1);
        set(handles.ChooseColormap,'Value',1);
        set(handles.Trans,'String',[]);

        axes(handles.TDPlane);
        imshow(StartImage,[0 1]);
        set(handles.TDEdit,'String',[]);
    end
    
end


% --- Executes on selection change in SaveDataType.
function SaveDataType_Callback(hObject, eventdata, handles)
% hObject    handle to SaveDataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SaveDataType contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SaveDataType

global AllResult NowResult NowResult_Label
AllResult=get(hObject,'String');
NowResult_Label=get(hObject,'Value');
NowResult=AllResult{NowResult_Label};


% --- Executes during object creation, after setting all properties.
function SaveDataType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SaveDataType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in DataSave.
function DataSave_Callback(hObject, eventdata, handles)
% hObject    handle to DataSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in DataReconstruction.

global NowResult NowResult_Label
if NowResult_Label>1
% [FileName,Path] = uiputfile('*.mat','Save Data',[NowResult,'.m']);
    if strcmp(NowResult,'MR_T1')
        global T1
        uisave('T1',NowResult);
    end
    if strcmp(NowResult,'CT')
        global CT
        uisave(NowResult,NowResult);
    end
    if strcmp(NowResult,'FsBrainmask')
        global fsmask
        uisave('fsmask',NowResult);
    end
    if strcmp(NowResult,'avg152T1')
        global avg152T1
        uisave(NowResult,NowResult);
    end
    if strcmp(NowResult,'regCT')
        global regCT
        uisave(NowResult,NowResult);
    end
    if strcmp(NowResult,'mask2regCT')
        global mask2regCT
        uisave(NowResult,NowResult);
    end
    if strcmp(NowResult,'T12avg')
        global T12avg
        uisave(NowResult,NowResult);
    end
    if strcmp(NowResult,'regCT2avg')
        global regCT2avg
        uisave(NowResult,NowResult);
    end
    if strcmp(NowResult,'T12avg_geomtform')
        global T12avg_geomtform
        uisave(NowResult,NowResult);
    end
    if strcmp(NowResult,'regCT2avg_geomtform')
        global regCT2avg_geomtform
        uisave(NowResult,NowResult);
    end
    if strcmp(NowResult,'e17EiContacts_sub')
        global e17EiContacts_sub
        uisave(NowResult,NowResult);
    end
    if strcmp(NowResult,'e17EiContacts_RASC')
        global e17EiContacts_RASC
        uisave(NowResult,NowResult);
    end
    if strcmp(NowResult,'ElecInfo_cell')
        global ElecInfo_cell
        uisave(NowResult,NowResult);
    end
else
    f=warndlg('Please choose data to save.','Warning');
end


% --- Executes on slider movement.
function AxialSlider_Callback(hObject, eventdata, handles)
% hObject    handle to AxialSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global NonReconJudge ReconstructionJudge RASTransformJudge RegistrationJudge NonChooseContact SegmentationJudge
global Images Images_recon IrAV nx1 nx2 ny2 xslice2 yslice2 zslice2
global Axial_WinMin Axial_WinMax AxialColormap AxialAlpha TDColormap TDAlpha

AV=get(hObject,'Value');
SV=get(handles.SagittalSlider,'Value');
CV=get(handles.CoronalSlider,'Value');
rAV=round(AV);
rSV=round(SV);
rCV=round(CV);
set(hObject,'Value',rAV);

if NonReconJudge==1
    IrAV=fliplr(squeeze(Images(nx1-rAV+1,:,:)));
    axes(handles.AxialPlane);
    imshow(IrAV,[Axial_WinMin Axial_WinMax]);
    colormap(handles.AxialPlane,AxialColormap);
    alpha(handles.AxialPlane,AxialAlpha);
end

if (ReconstructionJudge==1)||(RASTransformJudge==1)

    IrAV=fliplr(squeeze(Images_recon(nx2-rAV+1,:,:)));
    axes(handles.AxialPlane);
    imshow(IrAV,[Axial_WinMin Axial_WinMax]);
    colormap(handles.AxialPlane,AxialColormap);
    alpha(handles.AxialPlane,AxialAlpha);

    xslice1 = nx2-rAV+1; 
    yslice1 = ny2-rCV+1;
    zslice1 = rSV;

    global IrSampling nx3 ny3 nz3 ratio
    global RAS_IrSampling RAS_IrSampling_WorldLimits
    global XYZ_IrSampling XYZ_IrSampling_IntrinsicLimits

    xi=round(xslice1/ratio);
    yi=round(yslice1/ratio);
    zi=round(zslice1/ratio);
    xslice2 = xi+(xi==0)-(xi>nx3);
    yslice2 = yi+(yi==0)-(yi>ny3);
    zslice2 = zi+(zi==0)-(zi>nz3);
    axes(handles.TDPlane);
    CurrentView=get(handles.TDPlane,'View');

    if ReconstructionJudge==1
        SliceInRASC3(IrSampling,XYZ_IrSampling,1,xslice2,XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,XYZ_IrSampling,2,yslice2,XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,XYZ_IrSampling,3,zslice2,XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
        colormap(handles.TDPlane,TDColormap);
        alpha(handles.TDPlane,TDAlpha);
        axis equal;view(CurrentView);
        xlabel('X');ylabel('Y');zlabel('Z');
    end
    if RASTransformJudge==1
        SliceInRASC3(IrSampling,RAS_IrSampling,1,xslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,RAS_IrSampling,2,yslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,RAS_IrSampling,3,zslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        colormap(handles.TDPlane,TDColormap);
        alpha(handles.TDPlane,TDAlpha);
        axis equal;view(CurrentView);
        if NonChooseContact==1
            global ElecNum ElecLabel e17EiContacts_RASC NowElec_RASC
            axes(handles.TDPlane);
            for k=1:ElecNum
                if k~=ElecLabel
                    EiCrasc=e17EiContacts_RASC{k};
                    Scatter3InRASC5(EiCrasc,RAS_IrSampling_WorldLimits,9,'o',0.5,'b','b',0.75,0.75);hold on
                end
            end
            axes(handles.TDPlane);
            Scatter3InRASC5(NowElec_RASC,RAS_IrSampling_WorldLimits,9,'o',0.5,'g','g',0.75,0.75);hold on
            hold off;
        end
    end
    
    colormap(handles.TDPlane,TDColormap);
    alpha(handles.TDPlane,TDAlpha);

    set(handles.TDPlane,'XColor','w');
    set(handles.TDPlane,'YColor','w');
    set(handles.TDPlane,'ZColor','w');
    hold off
    axis equal;view(CurrentView);

end

if RegistrationJudge==1
    global fixedVolume regMoving nx4
    axes(handles.AxialPlane);
    imshowpair(fliplr(squeeze(fixedVolume(nx4-rAV+1,:,:))),fliplr(squeeze(regMoving(nx4-rAV+1,:,:))));
end

if SegmentationJudge==1
    global NowCon_sub
    set(hObject,'Value',NowCon_sub(1));
end


% --- Executes during object creation, after setting all properties.
function AxialSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AxialSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function SagittalSlider_Callback(hObject, eventdata, handles)
% hObject    handle to SagittalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


global NonReconJudge ReconstructionJudge RASTransformJudge RegistrationJudge NonChooseContact SegmentationJudge
global Images Images_recon IrSV nx2 ny2 xslice2 yslice2 zslice2
global Sagittal_WinMin Sagittal_WinMax SagittalColormap SagittalAlpha TDColormap TDAlpha

SV=get(hObject,'Value');
AV=get(handles.AxialSlider,'Value');
CV=get(handles.CoronalSlider,'Value');
rAV=round(AV);
rSV=round(SV);
rCV=round(CV);
set(hObject,'Value',rSV);

if NonReconJudge==1
    IrSV=squeeze(Images(:,:,rSV));
    axes(handles.SagittalPlane);
    imshow(IrSV,[Sagittal_WinMin Sagittal_WinMax]);
    colormap(handles.SagittalPlane,SagittalColormap);
    alpha(handles.SagittalPlane,SagittalAlpha);
end

if (ReconstructionJudge==1)||(RASTransformJudge==1)

    IrSV=squeeze(Images_recon(:,:,rSV));
    axes(handles.SagittalPlane);
    imshow(IrSV,[Sagittal_WinMin Sagittal_WinMax]);
    colormap(handles.SagittalPlane,SagittalColormap);
    alpha(handles.SagittalPlane,SagittalAlpha);

    xslice1 = nx2-rAV+1; 
    yslice1 = ny2-rCV+1;
    zslice1 = rSV;

    global IrSampling nx3 ny3 nz3 ratio
    global RAS_IrSampling RAS_IrSampling_WorldLimits
    global XYZ_IrSampling XYZ_IrSampling_IntrinsicLimits

    xi=round(xslice1/ratio);
    yi=round(yslice1/ratio);
    zi=round(zslice1/ratio);
    xslice2 = xi+(xi==0)-(xi>nx3);
    yslice2 = yi+(yi==0)-(yi>ny3);
    zslice2 = zi+(zi==0)-(zi>nz3);
    axes(handles.TDPlane);
    CurrentView=get(handles.TDPlane,'View');

    if ReconstructionJudge==1
        SliceInRASC3(IrSampling,XYZ_IrSampling,1,xslice2,XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,XYZ_IrSampling,2,yslice2,XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,XYZ_IrSampling,3,zslice2,XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
        colormap(handles.TDPlane,TDColormap);
        alpha(handles.TDPlane,TDAlpha);
        axis equal;view(CurrentView);
        xlabel('X');ylabel('Y');zlabel('Z');
    end
    if RASTransformJudge==1
        SliceInRASC3(IrSampling,RAS_IrSampling,1,xslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,RAS_IrSampling,2,yslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,RAS_IrSampling,3,zslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        colormap(handles.TDPlane,TDColormap);
        alpha(handles.TDPlane,TDAlpha);
        axis equal;view(CurrentView);
        if NonChooseContact==1
            global ElecNum ElecLabel e17EiContacts_RASC NowElec_RASC
            axes(handles.TDPlane);
            for k=1:ElecNum
                if k~=ElecLabel
                    EiCrasc=e17EiContacts_RASC{k};
                    Scatter3InRASC5(EiCrasc,RAS_IrSampling_WorldLimits,9,'o',0.5,'b','b',0.75,0.75);hold on
                end
            end
            axes(handles.TDPlane);
            Scatter3InRASC5(NowElec_RASC,RAS_IrSampling_WorldLimits,9,'o',0.5,'g','g',0.75,0.75);hold on
            hold off;
        end
    end
    
    colormap(handles.TDPlane,TDColormap);
    alpha(handles.TDPlane,TDAlpha);

    set(handles.TDPlane,'XColor','w');
    set(handles.TDPlane,'YColor','w');
    set(handles.TDPlane,'ZColor','w');
    hold off
    axis equal;view(CurrentView);

end

if RegistrationJudge==1
    global fixedVolume regMoving
    axes(handles.SagittalPlane);
    imshowpair(fixedVolume(:,:,rSV),regMoving(:,:,rSV));
end

if SegmentationJudge==1
    global NowCon_sub
    set(hObject,'Value',NowCon_sub(3));
end

% --- Executes during object creation, after setting all properties.
function SagittalSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SagittalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function CoronalSlider_Callback(hObject, eventdata, handles)
% hObject    handle to CoronalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global NonReconJudge ReconstructionJudge RASTransformJudge RegistrationJudge NonChooseContact SegmentationJudge
global Images Images_recon IrCV  ny1 nx2 ny2 xslice2 yslice2 zslice2
global Coronal_WinMin Coronal_WinMax CoronalColormap CoronalAlpha TDColormap TDAlpha

CV=get(hObject,'Value');
SV=get(handles.SagittalSlider,'Value');
AV=get(handles.AxialSlider,'Value');
rAV=round(AV);
rSV=round(SV);
rCV=round(CV);
set(hObject,'Value',rCV);

if NonReconJudge==1
    IrCV=fliplr(squeeze(Images(:,ny1-rCV+1,:)));
    axes(handles.CoronalPlane);
    imshow(IrCV,[Coronal_WinMin Coronal_WinMax]);
    colormap(handles.CoronalPlane,CoronalColormap);
    alpha(handles.CoronalPlane,CoronalAlpha);
end

if (ReconstructionJudge==1)||(RASTransformJudge==1)

    IrCV=fliplr(squeeze(Images_recon(:,ny2-rCV+1,:)));
    axes(handles.CoronalPlane);
    imshow(IrCV,[Coronal_WinMin Coronal_WinMax]);
    colormap(handles.CoronalPlane,CoronalColormap);
    alpha(handles.CoronalPlane,CoronalAlpha);

    xslice1 = nx2-rAV+1; 
    yslice1 = ny2-rCV+1;
    zslice1 = rSV;

    global IrSampling nx3 ny3 nz3 ratio
    global RAS_IrSampling RAS_IrSampling_WorldLimits
    global XYZ_IrSampling XYZ_IrSampling_IntrinsicLimits

    xi=round(xslice1/ratio);
    yi=round(yslice1/ratio);
    zi=round(zslice1/ratio);
    xslice2 = xi+(xi==0)-(xi>nx3);
    yslice2 = yi+(yi==0)-(yi>ny3);
    zslice2 = zi+(zi==0)-(zi>nz3);
    axes(handles.TDPlane);
    CurrentView=get(handles.TDPlane,'View');

    if ReconstructionJudge==1
        SliceInRASC3(IrSampling,XYZ_IrSampling,1,xslice2,XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,XYZ_IrSampling,2,yslice2,XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,XYZ_IrSampling,3,zslice2,XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
        colormap(handles.TDPlane,TDColormap);
        alpha(handles.TDPlane,TDAlpha);
        axis equal;view(CurrentView);
        xlabel('X');ylabel('Y');zlabel('Z');
    end
    if RASTransformJudge==1
        SliceInRASC3(IrSampling,RAS_IrSampling,1,xslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,RAS_IrSampling,2,yslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,RAS_IrSampling,3,zslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        colormap(handles.TDPlane,TDColormap);
        alpha(handles.TDPlane,TDAlpha);
        axis equal;view(CurrentView);
        if NonChooseContact==1
            global ElecNum ElecLabel e17EiContacts_RASC NowElec_RASC
            axes(handles.TDPlane);
            for k=1:ElecNum
                if k~=ElecLabel
                    EiCrasc=e17EiContacts_RASC{k};
                    Scatter3InRASC5(EiCrasc,RAS_IrSampling_WorldLimits,9,'o',0.5,'b','b',0.75,0.75);hold on
                end
            end
            axes(handles.TDPlane);
            Scatter3InRASC5(NowElec_RASC,RAS_IrSampling_WorldLimits,9,'o',0.5,'g','g',0.75,0.75);hold on
            hold off;
        end
    end
    
    colormap(handles.TDPlane,TDColormap);
    alpha(handles.TDPlane,TDAlpha);

    set(handles.TDPlane,'XColor','w');
    set(handles.TDPlane,'YColor','w');
    set(handles.TDPlane,'ZColor','w');
    hold off
    axis equal;view(CurrentView);

end

if RegistrationJudge==1
    global fixedVolume regMoving ny4
    axes(handles.CoronalPlane);
    imshowpair(fliplr(squeeze(fixedVolume(:,ny4-rCV+1,:))),fliplr(squeeze(regMoving(:,ny4-rCV+1,:))));
end

if SegmentationJudge==1
    global NowCon_sub
    set(hObject,'Value',NowCon_sub(2));
end

% --- Executes during object creation, after setting all properties.
function CoronalSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CoronalSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function AxialPlane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TDPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate TDPlane


% --- Executes during object creation, after setting all properties.
function SagittalPlane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SagittalPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate SagittalPlane


% --- Executes during object creation, after setting all properties.
function CoronalPlane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CoronalPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate CoronalPlane


% --- Executes during object creation, after setting all properties.
function TDPlane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TDPlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate TDPlane


function TDEdit_Callback(hObject, eventdata, handles)
% hObject    handle to TDEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TDEdit as text
%        str2double(get(hObject,'String')) returns contents of TDEdit as a double


% --- Executes during object creation, after setting all properties.
function TDEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TDEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Mouse.
function Mouse_Callback(hObject, eventdata, handles)
% hObject    handle to Mouse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotate3d off
zoom off

% --- Executes on button press in Rotate.
function Rotate_Callback(hObject, eventdata, handles)
% hObject    handle to Rotate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
rotate3d on


% --- Executes on button press in Zoom.
function Zoom_Callback(hObject, eventdata, handles)
% hObject    handle to Zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
zoom on


% --- Executes on selection change in ChoosePlane.
function ChoosePlane_Callback(hObject, eventdata, handles)
% hObject    handle to ChoosePlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChoosePlane contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChoosePlane

global AllPlaneName Plane
PlaneLabel=get(hObject,'Value')-1;
if PlaneLabel>0
    Plane=AllPlaneName{PlaneLabel};
    global Axial_WinMin Axial_WinMax AxialColormap_Label AxialAlpha
    global Sagittal_WinMin Sagittal_WinMax SagittalColormap_Label SagittalAlpha
    global Coronal_WinMin Coronal_WinMax CoronalColormap_Label CoronalAlpha
    global TDColormap_Label TDAlpha

    if strcmp(Plane,AllPlaneName{1})
        set(handles.WindowMin,'Value',Axial_WinMin);
        set(handles.WindowMax,'Value',Axial_WinMax);
        set(handles.ChooseColormap,'Value',AxialColormap_Label+1);
        set(handles.Trans,'String',AxialAlpha);
    end
    if strcmp(Plane,AllPlaneName{2})
        set(handles.WindowMin,'Value',Sagittal_WinMin);
        set(handles.WindowMax,'Value',Sagittal_WinMax);
        set(handles.ChooseColormap,'Value',SagittalColormap_Label+1);
        set(handles.Trans,'String',SagittalAlpha);
    end
    if strcmp(Plane,AllPlaneName{3})
        set(handles.WindowMin,'Value',Coronal_WinMin);
        set(handles.WindowMax,'Value',Coronal_WinMax);
        set(handles.ChooseColormap,'Value',CoronalColormap_Label+1);
        set(handles.Trans,'String',CoronalAlpha);
    end
    if strcmp(Plane,AllPlaneName{4})
        TD_WinMin=get(handles.WindowMin,'Min');
        set(handles.WindowMin,'Value',TD_WinMin);
        TD_WinMax=get(handles.WindowMax,'Max');
        set(handles.WindowMax,'Value',TD_WinMax);
        set(handles.ChooseColormap,'Value',TDColormap_Label+1);
        set(handles.Trans,'String',TDAlpha);
    end
end


% --- Executes during object creation, after setting all properties.
function ChoosePlane_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChoosePlane (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global AllPlaneName PlaneNum Plane
AllPlaneName={'Axial Plane','Sagittal Plane','Coronal Plane','3D View'};
PlaneNum=length(AllPlaneName);
Plane=AllPlaneName{1};


% --- Executes on slider movement.
function WindowMin_Callback(hObject, eventdata, handles)
% hObject    handle to WindowMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global WinMin WinMax 
WinMin=get(hObject,'Value');
WinMax=get(handles.WindowMax,'Value');

global AxialColormap SagittalColormap CoronalColormap
global AxialAlpha SagittalAlpha CoronalAlpha
global Axial_WinMin Sagittal_WinMin Coronal_WinMin
global Axial_WinMax Sagittal_WinMax Coronal_WinMax

global AllPlaneName Plane IrAV IrSV IrCV SegmentationJudge NowCon_sub nz2
if strcmp(Plane,AllPlaneName{1})
    Axial_WinMin=WinMin;
    Axial_WinMax=WinMax;
    axes(handles.AxialPlane);
    imshow(IrAV,[Axial_WinMin Axial_WinMax]);
    if SegmentationJudge==1
        hold on
        scatter(nz2-NowCon_sub(3)+1,NowCon_sub(2),16,'filled','o','Linewidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r',...
            'MarkerEdgeAlpha',0.75,'MarkerFaceAlpha',0.75);hold on;
        hold off
    end
    alpha(handles.AxialPlane,AxialAlpha);
    colormap(handles.AxialPlane,AxialColormap);
end
if strcmp(Plane,AllPlaneName{2})
    Sagittal_WinMin=WinMin;
    Sagittal_WinMax=WinMax;
    axes(handles.SagittalPlane);
    imshow(IrSV,[Sagittal_WinMin Sagittal_WinMax]);
    if SegmentationJudge==1
        hold on
        scatter(NowCon_sub(2),NowCon_sub(1),16,'filled','o','Linewidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r',...
            'MarkerEdgeAlpha',0.75,'MarkerFaceAlpha',0.75);hold on;
        hold off
    end
    alpha(handles.SagittalPlane,SagittalAlpha);
    colormap(handles.SagittalPlane,SagittalColormap);
end
if strcmp(Plane,AllPlaneName{3})
    Coronal_WinMin=WinMin;
    Coronal_WinMax=WinMax;
    axes(handles.CoronalPlane);
    imshow(IrCV,[Coronal_WinMin Coronal_WinMax]);
    if SegmentationJudge==1
        hold on
        scatter(nz2-NowCon_sub(3)+1,NowCon_sub(1),16,'filled','o','Linewidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r',...
            'MarkerEdgeAlpha',0.75,'MarkerFaceAlpha',0.75);hold on;
        hold off
    end
    alpha(handles.CoronalPlane,CoronalAlpha);
    colormap(handles.CoronalPlane,CoronalColormap);
end
if strcmp(Plane,AllPlaneName{4})
    TD_WinMin=get(handles.WindowMin,'Min');
    set(handles.WindowMin,'Value',TD_WinMin);
    TD_WinMax=get(handles.WindowMax,'Max');
    set(handles.WindowMax,'Value',TD_WinMax);
end


% --- Executes during object creation, after setting all properties.
function WindowMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function WindowMax_Callback(hObject, eventdata, handles)
% hObject    handle to WindowMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global WinMin WinMax
WinMin=get(handles.WindowMin,'Value');
WinMax=get(hObject,'Value');

global AxialColormap SagittalColormap CoronalColormap
global AxialAlpha SagittalAlpha CoronalAlpha
global Axial_WinMin Sagittal_WinMin Coronal_WinMin
global Axial_WinMax Sagittal_WinMax Coronal_WinMax

global AllPlaneName Plane IrAV IrSV IrCV SegmentationJudge NowCon_sub nz2
if strcmp(Plane,AllPlaneName{1})
    Axial_WinMin=WinMin;
    Axial_WinMax=WinMax;
    axes(handles.AxialPlane);
    imshow(IrAV,[Axial_WinMin Axial_WinMax]);
    if SegmentationJudge==1
        hold on
        scatter(nz2-NowCon_sub(3)+1,NowCon_sub(2),16,'filled','o','Linewidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r',...
            'MarkerEdgeAlpha',0.75,'MarkerFaceAlpha',0.75);hold on;
        hold off
    end
    alpha(handles.AxialPlane,AxialAlpha);
    colormap(handles.AxialPlane,AxialColormap);
end
if strcmp(Plane,AllPlaneName{2})
    Sagittal_WinMin=WinMin;
    Sagittal_WinMax=WinMax;
    axes(handles.SagittalPlane);
    imshow(IrSV,[Sagittal_WinMin Sagittal_WinMax]);
    if SegmentationJudge==1
        hold on
        scatter(NowCon_sub(2),NowCon_sub(1),16,'filled','o','Linewidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r',...
            'MarkerEdgeAlpha',0.75,'MarkerFaceAlpha',0.75);hold on;
        hold off
    end
    alpha(handles.SagittalPlane,SagittalAlpha);
    colormap(handles.SagittalPlane,SagittalColormap);
end
if strcmp(Plane,AllPlaneName{3})
    Coronal_WinMin=WinMin;
    Coronal_WinMax=WinMax;
    axes(handles.CoronalPlane);
    imshow(IrCV,[Coronal_WinMin Coronal_WinMax]);
    if SegmentationJudge==1
        hold on
        scatter(nz2-NowCon_sub(3)+1,NowCon_sub(1),16,'filled','o','Linewidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r',...
            'MarkerEdgeAlpha',0.75,'MarkerFaceAlpha',0.75);hold on;
        hold off
    end
    alpha(handles.CoronalPlane,CoronalAlpha);
    colormap(handles.CoronalPlane,CoronalColormap);
end
if strcmp(Plane,AllPlaneName{4})
    TD_WinMin=get(handles.WindowMin,'Min');
    set(handles.WindowMin,'Value',TD_WinMin);
    TD_WinMax=get(handles.WindowMax,'Max');
    set(handles.WindowMax,'Value',TD_WinMax);
end


% --- Executes during object creation, after setting all properties.
function WindowMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in ChooseColormap.
function ChooseColormap_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChooseColormap contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChooseColormap

global AllPlaneName Plane AllCmapValue Colormap
ColormapLabel=get(hObject,'Value')-1;
if ColormapLabel>0
    Colormap=AllCmapValue(:,(ColormapLabel*3-2):(ColormapLabel*3));
    global AxialColormap SagittalColormap CoronalColormap TDColormap
    global AxialColormap_Label SagittalColormap_Label CoronalColormap_Label TDColormap_Label
    if strcmp(Plane,AllPlaneName{1})
        AxialColormap=Colormap;
        AxialColormap_Label=ColormapLabel;
        colormap(handles.AxialPlane,Colormap);
    end
    if strcmp(Plane,AllPlaneName{2})
        SagittalColormap=Colormap;
        SagittalColormap_Label=ColormapLabel;
        colormap(handles.SagittalPlane,SagittalColormap);
    end
    if strcmp(Plane,AllPlaneName{3})
        CoronalColormap=Colormap;
        CoronalColormap_Label=ColormapLabel;
        colormap(handles.CoronalPlane,CoronalColormap);
    end
    if strcmp(Plane,AllPlaneName{4})
        TDColormap=Colormap;
        TDColormap_Label=ColormapLabel;
        colormap(handles.TDPlane,TDColormap);
    end
end


% --- Executes during object creation, after setting all properties.
function ChooseColormap_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChooseColormap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global AllCmapString AllCmapValue CmapNum Colormap
AllCmapString={'gray','parula','jet','hsv','hot','cool','spring','summer','autumn','winter','colorcube'};
AllCmapValue=[gray,parula,jet,hsv,hot,cool,spring,summer,autumn,winter,colorcube];
CmapNum=length(AllCmapString);
Colormap=gray;


function Trans_Callback(hObject, eventdata, handles)
% hObject    handle to Trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Trans as text
%        str2double(get(hObject,'String')) returns contents of Trans as a double

global AlphaValue AllPlaneName Plane
AlphaValue=str2double(get(hObject,'String'));

global AxialAlpha SagittalAlpha CoronalAlpha TDAlpha
if strcmp(Plane,AllPlaneName{1})
    AxialAlpha=AlphaValue;
    alpha(handles.AxialPlane,AxialAlpha);
end
if strcmp(Plane,AllPlaneName{2})
    SagittalAlpha=AlphaValue;
    alpha(handles.SagittalPlane,SagittalAlpha);
end
if strcmp(Plane,AllPlaneName{3})
    CoronalAlpha=AlphaValue;
    alpha(handles.CoronalPlane,CoronalAlpha);
end
if strcmp(Plane,AllPlaneName{4})
    TDAlpha=AlphaValue;
    alpha(handles.TDPlane,TDAlpha);
end


% --- Executes during object creation, after setting all properties.
function Trans_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Trans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global AlphaValue
AlphaValue=1;


% --- Executes on selection change in ChooseInterpolation.
function ChooseInterpolation_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseInterpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChooseInterpolation contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChooseInterpolation

global AllMethod IntpMethod IntpMethod_Label
IntpMethod_Label=get(hObject,'Value');
IntpMethod=AllMethod{IntpMethod_Label};


% --- Executes during object creation, after setting all properties.
function ChooseInterpolation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChooseInterpolation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global AllMethod
AllMethod=get(hObject,'String');


function DataReconstruction_Callback(hObject, eventdata, handles)
% hObject    handle to DataReconstruction (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global IntpMethod_Label
if IntpMethod_Label>1
    
    f = waitbar(0,'Data resampling. Please wait...');
    
    global NonReconJudge ReconstructionJudge RASTransformJudge RegistrationJudge NonChooseContact SegmentationJudge
    NonReconJudge=0;ReconstructionJudge=1;RASTransformJudge=0;RegistrationJudge=0;NonChooseContact=0;SegmentationJudge=0;

    set(handles.ChoosePlane,'Enable','on');
    set(handles.WindowMin,'Enable','on');
    set(handles.WindowMax,'Enable','on');
    set(handles.ChooseColormap,'Enable','on');
    set(handles.Trans,'Enable','on');

    global Images Images_PixelSpacing Images_SliceThickness Images_recon IntpMethod nx2 ny2 nz2
    Images_recon=VolumeReconstruction4(Images,Images_PixelSpacing,Images_SliceThickness,IntpMethod);
    waitbar(0.3,f);

    global FolderName T1_recon CT_recon fsmask_recon mask
    if strcmp(FolderName,'MR_T1')&&isempty(T1_recon)
        T1_recon=Images_recon;
    end
    if strcmp(FolderName,'CT')&&isempty(CT_recon)
        CT_recon=Images_recon;
    end
    if strcmp(FolderName,'FsBrainmask')&&isempty(fsmask_recon)
        fsmask_recon=Images_recon;
        mask=logical(fsmask_recon);
    end

    [nx2,ny2,nz2]=size(Images_recon);

    set(handles.AxialSlider,'Max',nx2);
    set(handles.SagittalSlider,'Max',nz2);
    set(handles.CoronalSlider,'Max',ny2);

    set(handles.AxialSlider,'Min',1);
    set(handles.SagittalSlider,'Min',1);
    set(handles.CoronalSlider,'Min',1);

    global graymin2 graymax2 SagCtr2 CorCtr2 AxiCtr2
    graymin2=min(min(min(Images_recon)));
    graymax2=max(max(max(Images_recon)));
    SagCtr2=Images_recon(:,:,round(nz2/2));
    CorCtr2=fliplr(squeeze(Images_recon(:,round(ny2/2),:)));
    AxiCtr2=fliplr(squeeze(Images_recon(round(nx2/2),:,:)));
    global IrAV IrSV IrCV
    IrAV=AxiCtr2;
    IrSV=SagCtr2;
    IrCV=CorCtr2;

    axes(handles.SagittalPlane);
    imshow(SagCtr2,[graymin2 graymax2]);
    axes(handles.CoronalPlane);
    imshow(CorCtr2,[graymin2 graymax2]);
    axes(handles.AxialPlane);
    imshow(AxiCtr2,[graymin2 graymax2]);

    set(handles.AxialSlider,'Value',round(nx2/2));
    set(handles.SagittalSlider,'Value',round(nz2/2));
    set(handles.CoronalSlider,'Value',round(ny2/2));

    set(handles.WindowMin,'Max',graymax2);
    set(handles.WindowMin,'Value',graymax2);
    set(handles.WindowMin,'Min',graymin2);
    set(handles.WindowMin,'Value',graymin2);

    set(handles.WindowMax,'Max',graymax2);
    set(handles.WindowMax,'Value',graymax2);
    set(handles.WindowMax,'Min',graymin2);

    waitbar(0.6,f);

    global Axial_WinMin Axial_WinMax AxialColormap AxialColormap_Label AxialAlpha
    global Sagittal_WinMin Sagittal_WinMax SagittalColormap SagittalColormap_Label SagittalAlpha
    global Coronal_WinMin Coronal_WinMax CoronalColormap CoronalColormap_Label CoronalAlpha
    global TDColormap TDColormap_Label TDAlpha

    Axial_WinMin=graymin2;Axial_WinMax=graymax2;AxialColormap=gray;AxialColormap_Label=1;AxialAlpha=1;
    Sagittal_WinMin=graymin2;Sagittal_WinMax=graymax2;SagittalColormap=gray;SagittalColormap_Label=1;SagittalAlpha=1;
    Coronal_WinMin=graymin2;Coronal_WinMax=graymax2;CoronalColormap=gray;CoronalColormap_Label=1;CoronalAlpha=1;
    TDColormap=gray;TDColormap_Label=1;TDAlpha=1;

    global Colormap AlphaValue Plane
    Plane='Axial Plane';Colormap=gray;AlphaValue=1;
    set(handles.ChoosePlane,'Value',1);
    set(handles.ChooseColormap,'Value',1);
    set(handles.Trans,'String',[]);

    global ratio IrSampling XYZ_IrSampling XYZ_IrSampling_IntrinsicLimits
    set(handles.TDEdit,'String','3D slices in native space');
    IrSampling=Images_recon(1:ratio:nx2,1:ratio:ny2,1:ratio:nz2);

    global Fname1
    if strcmp(Fname1,'avg152T1.nii')
        IrSampling=RotateAvg(IrSampling);
    end

    global nx3 ny3 nz3
    [nx3,ny3,nz3]=size(IrSampling);
    XYZ_IrSampling=zeros(4,ny3,nx3,nz3); %定义四维矩阵
    for x=1:ny3 %x为列索引,构建所有体素的X坐标（X坐标轴方向上是列的索引）
        XYZ_IrSampling(1,x,:,:)=(x-1)*ratio+1; 
    end
    for y=1:nx3 %y为行索引,构建所有体素的Y坐标（Y坐标轴方向上是行的索引）
        XYZ_IrSampling(2,:,y,:)=(y-1)*ratio+1; 
    end
    for z=1:nz3 %z为层索引,构建所有体素的Z坐标（Z坐标轴方向上是层的索引）
        XYZ_IrSampling(3,:,:,z)=(z-1)*ratio+1; 
    end
    XYZ_IrSampling(4,:,:,:)=1; %最后一维恒为1
    XYZ_IrSampling_IntrinsicLimits=[1,(ny3-1)*ratio+1,1,(nx3-1)*ratio+1,1,(nz3-1)*ratio+1];

    axes(handles.TDPlane);
    SliceInRASC3(IrSampling,XYZ_IrSampling,1,round(nx3/2),XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
    SliceInRASC3(IrSampling,XYZ_IrSampling,2,round(ny3/2),XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
    SliceInRASC3(IrSampling,XYZ_IrSampling,3,round(nz3/2),XYZ_IrSampling_IntrinsicLimits,handles.TDPlane);hold on;
    waitbar(1,f);

    xlabel('X');ylabel('Y');zlabel('Z');
    set(handles.TDPlane,'XColor','w');
    set(handles.TDPlane,'YColor','w');
    set(handles.TDPlane,'ZColor','w');
    hold off
    view(-5,-65);
    close(f);

else
    f=warndlg('Please choose the method of interpolation first.','Warning');    
end


% --- Executes on selection change in ChooseCoordinate.
function ChooseCoordinate_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseCoordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChooseCoordinate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChooseCoordinate

global AllCoordinate NowCoordinate NowCoordinate_Label
NowCoordinate_Label=get(hObject,'Value');
NowCoordinate=AllCoordinate{NowCoordinate_Label};


% --- Executes during object creation, after setting all properties.
function ChooseCoordinate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChooseCoordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global AllCoordinate
AllCoordinate=get(hObject,'String');


% --- Executes on button press in CoordinateTransformation.
function CoordinateTransformation_Callback(hObject, eventdata, handles)
% hObject    handle to CoordinateTransformation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global NowCoordinate_Label
if NowCoordinate_Label>1
    
    global NonReconJudge ReconstructionJudge RASTransformJudge RegistrationJudge NonChooseContact SegmentationJudge
    NonReconJudge=0;ReconstructionJudge=0;RASTransformJudge=1;RegistrationJudge=0;NonChooseContact=0;SegmentationJudge=0;

    set(handles.ChoosePlane,'Enable','on');
    set(handles.WindowMin,'Enable','on');
    set(handles.WindowMax,'Enable','on');
    set(handles.ChooseColormap,'Enable','on');
    set(handles.Trans,'Enable','on');

    global graymin2 graymax2 SagCtr2 CorCtr2 AxiCtr2 nx2 ny2 nz2
    axes(handles.SagittalPlane);
    imshow(SagCtr2,[graymin2 graymax2]);
    axes(handles.CoronalPlane);
    imshow(CorCtr2,[graymin2 graymax2]);
    axes(handles.AxialPlane);
    imshow(AxiCtr2,[graymin2 graymax2]);

    set(handles.AxialSlider,'Value',round(nx2/2));
    set(handles.SagittalSlider,'Value',round(nz2/2));
    set(handles.CoronalSlider,'Value',round(ny2/2));

    set(handles.WindowMin,'Max',graymax2);
    set(handles.WindowMin,'Value',graymax2);
    set(handles.WindowMin,'Min',graymin2);
    set(handles.WindowMin,'Value',graymin2);

    set(handles.WindowMax,'Max',graymax2);
    set(handles.WindowMax,'Value',graymax2);
    set(handles.WindowMax,'Min',graymin2);

    global Axial_WinMin Axial_WinMax AxialColormap AxialColormap_Label AxialAlpha
    global Sagittal_WinMin Sagittal_WinMax SagittalColormap SagittalColormap_Label SagittalAlpha
    global Coronal_WinMin Coronal_WinMax CoronalColormap CoronalColormap_Label CoronalAlpha
    global TDColormap TDColormap_Label TDAlpha

    Axial_WinMin=graymin2;Axial_WinMax=graymax2;AxialColormap=gray;AxialColormap_Label=1;AxialAlpha=1;
    Sagittal_WinMin=graymin2;Sagittal_WinMax=graymax2;SagittalColormap=gray;SagittalColormap_Label=1;SagittalAlpha=1;
    Coronal_WinMin=graymin2;Coronal_WinMax=graymax2;CoronalColormap=gray;CoronalColormap_Label=1;CoronalAlpha=1;
    TDColormap=gray;TDColormap_Label=1;TDAlpha=1;

    global Colormap AlphaValue Plane
    Plane='Axial Plane';Colormap=gray;AlphaValue=1;
    set(handles.ChoosePlane,'Value',1);
    set(handles.ChooseColormap,'Value',1);
    set(handles.Trans,'String',[]);

    global IrSampling ratio Images_IPP1 Images_IOP IrS_PixelSpacing IrS_SliceThickness nx3 ny3 nz3
    IrS_PixelSpacing=[0.6;0.6].*ratio;
    IrS_SliceThickness=0.6*ratio;
    coordinate_judge=1;

    global M_IrSampling RAS_IrSampling RAS_IrSampling_WorldLimits AllCoordinate NowCoordinate
    if strcmp(NowCoordinate,AllCoordinate{2})
        [M_IrSampling,RAS_IrSampling]=transform2LPH(IrSampling,IrS_PixelSpacing,IrS_SliceThickness,Images_IOP,Images_IPP1,coordinate_judge);
        set(handles.TDEdit,'String','3D slices in LPH coordinate system');
    end
    if strcmp(NowCoordinate,AllCoordinate{3})
        [M_IrSampling,RAS_IrSampling]=transform2RAS(IrSampling,IrS_PixelSpacing,IrS_SliceThickness,Images_IOP,Images_IPP1,coordinate_judge);
        set(handles.TDEdit,'String','3D slices in RAS coordinate system');
    end
    RAS_IrSampling_rmin=min(min(min(RAS_IrSampling(1,:,:,:))));
    RAS_IrSampling_rmax=max(max(max(RAS_IrSampling(1,:,:,:))));
    RAS_IrSampling_amin=min(min(min(RAS_IrSampling(2,:,:,:))));
    RAS_IrSampling_amax=max(max(max(RAS_IrSampling(2,:,:,:))));
    RAS_IrSampling_smin=min(min(min(RAS_IrSampling(3,:,:,:))));
    RAS_IrSampling_smax=max(max(max(RAS_IrSampling(3,:,:,:))));
    RAS_IrSampling_WorldLimits=[RAS_IrSampling_rmin RAS_IrSampling_rmax RAS_IrSampling_amin RAS_IrSampling_amax RAS_IrSampling_smin RAS_IrSampling_smax];

    axes(handles.TDPlane);
    SliceInRASC3(IrSampling,RAS_IrSampling,1,round(nx3/2),RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
    SliceInRASC3(IrSampling,RAS_IrSampling,2,round(ny3/2),RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
    SliceInRASC3(IrSampling,RAS_IrSampling,3,round(nz3/2),RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
    if strcmp(NowCoordinate,AllCoordinate{2})
        xlabel('L');ylabel('P');zlabel('H');
    end
    set(handles.TDPlane,'XColor','w');
    set(handles.TDPlane,'YColor','w');
    set(handles.TDPlane,'ZColor','w');
    hold off

else
    f=warndlg('Please choose coordinate before transformation.','Warning');
end


% --- Executes on selection change in ChooseFixed.
function ChooseFixed_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseFixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChooseFixed contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChooseFixed

global AllFixed NowFixed NowFixed_Label
NowFixed_Label=get(hObject,'Value');
NowFixed=AllFixed{NowFixed_Label};


% --- Executes during object creation, after setting all properties.
function ChooseFixed_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChooseFixed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global AllFixed NowFixed_Label
AllFixed=get(hObject,'String');
NowFixed_Label=1;


% --- Executes on selection change in ChooseMoving.
function ChooseMoving_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseMoving (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChooseMoving contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChooseMoving

global AllMoving NowMoving NowMoving_Label
NowMoving_Label=get(hObject,'Value');
NowMoving=AllMoving{NowMoving_Label};


% --- Executes during object creation, after setting all properties.
function ChooseMoving_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChooseMoving (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

global AllMoving NowMoving_Label
AllMoving=get(hObject,'String');
NowMoving_Label=1;


function InitialRadius_Callback(hObject, eventdata, handles)
% hObject    handle to InitialRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of InitialRadius as text
%        str2double(get(hObject,'String')) returns contents of InitialRadius as a double

global InitialRadius
InitialRadius=get(hObject,'String');
InitialRadius=str2double(InitialRadius);


% --- Executes during object creation, after setting all properties.
function InitialRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to InitialRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MaxInterations_Callback(hObject, eventdata, handles)
% hObject    handle to MaxInterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxInterations as text
%        str2double(get(hObject,'String')) returns contents of MaxInterations as a double

global MaxInterations
MaxInterations=get(hObject,'String');
MaxInterations=str2double(MaxInterations);


% --- Executes during object creation, after setting all properties.
function MaxInterations_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxInterations (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Registration.
function Registration_Callback(hObject, eventdata, handles)
% hObject    handle to Registration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global NowFixed NowMoving
if (strcmp(NowFixed,'MR_T1')&&strcmp(NowMoving,'CT'))||(strcmp(NowFixed,'regCT')&&strcmp(NowMoving,'FsBrainmask'))||...
        (strcmp(NowFixed,'avg152T1')&&strcmp(NowMoving,'MR_T1'))||(strcmp(NowFixed,'avg152T1')&&strcmp(NowMoving,'regCT'))
    
    global NonReconJudge ReconstructionJudge RASTransformJudge RegistrationJudge NonChooseContact SegmentationJudge
    NonReconJudge=0;ReconstructionJudge=0;RASTransformJudge=0;RegistrationJudge=1;NonChooseContact=0;SegmentationJudge=0;
    f=waitbar(0,'Image registration is in progress. Please wait...');

    set(handles.ChoosePlane,'Enable','off');
    set(handles.WindowMin,'Enable','off');
    set(handles.WindowMax,'Enable','off');
    set(handles.ChooseColormap,'Enable','off');
    set(handles.Trans,'Enable','off');
    set(handles.ChooseInterpolation,'Enable','off');
    set(handles.DataReconstruction,'Enable','off');
    set(handles.ChooseCoordinate,'Enable','off');
    set(handles.CoordinateTransformation,'Enable','off');

    global T1_path CT_path fsmask_path T1_recon CT_recon 
    global regCT mask2regCT mask InitialRadius MaxInterations
    global fixedVolume movingVolume regMoving

    if (strcmp(NowFixed,'MR_T1')&&strcmp(NowMoving,'CT')&&isempty(regCT))||(strcmp(NowFixed,'regCT')&&strcmp(NowMoving,'FsBrainmask')&&isempty(mask2regCT))
        
        if strcmp(NowFixed,'MR_T1')&&strcmp(NowMoving,'CT')&&isempty(regCT)
            MRI_picstr=dir([T1_path,'*.dcm']);
            MRI_info=dicominfo([T1_path,MRI_picstr(1).name]);
            CT_picstr=dir([CT_path,'*.dcm']);
            CT_info=dicominfo([CT_path,CT_picstr(1).name]);

            fixedHeader=MRI_info;
            movingHeader=CT_info;

            fixedVolume=single(T1_recon);
            movingVolume=single(CT_recon);
        end

        if strcmp(NowFixed,'regCT')&&strcmp(NowMoving,'FsBrainmask')&&isempty(mask2regCT)
            if isempty(regCT)
                f2=warndlg('regCT was not found. Please co-register MR_T1 and CT first.','Warning');
            else
                regCT_picstr=dir([T1_path,'*.dcm']);
                regCT_info=dicominfo([T1_path,regCT_picstr(1).name]);
                fsmask_picstr=dir([fsmask_path,'*.dcm']);
                fsmask_info=dicominfo([fsmask_path,fsmask_picstr(1).name]);

                fixedHeader=regCT_info;
                movingHeader=fsmask_info;

                fixedVolume=single(regCT);
                movingVolume=single(mask);
            end
        end
        waitbar(0.2,f);

        ratio=3;

        fixedHeader.PixelSpacing(1)=0.6;
        fixedHeader.PixelSpacing(2)=0.6;
        fixedHeader.SliceThickness=0.6;

        movingHeader.PixelSpacing(1)=0.6;
        movingHeader.PixelSpacing(2)=0.6;
        movingHeader.SliceThickness=0.6;

        [movingRegisteredVolume,optimizer_small,geomtform_small]=Registration(fixedHeader,fixedVolume,...
                  movingHeader,movingVolume,ratio,InitialRadius,MaxInterations,f);

        if strcmp(NowFixed,'MR_T1')&&strcmp(NowMoving,'CT')
            global regCT_optimizer_small regCT_geomtform_small
            regCT=movingRegisteredVolume;
            regCT_optimizer_small=optimizer_small;
            regCT_geomtform_small=geomtform_small;
            regMoving=regCT;

            SaveDataStr=get(handles.SaveDataType,'String');
            L=size(SaveDataStr,1);
            SaveDataStr2=cell(L+1,1);
            for i=1:L
                SaveDataStr2{i,1}=SaveDataStr{i,1};
            end
            SaveDataStr2{L+1,1}='regCT';
            set(handles.SaveDataType,'String',SaveDataStr2);
        end

        if strcmp(NowFixed,'regCT')&&strcmp(NowMoving,'FsBrainmask')
            global mask2regCT_optimizer_small mask2regCT_geomtform_small
            mask2regCT=movingRegisteredVolume;
            mask2regCT_optimizer_small=optimizer_small;
            mask2regCT_geomtform_small=geomtform_small;
            regMoving=mask2regCT;

            SaveDataStr=get(handles.SaveDataType,'String');
            L=size(SaveDataStr,1);
            SaveDataStr2=cell(L+1,1);
            for i=1:L
                SaveDataStr2{i,1}=SaveDataStr{i,1};
            end
            SaveDataStr2{L+1,1}='mask2regCT';
            set(handles.SaveDataType,'String',SaveDataStr2);
        end
        waitbar(0.8,f);

    end

    if strcmp(NowFixed,'MR_T1')&&strcmp(NowMoving,'CT')&&(~isempty(regCT))
        fixedVolume=single(T1_recon);
        regMoving=regCT;
        waitbar(0.8,f);
    end

    if strcmp(NowFixed,'regCT')&&strcmp(NowMoving,'FsBrainmask')&&(~isempty(mask2regCT))
        fixedVolume=single(regCT);
        regMoving=mask2regCT;
        waitbar(0.8,f);
    end

    global avg152T1 T12avg regCT2avg
    global T12avg_optimizer T12avg_geomtform
    global regCT2avg_optimizer regCT2avg_geomtform
    if (strcmp(NowFixed,'avg152T1')&&strcmp(NowMoving,'MR_T1')&&isempty(T12avg))||(strcmp(NowFixed,'avg152T1')&&strcmp(NowMoving,'regCT')&&isempty(regCT2avg))

        global avg152T1_IPP1 avg152T1_IOP avg152T1_PixelSpacing avg152T1_SliceThickness coordinate_judge M_avg152T1 RAS_avg152T1 RAS_avg152T1_WorldLimits
        global T1_PixelSpacing T1_SliceThickness T1_IOP T1_IPP1 M_T1_recon RAS_T1_recon RAS_T1_recon_WorldLimits
        coordinate_judge=1;

        if isempty(RAS_avg152T1_WorldLimits)
            [M_avg152T1,RAS_avg152T1]=transform2RAS(avg152T1,avg152T1_PixelSpacing,avg152T1_SliceThickness,avg152T1_IOP,avg152T1_IPP1,coordinate_judge);
            RAS_avg152T1_rmin=min(min(min(RAS_avg152T1(1,:,:,:))));
            RAS_avg152T1_rmax=max(max(max(RAS_avg152T1(1,:,:,:))));
            RAS_avg152T1_amin=min(min(min(RAS_avg152T1(2,:,:,:))));
            RAS_avg152T1_amax=max(max(max(RAS_avg152T1(2,:,:,:))));
            RAS_avg152T1_smin=min(min(min(RAS_avg152T1(3,:,:,:))));
            RAS_avg152T1_smax=max(max(max(RAS_avg152T1(3,:,:,:))));
            RAS_avg152T1_WorldLimits=[RAS_avg152T1_rmin RAS_avg152T1_rmax RAS_avg152T1_amin RAS_avg152T1_amax RAS_avg152T1_smin RAS_avg152T1_smax];
        end
        waitbar(0.1,f);

        if isempty(RAS_T1_recon_WorldLimits)
            [M_T1_recon,RAS_T1_recon]=transform2RAS(T1_recon,[0.6;0.6],0.6,T1_IOP,T1_IPP1,coordinate_judge);
            RAS_T1_recon_rmin=min(min(min(RAS_T1_recon(1,:,:,:))));
            RAS_T1_recon_rmax=max(max(max(RAS_T1_recon(1,:,:,:))));
            RAS_T1_recon_amin=min(min(min(RAS_T1_recon(2,:,:,:))));
            RAS_T1_recon_amax=max(max(max(RAS_T1_recon(2,:,:,:))));
            RAS_T1_recon_smin=min(min(min(RAS_T1_recon(3,:,:,:))));
            RAS_T1_recon_smax=max(max(max(RAS_T1_recon(3,:,:,:))));
            RAS_T1_recon_WorldLimits=[RAS_T1_recon_rmin RAS_T1_recon_rmax RAS_T1_recon_amin RAS_T1_recon_amax RAS_T1_recon_smin RAS_T1_recon_smax];
        end
        waitbar(0.3,f);

        Rmin=RAS_avg152T1_WorldLimits(1)-avg152T1_PixelSpacing(1)/2;
        Rmax=RAS_avg152T1_WorldLimits(2)+avg152T1_PixelSpacing(1)/2;
        Amin=RAS_avg152T1_WorldLimits(3)-avg152T1_PixelSpacing(2)/2;
        Amax=RAS_avg152T1_WorldLimits(4)+avg152T1_PixelSpacing(2)/2;
        Smin=RAS_avg152T1_WorldLimits(5)-avg152T1_SliceThickness/2;
        Smax=RAS_avg152T1_WorldLimits(6)+avg152T1_SliceThickness/2;

        rmin=RAS_T1_recon_WorldLimits(1)-T1_SliceThickness/2;
        rmax=RAS_T1_recon_WorldLimits(2)+T1_SliceThickness/2;
        amin=RAS_T1_recon_WorldLimits(3)-T1_PixelSpacing(2)/2;
        amax=RAS_T1_recon_WorldLimits(4)+T1_PixelSpacing(2)/2;
        smin=RAS_T1_recon_WorldLimits(5)-T1_PixelSpacing(1)/2;
        smax=RAS_T1_recon_WorldLimits(6)+T1_PixelSpacing(1)/2;

        global origVolume origVolume_ci3
        if strcmp(NowFixed,'avg152T1')&&strcmp(NowMoving,'MR_T1')&&isempty(T12avg)
            origVolume=T1_recon;
        end
        if strcmp(NowFixed,'avg152T1')&&strcmp(NowMoving,'regCT')&&isempty(regCT2avg)
            origVolume=regCT;
        end

        [nx5,ny5,nz5]=size(origVolume);
        origVolume_ci=origVolume;
        for i=1:nz5
            origVolume_ci(:,:,i)=fliplr(origVolume(:,:,i)); %左右翻转
        end
        origVolume_ci2=origVolume_ci;
        for i=1:nz5
            origVolume_ci2(:,:,i)=flipud(origVolume_ci(:,:,i)); %上下翻转
        end
        origVolume_ci3=permute(origVolume_ci2,[3,2,1]); %改变索引顺序
        [nx6,ny6,nz6]=size(origVolume_ci3);
        for i=1:nz6
            origVolume_ci3(:,:,i)=flipud(origVolume_ci3(:,:,i)); %上下翻转
        end

        Rfixed=imref3d(size(avg152T1),[Amin Amax],[Rmin Rmax],[Smin Smax]); %见文档《imref3d》
        Rmoving=imref3d(size(origVolume_ci3),[amin amax],[rmin rmax],[smin smax]); %见文档《imref3d》

        fixedVolume=single(avg152T1);
        movingVolume=single(origVolume_ci3);

        [optimizer,metric]=imregconfig('multimodal');
        optimizer.InitialRadius=InitialRadius;
        optimizer.MaximumIterations=MaxInterations;
        waitbar(0.4,f);

        if strcmp(NowFixed,'avg152T1')&&strcmp(NowMoving,'MR_T1')&&isempty(T12avg)
            if isempty(regCT2avg_geomtform)
                geomtform=imregtform(movingVolume,Rmoving,fixedVolume,Rfixed,'affine',optimizer,metric);
                waitbar(0.7,f);
        %         geomtform.T
            else
                geomtform=regCT2avg_geomtform;
            end
            movingRegisteredVolume=imwarp(movingVolume,...
                Rmoving,geomtform,'bicubic','OutputView',Rfixed);
            waitbar(0.9,f);

            T12avg=movingRegisteredVolume;
            T12avg_optimizer=optimizer;
            T12avg_geomtform=geomtform;
            regMoving=T12avg;

            SaveDataStr=get(handles.SaveDataType,'String');
            L=size(SaveDataStr,1);
            SaveDataStr2=cell(L+2,1);
            for i=1:L
                SaveDataStr2{i,1}=SaveDataStr{i,1};
            end
            SaveDataStr2{L+1,1}='T12avg';
            SaveDataStr2{L+2,1}='T12avg_geomtform';
            set(handles.SaveDataType,'String',SaveDataStr2);
        end

        if strcmp(NowFixed,'avg152T1')&&strcmp(NowMoving,'regCT')&&isempty(regCT2avg)
            if isempty(T12avg_geomtform)
                geomtform=imregtform(movingVolume,Rmoving,fixedVolume,Rfixed,'affine',optimizer,metric);
                waitbar(0.7,f);
        %         geomtform.T
            else
                geomtform=T12avg_geomtform;
            end
            movingRegisteredVolume=imwarp(movingVolume,...
                Rmoving,geomtform,'bicubic','OutputView',Rfixed);
            waitbar(0.9,f);

            regCT2avg=movingRegisteredVolume;
            regCT2avg_optimizer=optimizer;
            regCT2avg_geomtform=geomtform;
            regMoving=regCT2avg;

            SaveDataStr=get(handles.SaveDataType,'String');
            L=size(SaveDataStr,1);
            SaveDataStr2=cell(L+2,1);
            for i=1:L
                SaveDataStr2{i,1}=SaveDataStr{i,1};
            end
            SaveDataStr2{L+1,1}='regCT2avg';
            SaveDataStr2{L+2,1}='regCT2avg_geomtform';
            set(handles.SaveDataType,'String',SaveDataStr2);
        end

    end

    if strcmp(NowFixed,'avg152T1')&&strcmp(NowMoving,'MR_T1')&&(~isempty(T12avg))
        fixedVolume=single(avg152T1);
        regMoving=T12avg;
        waitbar(0.8,f);
    end

    if strcmp(NowFixed,'avg152T1')&&strcmp(NowMoving,'regCT')&&(~isempty(regCT2avg))
        fixedVolume=single(avg152T1);
        regMoving=regCT2avg;
        waitbar(0.8,f);
    end

    if strcmp(NowFixed,'avg152T1')
        fixedVolume=InvRotateAvg(fixedVolume);
        regMoving=InvRotateAvg(regMoving);
    end
    
    global nx4 ny4 nz4
    [nx4,ny4,nz4]=size(regMoving);
    graymin3=min(min(min(regMoving)));
    graymax3=max(max(max(regMoving)));

    set(handles.AxialSlider,'Max',nx4);
    set(handles.SagittalSlider,'Max',nz4);
    set(handles.CoronalSlider,'Max',ny4);

    set(handles.AxialSlider,'Min',1);
    set(handles.SagittalSlider,'Min',1);
    set(handles.CoronalSlider,'Min',1);

    set(handles.AxialSlider,'Value',round(nx4/2));
    set(handles.SagittalSlider,'Value',round(nz4/2));
    set(handles.CoronalSlider,'Value',round(ny4/2));

    set(handles.WindowMin,'Max',graymax3);
    set(handles.WindowMin,'Value',graymax3);
    set(handles.WindowMin,'Min',graymin3);
    set(handles.WindowMin,'Value',graymin3);

    set(handles.WindowMax,'Max',graymax3);
    set(handles.WindowMax,'Value',graymax3);
    set(handles.WindowMax,'Min',graymin3);

    global Colormap AlphaValue Plane
    Plane='Axial Plane';Colormap=gray;AlphaValue=1;
    set(handles.ChoosePlane,'Value',1);
    set(handles.ChooseColormap,'Value',1);
    set(handles.Trans,'String',[]);
    waitbar(1,f);

    axes(handles.AxialPlane);
    imshowpair(fliplr(squeeze(fixedVolume(round(nx4/2),:,:))),fliplr(squeeze(regMoving(round(nx4/2),:,:))));
    axes(handles.SagittalPlane);
    imshowpair(fixedVolume(:,:,round(nz4/2)),regMoving(:,:,round(nz4/2)));
    axes(handles.CoronalPlane);
    imshowpair(fliplr(squeeze(fixedVolume(:,round(ny4/2),:))),fliplr(squeeze(regMoving(:,round(ny4/2),:))));

    global StartImage
    axes(handles.TDPlane);
    imshow(StartImage,[0 1]);
    set(handles.TDEdit,'String',[]);
    close(f);
    
else
    f=warndlg('The selection of the fixed volume or moving volume is wrong, please choose the correct one.','Warning');
end


function ElecNum_Callback(hObject, eventdata, handles)
% hObject    handle to ElecNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ElecNum as text
%        str2double(get(hObject,'String')) returns contents of ElecNum as a double

global ElecNum
ElecNum=get(hObject,'String');
ElecNum=str2double(ElecNum);


% --- Executes during object creation, after setting all properties.
function ElecNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ElecNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB 
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ContactSegmentation_SYSU.
function ContactsSegmentation_Callback(hObject, eventdata, handles)
% hObject    handle to ContactSegmentation_SYSU (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global NonReconJudge ReconstructionJudge RASTransformJudge RegistrationJudge NonChooseContact SegmentationJudge
NonReconJudge=0;ReconstructionJudge=0;RASTransformJudge=1;RegistrationJudge=0;NonChooseContact=0;SegmentationJudge=0;

global T12avg_geomtform ElecNum T1_IOP T1_IPP1
if isempty(ElecNum)||isnan(ElecNum)||ElecNum==0
    f1=warndlg('ElecNum has not been input. Please imput it first.','Warning');
else
    
    if isempty(T12avg_geomtform)
        f2=warndlg('T12avg_geomtform has not been imported. Please import it first.','Warning');
    end

    if (~isempty(T12avg_geomtform))

        f3=warndlg('If you have already done the contact segmentation, you can import e17EiContacts_sub, e17EiContacts_RASC and ElecInfo_cell directly to view the contacts.','Prompt');

        set(handles.ChoosePlane,'Enable','on');
        set(handles.WindowMin,'Enable','on');
        set(handles.WindowMax,'Enable','on');
        set(handles.ChooseColormap,'Enable','on');
        set(handles.Trans,'Enable','on');
        set(handles.ChooseInterpolation,'Enable','off');
        set(handles.DataReconstruction,'Enable','off');
        set(handles.ChooseCoordinate,'Enable','off');
        set(handles.CoordinateTransformation,'Enable','off');

        f=waitbar(0,'Begin to segment electrode contacts. It may take a long time, please be patient...');

        global e17EiContacts_sub e17EiContacts_RASC ElecInfo_cell
        if isempty(e17EiContacts_sub)||isempty(e17EiContacts_RASC)||isempty(ElecInfo_cell)

            global M_avg152T1
            if isempty(M_avg152T1)
                global avg152T1 avg152T1_IPP1 avg152T1_IOP avg152T1_PixelSpacing avg152T1_SliceThickness
                [M_avg152T1,~]=transform2RAS(avg152T1,avg152T1_PixelSpacing,avg152T1_SliceThickness,avg152T1_IOP,avg152T1_IPP1,1);
            end

            global regCT mask2regCT
            global brainmask brainmask2 brainmask3 brainmask_edge brainmask_edge2 brainmask_edge3 T0 Te ...
                electrode1 electrode1B electrode2 electrode3 electrode4 imclose_judge imclose_r ...
                VIclose_nrec electrode5 electrode5_S
            [brainmask,brainmask2,brainmask3,brainmask_edge,brainmask_edge2,brainmask_edge3,T0,Te,...
                electrode1,electrode1B,electrode2,electrode3,electrode4,imclose_judge,imclose_r,...
                VIclose_nrec,electrode5,electrode5_S]=ContactSeg1(ElecNum,regCT,mask2regCT,f);

    %         path='D:\电脑资料2019.03.17开始\毕设\电极定位\subjects\wu chuyin\step2\';
    %         brainmask=load([path,'brainmask.mat']);brainmask=struct2cell(brainmask);brainmask=brainmask{1,1};
    %         brainmask2=load([path,'brainmask2.mat']);brainmask2=struct2cell(brainmask2);brainmask2=brainmask2{1,1};
    %         brainmask3=load([path,'brainmask3.mat']);brainmask3=struct2cell(brainmask3);brainmask3=brainmask3{1,1};
    %         brainmask_edge=load([path,'brainmask_edge.mat']);brainmask_edge=struct2cell(brainmask_edge);brainmask_edge=brainmask_edge{1,1};
    %         brainmask_edge2=load([path,'brainmask_edge2.mat']);brainmask_edge2=struct2cell(brainmask_edge2);brainmask_edge2=brainmask_edge2{1,1};
    %         brainmask_edge3=load([path,'brainmask_edge3.mat']);brainmask_edge3=struct2cell(brainmask_edge3);brainmask_edge3=brainmask_edge3{1,1};
    %         T0=load([path,'T0.mat']);T0=struct2cell(T0);T0=T0{1,1};
    %         Te=load([path,'Te.mat']);Te=struct2cell(Te);Te=Te{1,1};
    %         electrode1=load([path,'electrode1.mat']);electrode1=struct2cell(electrode1);electrode1=electrode1{1,1};
    %         electrode1B=load([path,'electrode1B.mat']);electrode1B=struct2cell(electrode1B);electrode1B=electrode1B{1,1};
    %         electrode2=load([path,'electrode2.mat']);electrode2=struct2cell(electrode2);electrode2=electrode2{1,1};
    %         electrode3=load([path,'electrode3.mat']);electrode3=struct2cell(electrode3);electrode3=electrode3{1,1};
    %         electrode4=load([path,'electrode4.mat']);electrode4=struct2cell(electrode4);electrode4=electrode4{1,1};
    %         imclose_judge=load([path,'imclose_judge.mat']);imclose_judge=struct2cell(imclose_judge);imclose_judge=imclose_judge{1,1};
    %         imclose_r=load([path,'imclose_r.mat']);imclose_r=struct2cell(imclose_r);imclose_r=imclose_r{1,1};
    %         VIclose_nrec=load([path,'VIclose_nrec.mat']);VIclose_nrec=struct2cell(VIclose_nrec);VIclose_nrec=VIclose_nrec{1,1};
    %         electrode5=load([path,'electrode5.mat']);electrode5=struct2cell(electrode5);electrode5=electrode5{1,1};
    %         electrode5_S=load([path,'electrode5_S.mat']);electrode5_S=struct2cell(electrode5_S);electrode5_S=electrode5_S{1,1};

            global electrode3_n electrode3_sub electrode8 electrode9 e9InnerConC electrode10 electrode11 T0B T0C electrode2B ...
                electrode12 electrode13 electrode14 electrode15 electrode16 electrode17
            [electrode3_n,electrode3_sub,electrode8,electrode9,e9InnerConC,electrode10,electrode11,T0B,T0C,electrode2B,...
                electrode12,electrode13,electrode14,electrode15,electrode16,electrode17]=ContactSeg2(ElecNum,T0,Te,regCT,brainmask,...
                brainmask2,brainmask3,brainmask_edge2,electrode1,electrode3,electrode5,electrode5_S,f);

    %         electrode3_n=load([path,'electrode3_n.mat']);electrode3_n=struct2cell(electrode3_n);electrode3_n=electrode3_n{1,1};
    %         electrode3_sub=load([path,'electrode3_sub.mat']);electrode3_sub=struct2cell(electrode3_sub);electrode3_sub=electrode3_sub{1,1};
    %         electrode8=load([path,'electrode8.mat']);electrode8=struct2cell(electrode8);electrode8=electrode8{1,1};
    %         electrode9=load([path,'electrode9.mat']);electrode9=struct2cell(electrode9);electrode9=electrode9{1,1};
    %         e9InnerConC=load([path,'e9InnerConC.mat']);e9InnerConC=struct2cell(e9InnerConC);e9InnerConC=e9InnerConC{1,1};
    %         electrode10=load([path,'electrode10.mat']);electrode10=struct2cell(electrode10);electrode10=electrode10{1,1};
    %         electrode11=load([path,'electrode11.mat']);electrode11=struct2cell(electrode11);electrode11=electrode11{1,1};
    %         electrode13=load([path,'electrode13.mat']);electrode13=struct2cell(electrode13);electrode13=electrode13{1,1};
    %         T0B=load([path,'T0B.mat']);T0B=struct2cell(T0B);T0B=T0B{1,1};
    %         T0C=load([path,'T0C.mat']);T0C=struct2cell(T0C);T0C=T0C{1,1};
    %         electrode2B=load([path,'electrode2B.mat']);electrode2B=struct2cell(electrode2B);electrode2B=electrode2B{1,1};
    %         electrode12=load([path,'electrode12.mat']);electrode12=struct2cell(electrode12);electrode12=electrode12{1,1};
    %         electrode13=load([path,'electrode13.mat']);electrode13=struct2cell(electrode13);electrode13=electrode13{1,1};
    %         electrode14=load([path,'electrode14.mat']);electrode14=struct2cell(electrode14);electrode14=electrode14{1,1};
    %         electrode15=load([path,'electrode15.mat']);electrode15=struct2cell(electrode15);electrode15=electrode15{1,1};
    %         electrode16=load([path,'electrode16.mat']);electrode16=struct2cell(electrode16);electrode16=electrode16{1,1};
    %         electrode17=load([path,'electrode17.mat']);electrode17=struct2cell(electrode17);electrode17=electrode17{1,1};

            global M_T1_recon RAS_T1_recon
            global e17Cendpoints e17Cendpoints_judge e17Centerps e17Ctargetps e17Cepsj_ras e17ElesLength ...
                e17EiContacts_cell e17Contacts      
            [e17Cendpoints,e17Cendpoints_judge,e17Centerps,e17Ctargetps,e17Cepsj_ras,e17ElesLength,e17EiContacts_RASC,...
                e17EiContacts_sub,e17EiContacts_cell,e17Contacts]=ContactSeg3(T0C,T1_IOP,T1_IPP1,regCT,brainmask2,brainmask3,...
                brainmask_edge2,electrode12,electrode17,M_T1_recon,RAS_T1_recon,f);
            SaveDataStr=get(handles.SaveDataType,'String');
            L=size(SaveDataStr,1);
            SaveDataStr2=cell(L+2,1);
            for i=1:L
                SaveDataStr2{i,1}=SaveDataStr{i,1};
            end
            SaveDataStr2{L+1,1}='e17EiContacts_RASC';
            SaveDataStr2{L+2,1}='e17EiContacts_sub';
            set(handles.SaveDataType,'String',SaveDataStr2);

            global ROI
            ROI=load('ROI_MNI_V4_List.mat');
            ROI=struct2cell(ROI);
            ROI=ROI{1,1};

            global AAL
            AAL=load('AAL.mat');
            AAL=struct2cell(AAL);
            AAL=AAL{1,1};
            
            global BN_Atlas
            BN_Atlas=load('BN_Atlas.mat');
            BN_Atlas=struct2cell(BN_Atlas);
            BN_Atlas=BN_Atlas{1,1};
            
    %         path='D:\电脑资料2019.03.17开始\毕设\电极定位\subjects\wu chuyin\step2\';
    %         e17EiContacts_sub=load([path,'e17EiContacts_sub.mat']);e17EiContacts_sub=struct2cell(e17EiContacts_sub);e17EiContacts_sub=e17EiContacts_sub{1,1};
    %         e17EiContacts_RASC=load([path,'e17EiContacts_RASC.mat']);e17EiContacts_RASC=struct2cell(e17EiContacts_RASC);e17EiContacts_RASC=e17EiContacts_RASC{1,1};
    
            global MR2avg_T
            [MR2avg_T,ElecInfo_cell]=ContactSeg4(M_avg152T1,e17EiContacts_RASC,T12avg_geomtform,ROI,AAL,BN_Atlas,f);
            SaveDataStr=get(handles.SaveDataType,'String');
            L=size(SaveDataStr,1);
            SaveDataStr2=cell(L+1,1);
            for i=1:L
                SaveDataStr2{i,1}=SaveDataStr{i,1};
            end
            SaveDataStr2{L+1,1}='ElecInfo_cell';
            set(handles.SaveDataType,'String',SaveDataStr2);

        end

        waitbar(0.95,f,'Complete contact segmentation...');

        global ElecInfo_cell2
        ElecInfo_cell2=ElecInfo_cell;
        for i=1:size(ElecInfo_cell,1)
            Contact=ElecInfo_cell{i,6};
            [MNIr,MNIc,MNIs]=ind2sub([91 109 91],Contact);
            NowConMNI=['[',num2str(MNIr),', ',num2str(MNIc),', ',num2str(MNIs),']'];
            ElecInfo_cell2{i,6}=NowConMNI;
        end

        global ChElStr AllContactsNum
        ChElStr=num2cell(0:ElecNum);
        ChElStr{1}='Please choose electrode：';
        set(handles.ChooseElectrode,'Value',1);
        set(handles.ChooseElectrode,'String',ChElStr);
        set(handles.ChooseContact,'Value',1);
        set(handles.RASCoordinate,'String',[]);
        set(handles.MNIIndex,'String',[]);
        set(handles.AALBrainArea,'String',[]);
        set(handles.BNBrainArea,'String',[]);
        axes(handles.TDPlane);hold off

        AllContactsNum=zeros(1,ElecNum);
        for i=1:ElecNum
            AllContactsNum(i)=size(e17EiContacts_sub{i},1);
        end

        global Images Images_recon T1 T1_recon graymin2 graymax2 SagCtr2 CorCtr2 AxiCtr2 nx2 ny2 nz2
        Images=T1;
        Images_recon=T1_recon;
        [nx2,ny2,nz2]=size(T1_recon);
        graymin2=min(min(min(T1_recon)));
        graymax2=max(max(max(T1_recon)));
        SagCtr2=T1_recon(:,:,round(nz2/2));
        CorCtr2=fliplr(squeeze(T1_recon(:,round(ny2/2),:)));
        AxiCtr2=fliplr(squeeze(T1_recon(round(nx2/2),:,:)));
        global IrAV IrSV IrCV
        IrAV=AxiCtr2;
        IrSV=SagCtr2;
        IrCV=CorCtr2;

        global ratio IrS_PixelSpacing IrS_SliceThickness nx3 ny3 nz3
        IrS_PixelSpacing=[0.6;0.6].*ratio;
        IrS_SliceThickness=0.6*ratio;
        coordinate_judge=1;

        global IrSampling M_IrSampling RAS_IrSampling RAS_IrSampling_WorldLimits
        IrSampling=T1_recon(1:ratio:nx2,1:ratio:ny2,1:ratio:nz2);
        [nx3,ny3,nz3]=size(IrSampling);
        [M_IrSampling,RAS_IrSampling]=transform2RAS(IrSampling,IrS_PixelSpacing,IrS_SliceThickness,T1_IOP,T1_IPP1,coordinate_judge);
        RAS_IrSampling_rmin=min(min(min(RAS_IrSampling(1,:,:,:))));
        RAS_IrSampling_rmax=max(max(max(RAS_IrSampling(1,:,:,:))));
        RAS_IrSampling_amin=min(min(min(RAS_IrSampling(2,:,:,:))));
        RAS_IrSampling_amax=max(max(max(RAS_IrSampling(2,:,:,:))));
        RAS_IrSampling_smin=min(min(min(RAS_IrSampling(3,:,:,:))));
        RAS_IrSampling_smax=max(max(max(RAS_IrSampling(3,:,:,:))));
        RAS_IrSampling_WorldLimits=[RAS_IrSampling_rmin RAS_IrSampling_rmax RAS_IrSampling_amin RAS_IrSampling_amax RAS_IrSampling_smin RAS_IrSampling_smax];

        axes(handles.SagittalPlane);
        imshow(SagCtr2,[graymin2 graymax2]);
        axes(handles.CoronalPlane);
        imshow(CorCtr2,[graymin2 graymax2]);
        axes(handles.AxialPlane);
        imshow(AxiCtr2,[graymin2 graymax2]);

        set(handles.AxialSlider,'Max',nx2);
        set(handles.SagittalSlider,'Max',nz2);
        set(handles.CoronalSlider,'Max',ny2);

        set(handles.AxialSlider,'Min',1);
        set(handles.SagittalSlider,'Min',1);
        set(handles.CoronalSlider,'Min',1);

        set(handles.AxialSlider,'Value',round(nx2/2));
        set(handles.SagittalSlider,'Value',round(nz2/2));
        set(handles.CoronalSlider,'Value',round(ny2/2));

        set(handles.WindowMin,'Max',graymax2);
        set(handles.WindowMin,'Value',graymax2);
        set(handles.WindowMin,'Min',graymin2);
        set(handles.WindowMin,'Value',graymin2);

        set(handles.WindowMax,'Max',graymax2);
        set(handles.WindowMax,'Value',graymax2);
        set(handles.WindowMax,'Min',graymin2);

        waitbar(1,f,'Complete contact segmentation...');

        global Axial_WinMin Axial_WinMax AxialColormap AxialColormap_Label AxialAlpha
        global Sagittal_WinMin Sagittal_WinMax SagittalColormap SagittalColormap_Label SagittalAlpha
        global Coronal_WinMin Coronal_WinMax CoronalColormap CoronalColormap_Label CoronalAlpha
        global TDColormap TDColormap_Label TDAlpha

        Axial_WinMin=graymin2;Axial_WinMax=graymax2;AxialColormap=gray;AxialColormap_Label=1;AxialAlpha=1;
        Sagittal_WinMin=graymin2;Sagittal_WinMax=graymax2;SagittalColormap=gray;SagittalColormap_Label=1;SagittalAlpha=1;
        Coronal_WinMin=graymin2;Coronal_WinMax=graymax2;CoronalColormap=gray;CoronalColormap_Label=1;CoronalAlpha=1;
        TDColormap=gray;TDColormap_Label=1;TDAlpha=1;

        global Colormap AlphaValue Plane
        Plane='Aixal Plane';Colormap=gray;AlphaValue=1;
        set(handles.ChoosePlane,'Value',1);
        set(handles.ChooseColormap,'Value',1);
        set(handles.Trans,'String',[]);

        global xslice2 yslice2 zslice2
        xslice2=round(nx3/2);
        yslice2=round(ny3/2);
        zslice2=round(nz3/2);
        axes(handles.TDPlane);
        SliceInRASC3(IrSampling,RAS_IrSampling,1,xslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,RAS_IrSampling,2,yslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        SliceInRASC3(IrSampling,RAS_IrSampling,3,zslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
        set(handles.TDEdit,'String','3D view of contacts in RAS coordinate system');
        set(handles.TDPlane,'XColor','w');
        set(handles.TDPlane,'YColor','w');
        set(handles.TDPlane,'ZColor','w');

        close(f);

    end
    
end

% --- Executes on button press in OutportContactsStats.
function OutportContactsStats_Callback(hObject, eventdata, handles)
% hObject    handle to OutportContactsStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ElecInfo_cell ElecInfo_cell2
[filename,pathname]=uiputfile('*.xlsx','Export the information of all contacts to an Excel format file');
if filename~=0
    f=waitbar(0,'The information of all contacts is being exported to an Excel format file, please be patient...');
    if isempty(ElecInfo_cell2)
        ElecInfo_cell2=ElecInfo_cell;
        for i=1:size(ElecInfo_cell,1)
            Contact=ElecInfo_cell{i,6};
            [MNIr,MNIc,MNIs]=ind2sub([91 109 91],Contact);
            NowConMNI=['[',num2str(MNIr),', ',num2str(MNIc),', ',num2str(MNIs),']'];
            ElecInfo_cell2{i,6}=NowConMNI;
        end
    end
    
    waitbar(0.2,f);
    filename=filename(1:end-5);
    filename1=[pathname,filename,'_ContactInformation.xlsx'];
    all_sheet={filename};
    sheet1name=char(all_sheet(1));
    
    waitbar(0.3,f);
    RenameSheet(filename1,sheet1name);
    
    waitbar(0.6,f);
    note1={'Electrode_ID','Contact_Num',...
        'R','A','S',...
        'MNI_Index',...
        'AAL_ID','AAL_Nom_C','AAL_Nom_L',...
        'BNA_ID','BNA_Lobe','BNA_Gyrus','BNA_Hemisphere','BNA_Descriptions'};
    sheet=sheet1name;
    xlRange1='A1';
    xlswrite(filename1,note1,sheet,xlRange1);
    waitbar(0.8,f);

    xlRange2='A2';
    note2=ElecInfo_cell2;
    xlswrite(filename1,note2,sheet,xlRange2);
    waitbar(1,f);
    close(f);
end

% --- Executes on selection change in ChooseElectrode.
function ChooseElectrode_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseElectrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChooseElectrode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChooseElectrode

global NonReconJudge ReconstructionJudge RASTransformJudge RegistrationJudge NonChooseContact SegmentationJudge
NonReconJudge=0;ReconstructionJudge=0;RASTransformJudge=1;RegistrationJudge=0;NonChooseContact=1;SegmentationJudge=0;

global e17EiContacts_sub e17EiContacts_RASC
global ElecLabel NowElec_sub NowElec_RASC ContactsNum ChConStr
ElecLabel=get(hObject,'Value')-1;

if ElecLabel>0
    
    set(handles.RASCoordinate,'String',[]);
    set(handles.MNIIndex,'String',[]);
    set(handles.AALBrainArea,'String',[]);
    set(handles.BNBrainArea,'String',[]);
    
    NowElec_sub=e17EiContacts_sub{ElecLabel};
    NowElec_RASC=e17EiContacts_RASC{ElecLabel};
    ContactsNum=size(NowElec_sub,1);
    ChConStr=num2cell(0:ContactsNum);
    ChConStr{1}='Please choose contact：';
    set(handles.ChooseContact,'Value',1);
    set(handles.ChooseContact,'String',ChConStr);

    global IrAV IrSV IrCV
    global Axial_WinMin Axial_WinMax AxialColormap AxialAlpha
    global Sagittal_WinMin Sagittal_WinMax SagittalColormap SagittalAlpha
    global Coronal_WinMin Coronal_WinMax CoronalColormap CoronalAlpha
    global TDColormap TDAlpha
    
    axes(handles.SagittalPlane);
    imshow(IrSV,[Sagittal_WinMin Sagittal_WinMax]);
    colormap(handles.SagittalPlane,SagittalColormap);
    alpha(handles.SagittalPlane,SagittalAlpha);
    axes(handles.CoronalPlane);
    imshow(IrCV,[Coronal_WinMin Coronal_WinMax]);
    colormap(handles.CoronalPlane,CoronalColormap);
    alpha(handles.CoronalPlane,CoronalAlpha);
    axes(handles.AxialPlane);
    imshow(IrAV,[Axial_WinMin Axial_WinMax]);
    colormap(handles.AxialPlane,AxialColormap);
    alpha(handles.AxialPlane,AxialAlpha);

    global IrSampling RAS_IrSampling RAS_IrSampling_WorldLimits ElecNum xslice2 yslice2 zslice2
    CurrentView=get(handles.TDPlane,'View');
    axes(handles.TDPlane);
    SliceInRASC3(IrSampling,RAS_IrSampling,1,xslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
    SliceInRASC3(IrSampling,RAS_IrSampling,2,yslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
    SliceInRASC3(IrSampling,RAS_IrSampling,3,zslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
    set(handles.TDPlane,'XColor','w');
    set(handles.TDPlane,'YColor','w');
    set(handles.TDPlane,'ZColor','w');
    colormap(handles.TDPlane,TDColormap);
    alpha(handles.TDPlane,TDAlpha);
    axis equal;view(CurrentView);

    axes(handles.TDPlane);
    for k=1:ElecNum
        if k~=ElecLabel
            EiCrasc=e17EiContacts_RASC{k};
            Scatter3InRASC5(EiCrasc,RAS_IrSampling_WorldLimits,9,'o',0.5,'b','b',0.75,0.75);hold on
        end
    end
    axes(handles.TDPlane);
    Scatter3InRASC5(NowElec_RASC,RAS_IrSampling_WorldLimits,9,'o',0.5,'g','g',0.75,0.75);hold on
    hold off;

end


% --- Executes during object creation, after setting all properties.
function ChooseElectrode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChooseElectrode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'Max',30);


% --- Executes on selection change in ChooseContact.
function ChooseContact_Callback(hObject, eventdata, handles)
% hObject    handle to ChooseContact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ChooseContact contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ChooseContact

global NonReconJudge ReconstructionJudge RASTransformJudge RegistrationJudge NonChooseContact SegmentationJudge
NonReconJudge=0;ReconstructionJudge=0;RASTransformJudge=0;RegistrationJudge=0;NonChooseContact=0;SegmentationJudge=1;

global NowElec_sub NowElec_RASC ElecLabel AllContactsNum ElecInfo_cell
global ContactLabel NowCon_sub NowCon_RASC InfoLabel NowCon_Info

ContactLabel=get(hObject,'Value')-1;

if (ElecLabel>0)&&(ContactLabel>0)
    
    NowCon_sub=NowElec_sub(ContactLabel,:);
    NowCon_RASC=NowElec_RASC(ContactLabel,:);
    InfoLabel=sum(AllContactsNum(1:(ElecLabel-1)))+ContactLabel;

    NowCon_Info=cell(1,14);
    for i=1:14
        NowCon_Info{1,i}=ElecInfo_cell{InfoLabel,i};
    end

    global NowConRAS NowConMNI 
    global NowConAAL NowCon_AALID NowCon_AALNomL
    global NowCon_BNID NowCon_BNLobe NowCon_BNGyrus NowCon_BNHemisphere NowCon_BNDescriptions NowConBN

    NCRASvpa=roundn(NowCon_RASC,-2);
    NowConRAS=['RAS：[',num2str(NCRASvpa(1)),', ',num2str(NCRASvpa(2)),', ',num2str(NCRASvpa(3)),']'];
    MNIind=NowCon_Info{1,6};
    [MNIr,MNIc,MNIs]=ind2sub([91 109 91],MNIind);
    NowConMNI=['MNI：[',num2str(MNIr),', ',num2str(MNIc),', ',num2str(MNIs),']'];

    NowCon_AALID=['AAL_ID：',num2str(NowCon_Info{1,7})];
    NowCon_AALNomL=['AAL_Nom_L：',NowCon_Info{1,9}];
    NowConAAL={NowCon_AALID,NowCon_AALNomL};

    NowCon_BNID=['BNA_ID：',num2str(NowCon_Info{1,10})];
    NowCon_BNLobe=['BNA_Lobe：',NowCon_Info{1,11}];
    NowCon_BNGyrus=['BNA_Lobe：',NowCon_Info{1,12}];
    NowCon_BNHemisphere=['BNA_Lobe：',NowCon_Info{1,13}];
    NowCon_BNDescriptions=['BNA_Lobe：',NowCon_Info{1,14}];
    NowConBN={NowCon_BNID,NowCon_BNLobe,NowCon_BNGyrus,NowCon_BNHemisphere,NowCon_BNDescriptions};

    set(handles.RASCoordinate,'String',NowConRAS);
    set(handles.MNIIndex,'String',NowConMNI);
    set(handles.AALBrainArea,'String',NowConAAL);
    set(handles.BNBrainArea,'String',NowConBN);

    global T1_recon IrAV IrSV IrCV nx2 ny2 nz2
    global Axial_WinMin Axial_WinMax AxialColormap AxialAlpha
    global Sagittal_WinMin Sagittal_WinMax SagittalColormap SagittalAlpha
    global Coronal_WinMin Coronal_WinMax CoronalColormap CoronalAlpha
    global TDColormap TDAlpha
    
    [~,~,nz2]=size(T1_recon);
    IrAV=fliplr(squeeze(T1_recon(NowCon_sub(1),:,:)));
    IrSV=T1_recon(:,:,NowCon_sub(3));
    IrCV=fliplr(squeeze(T1_recon(:,NowCon_sub(2),:)));
        
    axes(handles.SagittalPlane);
    imshow(IrSV,[Sagittal_WinMin Sagittal_WinMax]);hold on;
    colormap(handles.SagittalPlane,SagittalColormap);
    alpha(handles.SagittalPlane,SagittalAlpha);
    scatter(NowCon_sub(2),NowCon_sub(1),16,'filled','o','Linewidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r',...
        'MarkerEdgeAlpha',0.75,'MarkerFaceAlpha',0.75);hold on;
    hold off
    
    axes(handles.CoronalPlane);
    imshow(IrCV,[Coronal_WinMin Coronal_WinMax]);hold on;
    colormap(handles.CoronalPlane,CoronalColormap);
    alpha(handles.CoronalPlane,CoronalAlpha);
    scatter(nz2-NowCon_sub(3)+1,NowCon_sub(1),16,'filled','o','Linewidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r',...
        'MarkerEdgeAlpha',0.75,'MarkerFaceAlpha',0.75);hold on;
    hold off
    
    axes(handles.AxialPlane);
    imshow(IrAV,[Axial_WinMin Axial_WinMax]);hold on;
    colormap(handles.AxialPlane,AxialColormap);
    alpha(handles.AxialPlane,AxialAlpha);
    scatter(nz2-NowCon_sub(3)+1,NowCon_sub(2),16,'filled','o','Linewidth',0.5,'MarkerEdgeColor','r','MarkerFaceColor','r',...
        'MarkerEdgeAlpha',0.75,'MarkerFaceAlpha',0.75);hold on;
    hold off
    
    set(handles.AxialSlider,'Value',nx2-NowCon_sub(1)+1);
    set(handles.SagittalSlider,'Value',NowCon_sub(3));
    set(handles.CoronalSlider,'Value',ny2-NowCon_sub(2)+1);

    global IrSampling RAS_IrSampling RAS_IrSampling_WorldLimits ElecNum e17EiContacts_RASC ratio xslice2 yslice2 zslice2
    CurrentView=get(handles.TDPlane,'View');    
    xslice2=round(NowCon_sub(1)/ratio);
    yslice2=round(NowCon_sub(2)/ratio);
    zslice2=round(NowCon_sub(3)/ratio);
    axes(handles.TDPlane);
    SliceInRASC3(IrSampling,RAS_IrSampling,1,xslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
    SliceInRASC3(IrSampling,RAS_IrSampling,2,yslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
    SliceInRASC3(IrSampling,RAS_IrSampling,3,zslice2,RAS_IrSampling_WorldLimits,handles.TDPlane);hold on;
    set(handles.TDPlane,'XColor','w');
    set(handles.TDPlane,'YColor','w');
    set(handles.TDPlane,'ZColor','w');
    colormap(handles.TDPlane,TDColormap);
    alpha(handles.TDPlane,TDAlpha);
    axis equal;view(CurrentView);

    axes(handles.TDPlane);
    for k=1:ElecNum
        if k~=ElecLabel
            EiCrasc=e17EiContacts_RASC{k};
            Scatter3InRASC5(EiCrasc,RAS_IrSampling_WorldLimits,9,'o',0.5,'b','b',0.75,0.75);hold on
        end
    end
    
    axes(handles.TDPlane);
    NowElec_RASC_comp=setdiff(NowElec_RASC,NowCon_RASC,'Rows');
    Scatter3InRASC5(NowElec_RASC_comp,RAS_IrSampling_WorldLimits,9,'o',0.5,'g','g',0.75,0.75);hold on
    
    axes(handles.TDPlane);
    Scatter3InRASC5(NowCon_RASC,RAS_IrSampling_WorldLimits,25,'o',0.5,'r','r',0.75,0.75);hold on
    axis equal;view(CurrentView);hold off;
    
end

% --- Executes during object creation, after setting all properties.
function ChooseContact_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ChooseContact (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'Max',30);


function RASCoordinate_Callback(hObject, eventdata, handles)
% hObject    handle to RASCoordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of RASCoordinate as text
%        str2double(get(hObject,'String')) returns contents of RASCoordinate as a double


% --- Executes during object creation, after setting all properties.
function RASCoordinate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RASCoordinate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function MNIIndex_Callback(hObject, eventdata, handles)
% hObject    handle to MNIIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MNIIndex as text
%        str2double(get(hObject,'String')) returns contents of MNIIndex as a double


% --- Executes during object creation, after setting all properties.
function MNIIndex_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MNIIndex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function AALBrainArea_Callback(hObject, eventdata, handles)
% hObject    handle to AALBrainArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of AALBrainArea as text
%        str2double(get(hObject,'String')) returns contents of AALBrainArea as a double


% --- Executes during object creation, after setting all properties.
function AALBrainArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to AALBrainArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'Max',10);


function BNBrainArea_Callback(hObject, eventdata, handles)
% hObject    handle to BNBrainArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of BNBrainArea as text
%        str2double(get(hObject,'String')) returns contents of BNBrainArea as a double


% --- Executes during object creation, after setting all properties.
function BNBrainArea_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BNBrainArea (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'Max',20);
