function varargout = TNC_SS_GUI(varargin)
% TNC_SS_GUI MATLAB code for TNC_SS_GUI.fig
%      TNC_SS_GUI, by itself, creates a new TNC_SS_GUI or raises the existing
%      singleton*.
%
%      H = TNC_SS_GUI returns the handle to a new TNC_SS_GUI or the handle to
%      the existing singleton*.
%
%      TNC_SS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TNC_SS_GUI.M with the given input arguments.
%
%      TNC_SS_GUI('Property','Value',...) creates a new TNC_SS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before TNC_SS_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to TNC_SS_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help TNC_SS_GUI

% Last Modified by GUIDE v2.5 15-Nov-2013 12:01:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @TNC_SS_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @TNC_SS_GUI_OutputFcn, ...
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

% --- Executes just before TNC_SS_GUI is made visible.
function TNC_SS_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
    % This function has no output args, see OutputFcn.
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    % varargin   command line arguments to TNC_SS_GUI (see VARARGIN)

    % Choose default command line output for TNC_SS_GUI
    handles.output          = hObject;

    handles.numGraphs       = 3;
    handles.allFeatColor    = [0.65 0.65 0.65];
    handles.segList         = [3 2 4];
    handles.featureData     = [];
    handles.markStyle       = '.';
    handles.markSize        = 6;
    handles.path            = '';
    handles.shankNum        = 1;
    
    handles.avgWaves        = 1;
    handles.showXCorr       = 0;
    handles.updateWF        = 0;
    
    handles.xPlotName       = 'PC1'; 
    handles.xPlotNum        = 1; 
    handles.yPlotName       = 'PC2'; 
    handles.yPlotNum        = 2; 
    handles.zPlotName       = 'none'; 
    handles.zPlotNum        = 0; 
    
    handles.allDims         = zeros(1,5);    
    handles.lblTrue         = ones(1,5);
    handles.cntTrue         = zeros(1,5); 
    handles.conTrue         = zeros(1,5);
    handles.boundMethod     = 'ellipse';
    
    handles.numShanks       = 8;
    handles.numSegs         = 8;
    
    handles.addId           = 1;
    handles.delId           = 1;
    handles.cropId          = 1;
    handles.pickId          = 1;
    handles.reNumStart      = 0;
    handles.reNumEnd        = 0;
    handles.clustToProp     = 0;
    handles.confBound       = 2.5; 
    handles.scaler          = 0.25;
    
    handles.currentAxes     = 1;
    handles.cMap            = [ 92,200,200;
                               134,223,  4;
                               204,  3, 95;
                                14, 83,167;
                               247,254,  0;
                               220,  0, 85;
                               255,116,  0;
                               113,  9,170;
                                 0,153,153;
                               230, 60,137;
                                 2,146,146;
                               186,239,110;                               
                               249,175,114 ];
                           
    handles.cMap            = handles.cMap ./ 256;
    handles.numClusts       = ones(1,5).*7;
    
    % Update handles structure
    guidata(hObject, handles);

    % This sets up the initial plot - only do when we are invisible
    % so window can get raised using TNC_SS_GUI.
    if strcmp(get(hObject,'Visible'),'off')

        TNC_SS_UpdateGraphs(handles);

    end
    
% UIWAIT makes TNC_SS_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = TNC_SS_GUI_OutputFcn(hObject, eventdata, handles)
    % varargout  cell array for returning output args (see VARARGOUT);
    % hObject    handle to figure
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Get default command line output from handles structure
    varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
    % hObject    handle to pushbutton1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    axes(handles.axes1);
    cla;


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
    % hObject    handle to OpenMenuItem (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    file = uigetfile('*.fig');
    if ~isequal(file, 0)
        open(file);
    end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
    % hObject    handle to PrintMenuItem (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
    % hObject    handle to CloseMenuItem (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                         ['Close ' get(handles.figure1,'Name') '...'],...
                         'Yes','No','Yes');
    if strcmp(selection,'No')
        return;
    end

    delete(handles.figure1)


% --- Executes on selection change in selectData.
function selectData_Callback(hObject, eventdata, handles)
    % hObject    handle to selectData (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = get(hObject,'String') returns selectData contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from selectData

    contents = get(hObject,'String');
    handles.fileName = contents{get(hObject,'Value')};

    % Update handles structure
    guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function selectData_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to selectData (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
         set(hObject,'BackgroundColor','white');
    end

    clear listOfFiles
    d = dir('*_ft*');
    
    if numel(d)>0
        for i=1:numel(d)
            listOfFiles{i} = d(i).name;
        end
        handles.fileName = d(1).name;
        set(hObject, 'String', listOfFiles);
    else
        handles.fileName = 'none';
        set(hObject, 'String', 'none');
    end
    
    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in loadData.
function loadData_Callback(hObject, eventdata, handles)
    % hObject    handle to loadData (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    disp('Loading the feature data...');
    S = load(handles.fileName);
    forSort = TNC_SS_CreateSortStruct(S.featStruct);
    handles.featureData = forSort;

    disp('Loading the waveform data...');
    nameLength = numel(handles.fileName);
    spikeFile = [handles.fileName(1:nameLength-6) 'ss.mat'];
    S2 = load(spikeFile);
    handles.sessionStruct = S2.sessionStruct;
    
    % Update the menu to reflect the plotting options
    set(handles.xAxDisp, 'String', forSort.paramNames);
    set(handles.yAxDisp, 'String', forSort.paramNames);
    forSort.paramNames{1} = 'none';
    set(handles.zAxDisp, 'String', forSort.paramNames);

    handles.numShanks = numel(forSort.seg(1).shank);
    for i = 1:handles.numShanks
        shankNums{i} = num2str(i);
    end
    set(handles.selectShank, 'String', shankNums);
    
    handles.numSegs = numel(forSort.seg);
    for i = 1:handles.numSegs
        segList{i} = num2str(i);
    end
    set(handles.selectSlice, 'String', segList);

    if handles.numSegs<3

        set(handles.selectSlice, 'Value', 1);
            
        handles.segList = [0 -1 1] + 1;

        for i=1:numel(handles.segList)
            if handles.segList(i) < 1
                handles.segList(i) = 1;
            elseif handles.segList(i) > handles.numSegs
                handles.segList(i) = handles.numSegs;
            end            
        end
        
    end
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);

    % --- Executes on selection change in cMark.
function cMark_Callback(hObject, eventdata, handles)
    % hObject    handle to cMark (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns cMark contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from cMark
    contents = cellstr(get(hObject,'String'));
    chosenColor = contents{get(hObject,'Value')};

    switch lower(chosenColor)

        case 'white'
            handles.allFeatColor = [1 1 1];

        case 'grey'
            handles.allFeatColor = [0.65 0.65 0.65];

        case 'light grey'
            handles.allFeatColor = [0.75 0.75 0.75];

        case 'dark grey'
            handles.allFeatColor = [0.33 0.33 0.33];

        otherwise
            handles.allFeatColor = [0.5 0.5 0.5];

    end

    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);

    
% --- Executes during object creation, after setting all properties.
function cMark_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to cMark (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in dnMark.
function dnMark_Callback(hObject, eventdata, handles)
    % hObject    handle to dnMark (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    if handles.markSize > 1
        handles.markSize = handles.markSize-1;
    else
        handles.markSize = 1;    
    end

    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);

% --- Executes on button press in upMark.
function upMark_Callback(hObject, eventdata, handles)
    % hObject    handle to upMark (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    handles.markSize = handles.markSize+1;

    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);

% --- Executes on selection change in sMark.
function sMark_Callback(hObject, eventdata, handles)
    % hObject    handle to sMark (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns sMark contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from sMark

    contents = cellstr(get(hObject,'String'));
    chosenColor = contents{get(hObject,'Value')};

    switch lower(chosenColor)

        case 'dot'
            handles.markStyle = '.';

        case 'circle'
            handles.markStyle = 'o';

        case 'square'
            handles.markStyle = 's';

        case 'plus'
            handles.markStyle = '+';

        otherwise
            handles.markStyle = '.';

    end

    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);


% --- Executes during object creation, after setting all properties.
function sMark_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to sMark (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on selection change in xAxDisp.
function xAxDisp_Callback(hObject, eventdata, handles)
    % hObject    handle to xAxDisp (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns xAxDisp contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from xAxDisp

    contents = cellstr(get(hObject,'String'));
    handles.xPlotName = contents{get(hObject,'Value')};    
    handles.xPlotNum = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);

% --- Executes during object creation, after setting all properties.
function xAxDisp_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to xAxDisp (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on button press in dfMark.
function dfMark_Callback(hObject, eventdata, handles)
    % hObject    handle to dfMark (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    handles.markSize = 6;

    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);

% --- Executes on selection change in yAxDisp.
function yAxDisp_Callback(hObject, eventdata, handles)
    % hObject    handle to yAxDisp (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns yAxDisp contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from yAxDisp


    contents = cellstr(get(hObject,'String'));
    handles.yPlotName = contents{get(hObject,'Value')};
    handles.yPlotNum = get(hObject,'Value');

    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);    
    
% --- Executes during object creation, after setting all properties.
function yAxDisp_CreateFcn(hObject, eventdata, handles)
    % hObject    handle to yAxDisp (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    empty - handles not created until after all CreateFcns called

    % Hint: popupmenu controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end


% --- Executes on selection change in zAxDisp.
function zAxDisp_Callback(hObject, eventdata, handles)
    % hObject    handle to zAxDisp (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns zAxDisp contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from zAxDisp

    contents = cellstr(get(hObject,'String'));
    handles.zPlotName = contents{get(hObject,'Value')};
    handles.zPlotNum = get(hObject,'Value');
    
    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);  

% --- Executes during object creation, after setting all properties.
function zAxDisp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zAxDisp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in manClustAdd1.
function manClustAdd1_Callback(hObject, eventdata, handles)
    % hObject    handle to manClustAdd1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)
    
    % Let the user add a cluster
    newId = handles.addId;
    [newClustIds] = TNC_SS_AddCluster(handles,newId);
    handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id = newClustIds;
    
    % Update handles structure (again)
    guidata(hObject, handles);

    % Calculate centers for all projections
    [handles] = TNC_SS_UpdateClusterCenters(handles);

    % Update handles structure (again)
    guidata(hObject, handles);
    
    % Update all graphs based upon settings in the handles structure
    drawnow; pause(1);
    TNC_SS_UpdateGraphs(handles);   

% --- Executes on button press in manClustDelete1.
function manClustDelete1_Callback(hObject, eventdata, handles)
% hObject    handle to manClustDelete1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
        
    disp(['Deleting cluster ' num2str(handles.delId) '...']);

    if handles.delId==0
        [newClustIds] = TNC_SS_DelCluster(handles);
        handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id = newClustIds;
    else
        toDel = find(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id==handles.delId)
        handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id(toDel) = 0;
    end
    
    % Update handles structure (again)
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);   

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in selectSlice.
function selectSlice_Callback(hObject, eventdata, handles)
    % hObject    handle to selectSlice (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns selectSlice contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from selectSlice

    contents = cellstr(get(hObject,'String'));
    
    handles.segList = [0 -1 1] + get(hObject,'Value');
    
    for i=1:numel(handles.segList)
        if handles.segList(i) < 1
            handles.segList(i) = 1;
        elseif handles.segList(i) > handles.numSegs
            handles.segList(i) = handles.numSegs;
        end            
    end
    
    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);  

    
% --- Executes during object creation, after setting all properties.
function selectSlice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in selectShank.
function selectShank_Callback(hObject, eventdata, handles)
    % hObject    handle to selectShank (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: contents = cellstr(get(hObject,'String')) returns selectShank contents as cell array
    %        contents{get(hObject,'Value')} returns selected item from selectShank

    handles.shankNum = get(hObject,'Value');
        
    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);  

% --- Executes during object creation, after setting all properties.
function selectShank_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectShank (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in writeData.
function writeData_Callback(hObject, eventdata, handles)
% hObject    handle to writeData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Update all graphs based upon settings in the handles structure
    TNC_SS_WriteData(handles);
    TNC_SS_SaveSortedIds(handles);


% --- Executes on selection change in propClustMenu1.
function propClustMenu1_Callback(hObject, eventdata, handles)
% hObject    handle to propClustMenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns propClustMenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from propClustMenu1

 
    contents = cellstr(get(hObject,'String'));
    
    clustToProp = contents{get(hObject,'Value')};
    if strcmp(lower(clustToProp),'all')
        disp('Propagate all clusters.');
    else
        handles.clustToProp = str2double(clustToProp);
        disp(['Propagate cluster ' clustToProp '.']);
    end

    % Update handles structure (again)
    guidata(hObject, handles);
    
    
% --- Executes during object creation, after setting all properties.
function propClustMenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to propClustMenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in propClustRight1.
function propClustRight1_Callback(hObject, eventdata, handles)
% hObject    handle to propClustRight1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)\

    % Recompute the cluster centers and eccentricities
    [handles] = TNC_SS_UpdateClusterCenters(handles);
    handles.clusToProp = 0;

    % Update handles structure (again)
    guidata(hObject, handles);

    % Propagate the cluster centers for currently selected cluster number
    [handles] = TNC_SS_PropClust(handles.segList(1),handles.segList(1)+1,handles);

    % Update handles structure (again)
    guidata(hObject, handles);
    
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);  
    

% --- Executes on button press in propClustLeft1.
function propClustLeft1_Callback(hObject, eventdata, handles)
% hObject    handle to propClustLeft1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % Recompute the cluster centers and eccentricities
    [handles] = TNC_SS_UpdateClusterCenters(handles);
    handles.clusToProp = 0;

    % Update handles structure (again)
    guidata(hObject, handles);

    % Propagate the cluster centers for currently selected cluster number
    [handles] = TNC_SS_PropClust(handles.segList(1),handles.segList(1)-1,handles);

    % Update handles structure (again)
    guidata(hObject, handles);
    
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles); 
    
% --- Executes on button press in manClustMerge1.
function manClustMerge1_Callback(hObject, eventdata, handles)
% hObject    handle to manClustMerge1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function cropClustEdit_Callback(hObject, eventdata, handles)
% hObject    handle to cropClustEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cropClustEdit as text
%        str2double(get(hObject,'String')) returns contents of cropClustEdit as a double

    % update variable
    handles.cropId = str2double(get(hObject,'String'));
    
    % Update handles structure
    guidata(hObject, handles);
    
    
% --- Executes during object creation, after setting all properties.
function cropClustEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cropClustEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in labelChk.
function labelChk_Callback(hObject, eventdata, handles)
% hObject    handle to labelChk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of labelChk


% --- Executes during object creation, after setting all properties.
function axes5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes5


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes on button press in setPath.
function setPath_Callback(hObject, eventdata, handles)
    % hObject    handle to setPath (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    dirname = uigetdir;
    handles.path = dirname;
    
    cd(dirname); 
    
    clear listOfFiles
    d = dir('*_ft*');
    
    if numel(d)>0
        for i=1:numel(d)
            listOfFiles{i} = d(i).name;
        end
        handles.fileName = d(1).name;
        set(handles.selectData, 'String', listOfFiles);
    else
        handles.fileName = 'none';
        set(handles.selectData, 'String', 'none');
    end
    

    % Update handles structure
    guidata(hObject, handles);

% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in boundChk.
function boundChk_Callback(hObject, eventdata, handles)
% hObject    handle to boundChk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of boundChk


% --- Executes on button press in dispLabels1.
function dispLabels1_Callback(hObject, eventdata, handles)
% hObject    handle to dispLabels1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dispLabels1

    handles.lblTrue(1) = get(hObject,'Value');

    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);  
    
% --- Executes on button press in dispCenters1.
function dispCenters1_Callback(hObject, eventdata, handles)
% hObject    handle to dispCenters1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dispCenters1

    handles.cntTrue(1) = get(hObject,'Value');

    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);  
    
% --- Executes on button press in dispBounds1.
function dispBounds1_Callback(hObject, eventdata, handles)
% hObject    handle to dispBounds1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dispBounds1

    handles.conTrue(1) = get(hObject,'Value');

    % Update handles structure
    guidata(hObject, handles);

    [handles] = TNC_SS_UpdateClusterCenters(handles);

    % Update handles structure
    guidata(hObject, handles);
    
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);  
    
% --- Executes on selection change in ColorPaletteMenu1.
function ColorPaletteMenu1_Callback(hObject, eventdata, handles)
% hObject    handle to ColorPaletteMenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ColorPaletteMenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ColorPaletteMenu1


% --- Executes during object creation, after setting all properties.
function ColorPaletteMenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ColorPaletteMenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in kMeansButton1.
function kMeansButton1_Callback(hObject, eventdata, handles)
    % hObject    handle to kMeansButton1 (see GCBO).cnt
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    disp(['Running kmeans clustering with ' num2str(handles.numClusts(1)) ' clusters...']);
    
    if handles.allDims==0
        [unit] = TNC_EventCluster(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params(:,[handles.xPlotNum,handles.yPlotNum]),handles.numClusts(1),'km',0);
        handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id  = unit.clustIds;
        handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).cnt = zeros(handles.numClusts(1),size(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params,2));
        handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).cnt(:,handles.xPlotNum) = unit.clustCenters(:,1);
        handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).cnt(:,handles.yPlotNum) = unit.clustCenters(:,2);
        handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).cnt
    else
        [unit] = TNC_EventCluster(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).params,handles.numClusts(1),'km',0);
        handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id  = unit.clustIds;
        handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).cnt = unit.clustCenters;  
        handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).cnt
    end

    % Update handles structure
    guidata(hObject, handles);
    
    % Calculate centers for all projections
    [handles] = TNC_SS_UpdateClusterCenters(handles);

    % Update handles structure (again)
    guidata(hObject, handles);
    
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);  
    
    
function kMeansNumber1_Callback(hObject, eventdata, handles)
    % hObject    handle to kMeansNumber1 (see GCBO)
    % eventdata  reserved - to be defined in a future version of MATLAB
    % handles    structure with handles and user data (see GUIDATA)

    % Hints: get(hObject,'String') returns contents of kMeansNumber1 as text
    %        str2double(get(hObject,'String')) returns contents of kMeansNumber1 as a double

    handles.numClusts(1) = str2double(get(hObject,'String'));
    
    % Update handles structure
    guidata(hObject, handles);
    
% --- Executes during object creation, after setting all properties.
function kMeansNumber1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kMeansNumber1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in JumpMethodButton.
function JumpMethodButton_Callback(hObject, eventdata, handles)
% hObject    handle to JumpMethodButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function JumpMethodTxt_Callback(hObject, eventdata, handles)
% hObject    handle to JumpMethodTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of JumpMethodTxt as text
%        str2double(get(hObject,'String')) returns contents of JumpMethodTxt as a double


% --- Executes during object creation, after setting all properties.
function JumpMethodTxt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to JumpMethodTxt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ManSeedButton.
function ManSeedButton_Callback(hObject, eventdata, handles)
% hObject    handle to ManSeedButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function ClustNumDisp1_Callback(hObject, eventdata, handles)
% hObject    handle to ClustNumDisp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ClustNumDisp1 as text
%        str2double(get(hObject,'String')) returns contents of ClustNumDisp1 as a double


% --- Executes during object creation, after setting all properties.
function ClustNumDisp1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClustNumDisp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function kMeansButton1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to kMeansButton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in allDims1.
function allDims1_Callback(hObject, eventdata, handles)
% hObject    handle to allDims1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allDims1

    handles.allDims(1) = get(hObject,'Value'); 
    
    % Update handles structure
    guidata(hObject, handles);


% --- Executes on button press in clearClust1.
function clearClust1_Callback(hObject, eventdata, handles)
% hObject    handle to clearClust1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    resetSize = size(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id);
    handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id = zeros(resetSize);
    
    % Update handles structure
    guidata(hObject, handles);
    
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);  

% --- Executes on button press in useConfidence1.
function useConfidence1_Callback(hObject, eventdata, handles)
% hObject    handle to useConfidence1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    % For points that are not in, leave them alone unless they have the
    % number of the current cluster; then they should be zeroed
    [newIds] = TNC_SS_CropCluster(handles,0);

    % Update handles structure
    handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id = newIds;
    guidata(hObject, handles);
    
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);  

% --- Executes on button press in hullCheck1.
function hullCheck1_Callback(hObject, eventdata, handles)
% hObject    handle to hullCheck1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hullCheck1

    if get(hObject,'Value')==1
        handles.boundMethod = 'hull';
    else
        handles.boundMethod = 'ellipse';
    end
    
    % Update handles structure
    guidata(hObject, handles);        


% --- Executes on button press in reNum1.
function reNum1_Callback(hObject, eventdata, handles)
% hObject    handle to reNum1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if handles.reNumStart==0

        clustNums = unique(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id);

        for p=2:numel(clustNums) % first number will be 0 for unsorted clusters
            theseIds = find(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id==clustNums(p));
            handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id(theseIds) = p-1;
        end

    else

            theseIds = find(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id==handles.reNumStart);
            handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id(theseIds) = handles.reNumEnd;

    end

    % Update handles structure
    guidata(hObject, handles);
    
     
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);     
    

function delClustEdit_Callback(hObject, eventdata, handles)
% hObject    handle to delClustEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of delClustEdit as text
%        str2double(get(hObject,'String')) returns contents of delClustEdit as a double

    % update variable
    handles.delId = str2double(get(hObject,'String'));
    
    % Update handles structure
    guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function delClustEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to delClustEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

    

function addClustEdit_Callback(hObject, eventdata, handles)
% hObject    handle to addClustEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of addClustEdit as text
%        str2double(get(hObject,'String')) returns contents of addClustEdit as a double

    % update variable
    handles.addId = str2double(get(hObject,'String'));

    % Update handles structure
    guidata(hObject, handles);

    
% --- Executes during object creation, after setting all properties.
function addClustEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to addClustEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function reNumStart_Callback(hObject, eventdata, handles)
% hObject    handle to reNumStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reNumStart as text
%        str2double(get(hObject,'String')) returns contents of reNumStart as a double

    handles.reNumStart = str2double(get(hObject,'String'));

    % Update handles structure
    guidata(hObject, handles);
    


% --- Executes during object creation, after setting all properties.
function reNumStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reNumStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function reNumEnd_Callback(hObject, eventdata, handles)
% hObject    handle to reNumEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of reNumEnd as text
%        str2double(get(hObject,'String')) returns contents of reNumEnd as a double


    handles.reNumEnd = str2double(get(hObject,'String'));

    % Update handles structure
    guidata(hObject, handles);
    

% --- Executes during object creation, after setting all properties.
function reNumEnd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to reNumEnd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in allShankSave.
function allShankSave_Callback(hObject, eventdata, handles)
% hObject    handle to allShankSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of allShankSave


% --- Executes on button press in loadSorted.
function loadSorted_Callback(hObject, eventdata, handles)
% hObject    handle to loadSorted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    fileToLoad = [handles.fileName(1:numel(handles.fileName)-4) '_tns.mat'];    
    
    [handles] = TNC_SS_LoadSortedIds(handles,fileToLoad);
    
    % Update handles structure
    guidata(hObject, handles);
     
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);     
    


% --- Executes on button press in saveSorted.
function saveSorted_Callback(hObject, eventdata, handles)
% hObject    handle to saveSorted (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    TNC_SS_SaveSortedIds(handles);

function confSD_Callback(hObject, eventdata, handles)
% hObject    handle to confSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of confSD as text
%        str2double(get(hObject,'String')) returns contents of confSD as a double
    
    handles.confBound = str2double(get(hObject,'String'));

    disp(['Set the confidence bound to ' num2str(str2double(get(hObject,'String'))) ' x s.d.']);
    
    % Update handles structure
    guidata(hObject, handles);
        
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles); 
    

% --- Executes during object creation, after setting all properties.
function confSD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to confSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox9.
function checkbox9_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox9


% --- Executes on button press in autoPickButton.
function autoPickButton_Callback(hObject, eventdata, handles)
% hObject    handle to autoPickButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    % Take a wild guess at the cluster boundaries from the choice
    if handles.pickId == 0
        handles.pickId = max(unique(handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id))+1;
    end
    
    % Update handles structure
    guidata(hObject, handles);

    [newIds]        = TNC_SS_GrowClusterBounds(handles);
    handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id = newIds;
    
    % Update handles structure
    guidata(hObject, handles);
 
    [handles]   = TNC_SS_UpdateClusterCenters(handles);

    % Update handles structure
    guidata(hObject, handles);
    
    [newIds]        = TNC_SS_FindOptimalBoundary(handles);
    handles.featureData.seg(handles.segList(1)).shank(handles.shankNum).id = newIds;

    % Update handles structure
    guidata(hObject, handles);   
    
    % Update handles structure
    guidata(hObject, handles);   
    
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles); 
    
    
function pickClustEdit_Callback(hObject, eventdata, handles)
% hObject    handle to pickClustEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pickClustEdit as text
%        str2double(get(hObject,'String')) returns contents of pickClustEdit as a double

    handles.pickId = str2double(get(hObject,'String'));
     
    % Update handles structure
    guidata(hObject, handles); 
    
    
% --- Executes during object creation, after setting all properties.
function pickClustEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pickClustEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mixModelButton.
function mixModelButton_Callback(hObject, eventdata, handles)
% hObject    handle to mixModelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function mixModelEdit_Callback(hObject, eventdata, handles)
% hObject    handle to mixModelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mixModelEdit as text
%        str2double(get(hObject,'String')) returns contents of mixModelEdit as a double


% --- Executes during object creation, after setting all properties.
function mixModelEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mixModelEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox10.
function checkbox10_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox10


% --- Executes on button press in checkbox11.
function checkbox11_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox11



function setBoundScale_Callback(hObject, eventdata, handles)
% hObject    handle to setBoundScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of setBoundScale as text
%        str2double(get(hObject,'String')) returns contents of setBoundScale as a double

    handles.scaler = str2double(get(hObject,'String'));
     
    % Update handles structure
    guidata(hObject, handles); 
        
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles); 

% --- Executes during object creation, after setting all properties.
function setBoundScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setBoundScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkDispWaves.
function checkDispWaves_Callback(hObject, eventdata, handles)
% hObject    handle to checkDispWaves (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkDispWaves

    handles.avgWaves = get(hObject,'Value');

    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);   

    
% --- Executes on button press in undoButton.
function undoButton_Callback(hObject, eventdata, handles)
% hObject    handle to undoButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    fileToLoad = [handles.fileName(1:numel(handles.fileName)-4) '_tns.mat'];    
    
    [handles] = TNC_SS_LoadSortedIds(handles,fileToLoad);
    
    % Update handles structure
    guidata(hObject, handles);
     
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);   


% --- Executes on button press in shankDown.
function shankDown_Callback(hObject, eventdata, handles)
% hObject    handle to shankDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if handles.shankNum>1
        handles.shankNum = handles.shankNum - 1;
    end
    
    set(handles.selectShank,'Value',handles.shankNum);
    
    % Update handles structure
    guidata(hObject, handles);
     
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles); 
    
% --- Executes on button press in shankUp.
function shankUp_Callback(hObject, eventdata, handles)
% hObject    handle to shankUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    if handles.shankNum < handles.numShanks
        handles.shankNum = handles.shankNum + 1;
    end
    
    set(handles.selectShank,'Value',handles.shankNum);

    % Update handles structure
    guidata(hObject, handles);
     
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles); 
    


% --- Executes on button press in segDown.
function segDown_Callback(hObject, eventdata, handles)
% hObject    handle to segDown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.segList = (handles.segList(1) + [0 -1 1]) - 1;
    
    for i=1:numel(handles.segList)
        if handles.segList(i) < 1
            handles.segList(i) = 1;
        elseif handles.segList(i) > handles.numSegs
            handles.segList(i) = handles.numSegs;
        end            
    end
    
    set(handles.selectSlice,'Value',handles.segList(1));

    % Update handles structure
    guidata(hObject, handles);
     
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles); 
    


% --- Executes on button press in segUp.
function segUp_Callback(hObject, eventdata, handles)
% hObject    handle to segUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    handles.segList = (handles.segList(1) + [0 -1 1]) + 1;
    
    for i=1:numel(handles.segList)
        if handles.segList(i) < 1
            handles.segList(i) = 1;
        elseif handles.segList(i) > handles.numSegs
            handles.segList(i) = handles.numSegs;
        end            
    end
    
    set(handles.selectSlice,'Value',handles.segList(1));

    % Update handles structure
    guidata(hObject, handles);
     
    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles); 
    


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over segUp.
function segUp_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to segUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over selectData.
function selectData_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to selectData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on setPath and none of its controls.
function setPath_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to setPath (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkDispXcorr.
function checkDispXcorr_Callback(hObject, eventdata, handles)
% hObject    handle to checkDispXcorr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkDispXcorr

    handles.showXCorr = get(hObject,'Value');

    % Update handles structure
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);   


% --- Executes during object creation, after setting all properties.
function setPath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setPath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in updateWFbutton.
function updateWFbutton_Callback(hObject, eventdata, handles)
% hObject    handle to updateWFbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    % Update handles structure
    handles.updateWF = 1;
    guidata(hObject, handles);

    % Update all graphs based upon settings in the handles structure
    TNC_SS_UpdateGraphs(handles);      
    
    % Revert handles structure
    handles.updateWF = 0;
    guidata(hObject, handles);
