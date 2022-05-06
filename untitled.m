function varargout = untitled(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @untitled_OpeningFcn, ...
                   'gui_OutputFcn',  @untitled_OutputFcn, ...
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

function untitled_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;

guidata(hObject, handles);

set(gcf,'outerposition',get(0,'screensize'));

if strcmp(get(hObject,'Visible'),'off')
    plot(rand(5));
end

function varargout = untitled_OutputFcn(hObject, eventdata, handles)

varargout{1} = handles.output;



function pushbutton1_Callback(hObject, eventdata, handles)

axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        global im           

        im=imread('D:\MATLABR2016a\bin\GUIsimulate\image\fig_1.jpg');   

        axes(handles.axes1);  

        imshow(im);    
        
        global in   

        in=imread('D:\MATLABR2016a\bin\GUIsimulate\image\fig_2.jpg');   

        axes(handles.axes2);  

        imshow(in);   
    case 2
        global ip           

        ip=imread('D:\MATLABR2016a\bin\GUIsimulate\image\fig_3.jpg');   

        axes(handles.axes1);  

        imshow(ip);    
        
        global ir   

        ir=imread('D:\MATLABR2016a\bin\GUIsimulate\image\fig_4.jpg');   

        axes(handles.axes2);  

        imshow(ir);  
    case 3
        global is           

        is=imread('D:\MATLABR2016a\bin\GUIsimulate\image\fig_5.jpg');   

        axes(handles.axes1);  

        imshow(is);    
        
        global it   

        it=imread('D:\MATLABR2016a\bin\GUIsimulate\image\fig_6.jpg');   

        axes(handles.axes2);  

        imshow(it);          
    case 4
        global is           

        is=imread('D:\MATLABR2016a\bin\GUIsimulate\image\fig_7.jpg');   

        axes(handles.axes1);  

        imshow(is);    
        
        global it   

        it=imread('D:\MATLABR2016a\bin\GUIsimulate\image\fig_8.jpg');   

        axes(handles.axes2);  

        imshow(it);          
    case 5
        global is           

        is=imread('D:\MATLABR2016a\bin\GUIsimulate\image\fig_9.jpg');   

        axes(handles.axes1);  

        imshow(is);    
        
        global it   

        it=imread('D:\MATLABR2016a\bin\GUIsimulate\image\fig_10.jpg');   

        axes(handles.axes2);  

        imshow(it);          
end

function FileMenu_Callback(hObject, eventdata, handles)

function OpenMenuItem_Callback(hObject, eventdata, handles)

file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

function PrintMenuItem_Callback(hObject, eventdata, handles)

printdlg(handles.figure1)

function CloseMenuItem_Callback(hObject, eventdata, handles)

selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

function popupmenu1_Callback(hObject, eventdata, handles)

function popupmenu1_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'2D UAV path planning', '3D UAV path planning', 'delivery center generation FCD', 'drone delivery of IACO and ACO','drone delivery of GA and PSO'});

if ~exist('windowMaximize.dll','file')
    error('windowMaximize.dll not found.');
end

if nargin == 0
    h = gcf;
end

if ischar(h)

    if strcmpi(h,'all')

        h = findobj('Type','figure');            
    else

        error('Argument must be the correct string.');
    end
else

    for n=1:length(h)

        if ~ishandle(h(n)) || (~strcmp(get(h(n),'Type'),'root') && ...
                               ~strcmp(get(h(n),'Type'),'figure'))

            error('Argument(s) must be a correct handle(s).');
        end
    end
end

if h ~= 0

    for n=length(h):-1:1

        windowname = ['maximize_',num2str(h(n))];
 
        numTitle = get(h(n),'NumberTitle');
        figName = get(h(n),'Name');

        set(h(n),'Name',windowname,'NumberTitle','off');

        drawnow;

        windowMaximize(windowname,get(h(n),'Resize'));

        set(h(n),'Name',figName,'NumberTitle',numTitle);
    end
else

    windowMaximize('MATLAB');
end