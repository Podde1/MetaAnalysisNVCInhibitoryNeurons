function [strHandle] = addSubfigureNotations(figHandle,figLabelFontsize,varargin)
%% Adds  sub figure notations A, B, C etc. to the subfigures of a matlab plot. 
% Subfigures must be axes-objects. Notations will be added in alphanumerical order starting at "A", unless otherwise specified. 
%
% Inputs - figHandle - the matlab handle of the parent figure object. 
%          figLabelFontsize - The desired font size of the sub figur notations.  
%          startingNotation (optional) - a single case-sensitive letter specifying starting notation. 

    if nargin ==3 
        startingNotation = varargin{1};
        startingNotationIdx = double(startingNotation);
    else
        startingNotationIdx = 65; % set the starting notation to "A" if not otherwies specified. "A" has ascii code 65. 
    end
    figure_structure_type = get(figHandle.Children,'type');


    % differ between TiledChartLayouts and regular subplot layouts since matlab treats these differently. 
    % TiledChatLayouts adds another graphical object layer between the parent figure and the axes objects. 

    if strcmp(figure_structure_type,'tiledlayout') % if a tiledLayout is used
        axes_objects = flipud(figHandle.Children.Children(strcmp(get(figHandle.Children.Children,'type'),'axes')));% get only the axes objects in the correct order
       
    elseif iscell(figure_structure_type) % if a subplot structure is used
        axes_objects = flipud(figHandle.Children(strcmp(get(figHandle.Children,'type'),'axes'))); % get only the axes objects in the corect order
    end

    sgTitles = cellfun(@char,num2cell(startingNotationIdx + [0:length(axes_objects)-1]),'uni',false); % get the labels from A to the number of sub figures
        
    for i = 1:length(axes_objects)% loop over the axes objects 
        %add annotation-object
        strHandel(i) = annotation(figHandle,'textbox',[.01 .1 .3 .3],'String',sgTitles{i},'FitBoxToText','on','FontSize',figLabelFontsize,'FontWeight','bold','LineStyle','none');
        
        %set position of annotation to match axes-object. 
        strHandel(i).Position(1:2)=[axes_objects(i).Position(1),axes_objects(i).Position(2)+axes_objects(i).Position(4)-strHandel(i).Position(4)];
    end
end
