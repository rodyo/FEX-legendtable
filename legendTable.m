% TODO
%{
WindowButtonMotionFcn to couple together table+legend when moving them with the mouse

make header optional

location of legend determines direction for resizing/repositioning

%}


% 1) legend handle is given
% function LegendTable(legendHandle, tableHead, tableContent)
% function LegendTable(legendHandle, table)

% 2) find legend handle automatically%
% function LegendTable(tableHead, tableContent)
% function LegendTable(table)


% ans =
%     'North'
%     'South'
%     'East'
%     'West'
%     'NorthEast'
%     'SouthEast'
%     'NorthWest'
%     'SouthWest'
%     'NorthOutside'
%     'SouthOutside'
%     'EastOutside'
%     'WestOutside'
%     'NorthEastOutside'
%     'SouthEastOutside'
%     'NorthWestOutside'
%     'SouthWestOutside'
%     'Best'                % UNSUPPORTED
%     'BestOutside'         % UNSUPPORTED
%     'none'                % UNSUPPORTED


function legendTable(varargin)
        
    %% Initialize
    
    argc = nargin;
    
    % Legend handle
    if ishandle(varargin{1})
        legendHandle = varargin{1};
        if ~strcmpi(get(legendHandle, 'tag'), 'legend')
            warning('LegendTable:invalid_legend_handle',...
                'Invalid legend handle given. Not drawing table...');
            return;
        end
        argOffset = 1;
    else
        figChildren = get(gcf, 'children');
        legendHandle = figChildren(strcmpi(get(figChildren, 'tag'), 'legend'));
        if isempty(legendHandle)
            warning('LegendTable:no_legend_handle',...
                'No legend found. Not drawing table...');
            return;
        end
        argOffset = 0;
    end
    
    % Table contents
    if argc == 2+argOffset
        tableHead = varargin{1+argOffset};
        tableContent = varargin{2+argOffset};
    elseif argc == 1+argOffset
        table = varargin{1+argOffset};
        tableHead    = table(1,:);
        tableContent = table(2:end,:);
    else
        warning('LegendTable:no_table_data',...
            'No table data received. Not drawing table...');
        return;
    end
    
    %% Prepare data
    
    % Convert header & content to cell-array when necessary
    if isnumeric(tableContent)
        tableContent = cellfun(@num2str, ...
            num2cell(tableContent), 'UniformOutput', false);
    end
    if isnumeric(tableHead)
        tableHead = cellfun(@num2str, ...
            num2cell(tableHead), 'UniformOutput', false);
    end
    
    % Get() all relevant info
    legendPosition = get(legendHandle, 'position');
    plotHandle     = get(legendHandle, 'parent');
    legendlocation = get(legendHandle, 'Location');
    children = get(legendHandle, 'children');
    labels   = children(strcmp(get(children, 'type'), 'text'));
    
    
     % Basic error traps
    if size(tableContent,1) ~= numel(labels)
        error('LegendTable:dimension_mismatch',...
            'Each legend entry must have a corresponding row in the table.')
    end
    
    if size(tableHead,2) ~= size(tableContent,2)
        error('LegendTable:dimension_mismatch',...
            'Table header dimensions are inconsistent with table data.');
    end
    
    
    
    % Resize/relocate direction depends on legend location string    
    switch lower(regexprep(legendlocation, 'inside|outside$', '', 'ignorecase'))
        case 'north'
            legendMove = [-1/2 -1 0 0]
            tableMove  = [+1/2 -1 0 0];
            
        case 'south'
            legendMove = [-1/2 0 0 0]
            tableMove  = [+1/2 0 0 0];
            
        case 'east'
            legendMove = [-1 -1 0 0]
            tableMove  = [ 0 -1 0 0];
            
        case 'west'
            legendMove = [ 0 -1 0 0]
            tableMove  = [+1 -1 0 0];
            
        case 'northeast'
            % HERE
            legendMove = [ 0 -1 0 0]
            tableMove  = [+1 -1 0 0];
            
        case 'northwest'
            legendMove = [ 0 -1 0 0]
            tableMove  = [+1 -1 0 0];
                        
        case 'southeast'
            legendMove = [ 0 -1 0 0]
            tableMove  = [+1 -1 0 0];
            
        case 'southwest'
            legendMove = [ 0 -1 0 0]
            tableMove  = [+1 -1 0 0];
                        
        otherwise
            warning('LegendTable:unsupported_legend_location',...
                'Unsupported legend location: ''%s'' Using default directions...', legendlocation)
    end
        
    % initial location of the axes gridlines ("table cells")
    xticks = linspace(0, 1, numel(tableHead)+1);
    yticks = linspace(0, 1, numel(labels)+2);
    
    % Locations of the text are in the centers of the cells
    txt_xPositions = xticks(1:end-1) + (xticks(2)-xticks(1))/2;
    txt_yPositions = fliplr(yticks(1:end-1) + (yticks(2)-yticks(1))/2);
    
    % Set initial position/size of table/legend
    headerHeight  = legendPosition(4)/numel(labels);
    tablePosition = legendPosition + [0 -headerHeight 0 headerHeight];
    
    % Shift position of legend
    % TODO: generalize by looking at legend positions; this is only valid for
    % 'NorthEast'
    set(legendHandle, 'position', legendPosition + [-tablePosition(3) -headerHeight 0 0])
    
    
    %% Draw initial table
    
    % Create table
    table = axes(...
        'position', tablePosition,...
        'xtick', xticks,...
        'ytick', yticks,...
        'xticklabel', [],...
        'yticklabel', [],...
        'gridlinestyle', '-',...
        'box', 'on',...
        'tag', 'LegendTable');
    grid on
    
    
    % Print table header & table entries
    kk = 1;
    tableTexts = zeros(numel(tableHead)+numel(tableContent),1);
    for ii = 1:numel(txt_xPositions)
        
        % Column header
        tableTexts(kk) = text(txt_xPositions(ii), txt_yPositions(1), tableHead{ii},...
            'parent', table,...
            'HorizontalAlignment', 'center');
        kk = kk + 1;
        
        % Column content
        for jj = 1:numel(txt_yPositions)-1
            tableTexts(kk) = text(txt_xPositions(ii), txt_yPositions(jj+1), tableContent{jj,ii},...
                'parent', table,...
                'HorizontalAlignment', 'center',...
                'VerticalAlignment'  , 'middle');
            kk = kk + 1;
        end
    end
    
    
    %% Corrections
    
    % Adjust table & legend size based on table content
    maxExtents = get(tableTexts, 'extent');
    maxExtents = max(cat(1,maxExtents{:}));
    
    % X-direction
    if maxExtents(3) > 1/size(tableContent,2)
        
    end
    
    % Y-direction
    if maxExtents(4) > 1/(size(tableContent,1)+1) %header
        
    end
    
    
    
end
