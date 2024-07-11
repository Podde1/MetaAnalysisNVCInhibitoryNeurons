function [Figure] = maximizeAxes(Figure,nRows,nColumns,axOrder)
% [Figure] = untitled(Figure,n,m) - adjusts the outerpositions of the axes
% in the figure such that the figure margin is minimised. 
% Figure is a matlab figure handel
% nRows and nColumns are the number of rows and columnd of subplots presens. 
% axOrder is an index matrix which is used if the order of the subplots
% differ from the matlab default order. e.i. left to right, top to bottom.

if ~isa(Figure,'matlab.ui.Figure')
error('%s is NOT a valid figure handle!',inputname(1));
return
else

    if length(Figure.Children)<nRows*nColumns% if the number of subplots are less then nRows*nColumns i.e. if the lastrow are missing some plots
        figure(Figure.Number)% set the current figure to be the active figure
        for i = [length(Figure.Children)+1:nRows*nColumns]% for the subplots thtat are missing to fill out the last row
            subplot(nRows,nColumns,i)% add each subplot
        end
    end
   axesIndex = false(length(Figure.Children),1);
   for k = 1:length(Figure.Children)
       axesIndex(k) = isa(Figure.Children(k),'matlab.graphics.axis.Axes');
   end
   axes = flipud(Figure.Children(axesIndex));
   
   if  length(axes) > nRows*nColumns
       error('inputs n=%d and m=%d does not match the number of axes in figure %s',nRows,nColumns,inputname(1));
       error('Make sure that %s corresopnds to the number of rows of subplots and %s corresopnds to the number of columns of subplots',inputname(2),inputname(3))
       return
   else
       if nargin < 4 || isempty(axOrder)
           axOrder = [1:length(axes)];
       elseif nargin == 4 && ~isempty(axOrder)
           if length(unique(axOrder)) ~= length(axes)
               error('There are %d axes in the figure but only %d unique indecies in the vector axOrder! Ensure that the two match in length and that all elements in axOrder are uique.',length(axes),length(unique(axOrder)));
               return
           end
       end
       if length(axes) < nRows*nColumns
           for i = (nRows*nColumns)-length(axes):-1:1
              subplot(nRows,nColumns,(nColumns*nRows+1)-1)
          end
       end
       
       AxPositions=cell(nRows*nColumns,1);
       for i = 1:nRows
           for j = 0:(nColumns-1)
               AxPositions(i+j+(i-1)*(nColumns-1))={[(j*1/nColumns),((nRows-i)/nRows),(1/nColumns),(1/nRows)]};
           end
       end
       set(axes,{'OuterPosition'},AxPositions(axOrder));
   end
   
   %Aling the x.positions and the widths of the columns
   for i = 1:nColumns
       for j = 0:nRows-1
           axes(i+j*nColumns).Position([1,3])= axes(i).Position([1,3]);
       end
   end
   
   %Aling the y.positions and the hights of the rows
   for i = 1:nRows
       for j = (1:nColumns)+(i-1)*nColumns
           axes(j).Position([2,4])= axes(1+nColumns*(i-1)).Position([2,4]);
       end
   end
   
end

end

