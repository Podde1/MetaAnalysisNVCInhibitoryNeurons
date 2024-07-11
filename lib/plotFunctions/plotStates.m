function [f] = plotStates(sol, varargin)
% sol = AMICI solution struct
% p = parmter vector

f=figure();
numberOfPlots = size(sol.x,2)+size(sol.y,2); 

nRows = floor(sqrt(numberOfPlots));
nCols =  ceil(sqrt(numberOfPlots));

if nRows*nCols < numberOfPlots
    nCols = nCols + 1;
end

%% construct the titles 
if nargin == 1 
    TitleText = [strcat('x',cellfun(@num2str,num2cell(1:size(sol.x,2)),'Uni',0)), strcat('y',cellfun(@num2str,num2cell(1:size(sol.y,2)),'Uni',0)),strcat('var',cellfun(@num2str,num2cell(1:numberOfPlots-(size(sol.x,2)+size(sol.y,2))),'Uni',0))];
else 
    TitleText  = varargin{1};
end

%%
for i=1:size(sol.x,2) % loop over each state x 
subplot(nRows,nCols,i);
plot(sol.t,sol.x(:,i),'LineWidth',2)
title(TitleText{i});    
end

k=1;
for j = i+1:size(sol.y,2)+i % loop over each obserable y
    subplot(nRows,nCols,j);
    plot(sol.t,sol.y(:,k),'LineWidth',2)
    title(TitleText{j});
    k=k+1;
end


end
