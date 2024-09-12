%% ViolinPlot
% 
% Creates a violin plot to present data distributions.
% 
% Description:
%   This function generates a violin plot with overlaid swarm points to visualize the distribution of data across multiple groups. 
%   Each row of the input data is shown in an independent swarmchart.
% 
% Inputs:
%   data    - [N, L] matrix where N is the number of groups and L is the number of observations per group. Each row represents a different group.
%   labs    - 1xN cell array of strings containing labels for each group. Length must match the number of rows in `data`.
%   ylab    - (Optional) String specifying the label for the y-axis. Default is an empty string.
%   tit     - (Optional) String specifying the title of the plot. Default is an empty string.
%   save    - (Optional) String specifying the path and filename to save the plot. If empty, the plot will not be saved.
% 
% Outputs:
%   None. The function does not return any outputs. 
%   It creates and optionally saves a figure showing the violin plot with swarm points.
% 
% Alberto Liardi, 2024


function [] = ViolinPlot(data,labs,ylab,tit,save)

    assert(size(data,1)==length(labs), "Mismatch in the dimesions of data and labels!");
    if ~exist('ylab','var') || isempty(ylab), ylab=""; end
    if ~exist('tit','var') || isempty(tit), tit=""; end
    if ~exist('save','var') || isempty(save), save=""; end

    % prepare the labels
    [N,L] = size(data);
    labels = strings(1,N*L);
    for i=1:N, labels(i,1+(i-1)*L:i*L) = labs(i); end
    
    % prepare the colormap
    cmap = colormap("jet");
    map = cmap(1:floor(size(cmap,1)/N):end-1,:);
    map2 = cell(N*L,1);
    for l = 1:L, for n = 1:N, map2{(n-1)*L+l} = map(n,:); end, end
    map = cell2mat(map2);
    
    % prepare the data
    x = categorical(labels, labs);
    y = reshape(data, 1, N*L);
    
    % finally create the plot
    fig = figure(); swarmchart(x,y,'.', 'SizeData',200, 'CData', map);
    title(tit); 
    ax = gca;
    set(gca,'FontName','CMU serif','FontSize',15);
    yline(0, 'k--'); ylabel(ylab,'FontSize',15,'interpreter','latex');
    if save~="", exportgraphics(fig,save,'Resolution',300); end

end