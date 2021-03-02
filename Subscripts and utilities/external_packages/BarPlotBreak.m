function b=BarPlotBreak(Y,y_break_start,y_break_end,break_type,scale,x_lim)

% BarBreakPlot(y,y_break_start,y_break_end,break_type,scale)
% Produces a plot who's y-axis skips to avoid unecessary blank space
% 
% INPUT
% y
% y_break_start
% y_break_end
% break_type
%    if break_type='RPatch' the plot will look torn
%       in the broken space
%    if break_type='Patch' the plot will have a more
%       regular, zig-zag tear
%    if break_plot='Line' the plot will merely have
%       some hash marks on the y-axis to denote the
%       break
% scale = between 0-1, the % of max Y value that needs to subtracted from
% the max value bars
% USAGE:
% figure;
% BarPlotBreak([10,40,1000], 50, 60, 'RPatch', 0.85);
%
% Original Version developed by:
% Michael Robbins
% robbins@bloomberg.net
% michael.robbins@bloomberg.net
%
% Modified by: 
% Chintan Patel
% chintan.patel@dbmi.columbia.edu

% data
if nargin<5; break_type='RPatch'; end
if nargin<4; y_break_end=40; end
if nargin<3; y_break_start=10; end
if nargin<2; Y=[1:10,40:50]; end
% if nargin<1; x=rand(1,21); end

% y_break_mid  = (y_break_end-y_break_start)./2+y_break_start;

Y2 = Y;
Y2(Y2>=y_break_end) = Y2(Y2>=y_break_end)-(Y2(Y2>=y_break_end)*scale);

%find the max and min and cut max to 1.5 times the min
% bar(Y2, 0.5);
bar(Y2,'FaceColor',[0.8 0.8 0.8]);
hold on
axis tight
axis manual

xlim = get(gca,'xlim');
ytick = get(gca,'YTick');
[~,i] = min(ytick<=y_break_start);
% if all elements are the same, pick the last index
if all(ytick<=y_break_start)
    i = length(ytick);
end
y = (ytick(i)-ytick(i-1))./2+ytick(i-1);
dy = (ytick(2)-ytick(1))./5;
% xtick = get(gca,'XTick');
% % x = xtick(1);
x = 1;
% % dx = (xtick(2)-xtick(1))./2;
dx = 1;
switch break_type
    case 'Patch'
		% this can be vectorized
        dx = (xlim(2)-xlim(1))./10;
%         dx = (x_lim-xlim(1))./x_lim;
        yy = repmat([y-2.*dy y-dy],1,6);
        xx = xlim(1)+dx.*[0:11]-1;
		patch([xx(:);flipud(xx(:))], ...
            [yy(:);flipud(yy(:)-2.*dy)], ...
            [1 1 1],'EdgeColor','w')
        
        plot(xx,yy,'k')
        plot(xx,yy-2.*dy,'k')
%     case 'RPatch'
% 		% this can be vectorized
%         dx = (xlim(2)-xlim(1))./100;
%         yy = y+rand(101,1).*2.*dy;
%         xx = xlim(1)+dx.*(0:100);
% 		patch([xx(:);flipud(xx(:))], ...
%             [yy(:);flipud(yy(:)-2.*dy)], ...
%             [.8 .8 .8])
    case 'Line'
		line([x-dx x+dx   ],[y-2.*dy y-dy   ],'Color','k');
% 		line([x    x+dx],[y+dy    y+2.*dy]);
		line([x-dx x+dx   ],[y-3.*dy y-2.*dy],'Color','k');
% 		line([x    x+dx],[y+2.*dy y+3.*dy]);
    case 'Across'
        dx = (xlim(2)-xlim(1))./10;
%         dx = (x_lim-xlim(1))./x_lim;
        yy = repmat([y-2.*dy y-dy],1,6);
        xx = xlim(1)+dx.*[0:11];
        plot(xx,yy,'k')
        plot(xx,yy-2.*dy,'k')
end

%ytick(ytick>y_break_start)=ytick(ytick>y_break_start)+y_break_mid;

ytick(ytick>y_break_start) = ytick(ytick>y_break_start)+(Y(Y>=y_break_end)*scale);

for i=1:length(ytick)
   yticklabel{i}=sprintf('%d',ytick(i));
end
set(gca,'yticklabel',yticklabel);