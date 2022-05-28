% ohlc.m - make Open-High-Low-Close bar chart
%
% Usage: h = ohlc(t,Y,varargin)
%        h = ohlc(t,Y,'s',scale,varargin)
%
% t = time axis, column or row of length(t) = size(Y,1)
% Y = matrix with 4 columns [Open,High,Low,Close], or,
%            with 3 columns [High,Low,Close], or,
%            with 2 columns [High,Low]
% varargin = specifies linestyle and/or color     
% scale    = scale factor of the forward/backward horizontal close/open prices, 
%            default=1/3, must be in the range 0 <= scale <= 0.5
%
% h = column vector of line handles, useful for legends
%     h(1) is for Y, h(2:end) for the other lines
%
% examples:  ohlc(t,[O,H,L,C]);  
%            ohlc(t,[H,L,C]);
%            ohlc(t,[H,L]);
%            ohlc(t,[O,H,L,C], 's',1/2); 
%            ohlc(t,[O,H,L,C], 'color','r');
%            ohlc(t,[O,H,L,C], 'linestyle','--');
%            ohlc(t,[O,H,L,C], 's',1/5, 'linestyle','--');
%           
% example: the file MSFT.dat has column structure: [O, H, L, C, Vol, AdjClose, serial_date, string_date] 
%    A = loadfile('MSFT.dat');                  % load file from OSP toolbox   
%    i = find(A(:,7) == datenum('01/02/09'));   % column 7 holds the serial dates, start at Jan.02,2009
%    Y = A(i:i+59, 1:4);                        % extract 60 rows, keep first 4 columns [O,H,L,C]
%    t = (0:size(Y,1)-1);                       % for real date use t = A(i:i+29,7);                   
%    figure; ohlc(t,Y);         xlim([-1,60]); yaxis(14,22, 14:2:22);  % Y = [O,H,L,C]
%    figure; ohlc(t,Y,'s',1/2); xlim([-1,60]); yaxis(14,22, 14:2:22);  % scale = 1/2
%    figure; ohlc(t,Y(:,2:4));  xlim([-1,60]); yaxis(14,22, 14:2:22);  % Y = [H,L,C]
%    figure; ohlc(t,Y(:,2:3));  xlim([-1,60]); yaxis(14,22, 14:2:22);  % Y = [H,L]
%    [L,U] = pbands(Y, 7);                                                     % 7-day projection bands
%    figure; h = ohlc(t,Y, t,L,'r--', t,U,'g--');  
%    xaxis(-1,60, 0:10:60); yaxis(14,22, 14:2:22);              % OHLC bars with bands
%    legend(h, ' data', ' lower', ' upper', 'location','ne');   % must use h in legend to get the right colors
%
% see also OHLCYY, YYLIM

% S. J. Orfanidis - 2009-2018

function h = ohlc(t,Y,varargin)

if nargin==0, help ohlc; return; end

scale = 1/3;                                % default horizontal bar scale parameter

if nargin>3 & varargin{1}=='s',                        % read another scale value
   scale = varargin{2};  
   varargin{1}=[]; varargin{2} = [];        
end

dt = (max(t)-min(t))/(length(t)-1);         % determine bar spacing - assumes equally-spaced data
tf = t + scale*dt;                          % forward bar for closing prices
tb = t - scale*dt;                          % backward bar for open prices

[N,K] = size(Y); K = min(4,K);              % restrict K to 4,3,2,1           

switch K,
   case 4
      O = Y(:,1); H = Y(:,2); L = Y(:,3); C = Y(:,4); 
   case 3, 
      H = Y(:,1); L = Y(:,2); C = Y(:,3); 
   case 2,
      H = Y(:,1); L = Y(:,2); 
   otherwise,
      disp(' '); disp('Y must be [Open,High,Low,Close,...], or [High,Low,Close], or [High,Low]'); return;
end

T = NaN(3,N); T(1,:) = t; T(2,:) = t;  T = T(:);              % high/low vertical bars, T = [t1,t1,nan, t2,t2,nan, ...]'
Z = NaN(3,N); Z(1,:) = L; Z(2,:) = H;  Z = Z(:);              %                         Z = [L1,H1,nan, L2,H2,nan, ...]'

if K>=3,
   Tf = NaN(3,N); Tf(1,:) = t; Tf(2,:) = tf;  Tf = Tf(:);     % closing horizontal bars, Tf = [t1,tf1,nan, t2,tf2,nan, ...]'
   F  = NaN(3,N);  F(1,:) = C;  F(2,:) = C;    F = F(:);      %                           F = [C1, C1,nan, C2, C2,nan, ...]'
end

if K==4,
   Tb = NaN(3,N); Tb(1,:) = t; Tb(2,:) = tb;   Tb = Tb(:);    % opening horizontal bars, Tb = [t1,tb1,nan, t2,tb2,nan, ...]'
   B  = NaN(3,N);  B(1,:) = O;  B(2,:) = O;     B = B(:);     %                           B = [O1, O1,nan, O2, O2,nan, ...]'
end

switch K,
   case 4,
      h = plot(T,Z,'b', Tf,F,'b-', Tb,B,'b-', varargin{:});       % include opening and closing bars
      h(1:2)=[];                                                  % helps with legend
   case 3, 
      h = plot(T,Z,'b', Tf,F,'b-', varargin{:});                  % include closing bars only
      h(1)=[];                                                    % helps with legend
   case 2,
      h = plot(T,Z,'b', varargin{:});                             % high/low vertical bars only
end

xlim([min(t)-dt, max(t)+dt]);             % makes end-point bars more visible


