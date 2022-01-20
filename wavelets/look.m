function y = look(im,sc)
%look - View an image.
%------------------------------------------------------------------------------
%SYNOPSIS       look(I)
%                 Shows an image I in grayscale in a new window. A scaling
%                 0 to 255 is assumed.
%
%               look(I,[min max])
%                 Adjust dynamic range of colortable to [min max]. If value
%                 is NaN the value will be taken from the image itself.
%
%COMMENTS       If images toolbox is available imzoom is activated.
%
%SEE ALSO       image. (imzoom in images toolbox)
%------------------------------------------------------------------------------
%Jonas Svanberg                        (C) 1994 Image Coding Group. LiU, SWEDEN

% 960514 Added imzoom functionality. /svan
% 970905 Now uses gca() to create axes.

FRAME=6;

width=size(im,2);
height=size(im,1);

global ICG_BUTTON_X ICG_BUTTON_Y ICG_BUTTON_HEIGHT
ICG_BUTTON_HEIGHT=25;
%ICG_BUTTON_Y=height+3;
ICG_BUTTON_Y=1;

ss = get(0,'screensize');

% Create new figure

newfig = figure;
axes=gca;
% axes=get(newfig,'currentaxes'); % worked in ver.4 but not any longer
set(axes,'units','pixels','position', ...
    [FRAME ICG_BUTTON_Y+ICG_BUTTON_HEIGHT+FRAME width height]);
set(newfig,'units','pixels','color',[0 0 0]);

% Create buttons

ICG_BUTTON_X=0;
icgmakebutton('Dismiss', 60, 'lookcallback', 0);
icgmakebutton('Nrml', 40, 'lookcallback', 5);
icgmakebutton('2x', 25, 'lookcallback', 6);
icgmakebutton('Fit', 25, 'lookcallback', 4);
icgmakebutton('Cur', 30, 'lookcallback', 1);
icgmakebutton('Gray', 44, 'lookcallback', 2);
icgmakebutton('Pseudo', 56, 'lookcallback', 3);



% Set dimensions

scr_width=max(ICG_BUTTON_X , width + 2*FRAME);
scr_height=height+ICG_BUTTON_HEIGHT + 2*FRAME;
pos=get(newfig,'position');
pos(2) = min(pos(2),ss(4)-scr_height-60);
set(newfig, ...
    'position', [pos(1:2) scr_width scr_height], ...
    'userdata', [width height scr_width scr_height]);


% Display image etc.

colormap(gray(256));
if nargin < 2
  sc=[0 255];
end

low = sc(1);
high = sc(2);
if isnan(low)
  low = min(min(im));
end

if isnan(high)
  high = max(max(im))+1;
end

% We make a one pixel frame around image. This together with
% set(axes,'visible','off') will guarantee a one-color frame
% around the image making it suitable for xv grab and auto
% crop.

%im2 = zeros(size(im,1)+2, size(im,2)+2);
%im2(2:size(im,1)+1,2:size(im,2)+1)=
image((im-low)*255.999999/(high-low));
set(axes,'visible','off'); 

if which('imzoom')
  imzoom on
  fprintf(2,'imzoom activated.\n')
end


