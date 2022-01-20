function lookcallback(icg_fig,icg_key)
% lookcallback - Script called from the look utility buttons.
%

FRAME=6;

%fprintf(1, 'lookcallback called.\n');
%icg_fig
%icg_key

if icg_key == 0		% Kill window
  delete(icg_fig)
end

if icg_key == 1		% Make figure current
  gcf=icg_fig;
end

if icg_key == 2		% Colormap gray
  figure(icg_fig);
  colormap('gray');
end

if icg_key == 3		% Colormap pseudo
  figure(icg_fig);
  colormap((1+sin((1:256)'*rand(1,3)*0.3))/2);
end

if icg_key == 4         % Fit
  pos = get(icg_fig,'position');
  ax = get(icg_fig,'currentaxes');
  axp = get(ax,'pos');
  siz = pos(3:4)-axp(1:2)-[FRAME FRAME]+[0 1];
  set(ax,'position',[axp(1:2) siz]);
end

if icg_key == 5         % Normal size
  if which('imzoom')
    imzoom(-1000)
  end
  pos = get(icg_fig,'position');
  udat = get(icg_fig,'userdata');
  ax = get(icg_fig,'currentaxes');
  axp = get(ax,'pos');
  set(icg_fig, 'position', [pos(1:2),udat(3:4)]);
  set(ax, 'position', [axp(1:2) udat(1:2)]);
end  

if icg_key == 6         % Double size
  pos = get(icg_fig,'position');
  udat = get(icg_fig,'userdata');
  ax = get(icg_fig,'currentaxes');
  axp = get(ax,'pos');
  siz = 2*axp(3:4);
  set(icg_fig, 'position', ...
      [pos(1:2), max(pos(3:4), axp(1:2) + [FRAME FRAME] + siz)]);
  set(ax,'position',[axp(1:2) siz]);
end
  
