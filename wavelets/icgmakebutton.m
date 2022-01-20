function icgmakebutton(name, width, cb, key)
%
%
%

global ICG_BUTTON_X ICG_BUTTON_Y ICG_BUTTON_HEIGHT

pos = [ICG_BUTTON_X ICG_BUTTON_Y width ICG_BUTTON_HEIGHT];
button = uicontrol('style', 'push', ... 
	'units', 'pixels', 'position', pos, ...
	'callback', ...
	sprintf('%s(%f,%f)', cb,gcf,key));

set(button, 'string', name);

ICG_BUTTON_X = ICG_BUTTON_X + width;


