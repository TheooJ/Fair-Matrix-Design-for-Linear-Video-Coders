function im = pgmload(name)
%pgmload - Read an PGM image from file.
%------------------------------------------------------------------------------
%SYNOPSIS       I = pgmload('STR')
%                 Opens a PGM image file named STR and returns a matrix
%                 consisting of values 0:255. It first searches current
%                 and then makes use of PICPATH global variable.
%
%SEE ALSO       pgmsave, the pbmplus package.
%------------------------------------------------------------------------------
%Jonas Svanberg                        (C) 1994 Image Coding Group. LiU, SWEDEN

global PICPATH

if nargin < 1
  error('Not enough input arguments.')
end

if ~isstr(name)
  error('Filename argument must be a string.')
end

fullname=name;

fid = fopen(fullname,'r');

if fid == -1 & exist('PICPATH')
  picpath=PICPATH;
  while prod(size(picpath))>0 & fid==-1
    [p,picpath] = strtok(picpath);
    fullname=[p '/' name];
    fid = fopen([p '/' name], 'r');
  end
end

if fid == -1
  error(['Can''t find ' name '.'])
end


% read and check header

typ = fscanf(fid, '%s', 1);

if ~ (size(typ,2) == 2)
  error(['File ' fullname ': Unknown fileformat.']);
end

if typ == 'P5'
  t = 5;
elseif typ == 'P2'
  t = 2;
else
  error(['Unknown fileformat ' typ '.'])
end

dm = [];

while length(dm) < 3

  [str,n] = fscanf(fid, '%s', 1);
  if n ~= 1
    error('Could not read header information.');
  end
  if str(1) ~= '#';
    a = str2num(str);
    if length(a) == 0
      error([ 'Strangeness in header: ' a '.' ]);
    end
    dm = [dm a];
  else % skip rest of line. Not quite sure what happens if empty comment
    c = 0;
    fprintf(1,'  Comment: ');
    while c ~= 10;
      c = fread(fid,1,'uchar');
      fprintf(1,'%c',c);
    end
  end

end

% Read pass the last newline
c = 0;
while c ~= 10;
  c = fread(fid,1,'uchar');
end


% [dm,c] = fscanf(fid, '%d %d %d', 3);

if  dm(1) < 1 | dm(2) < 1
  error(['Shitty dimensions of picture: ' dm(1) 'x' dm(2) '.'])
end

fprintf(1, '  Reading %s (%s ', fullname, typ);
fprintf(1, '%dx%d) - please wait.\n', [dm(1) dm(2)]);

% allocate image matrix

im = zeros(dm(2),dm(1));

% now read rows

for r=1:dm(2)
  if t == 5
    [row,c] = fread(fid, dm(1), 'uchar');
  else
    [row,c] = fscanf(fid, '%s ', dm(1));
  end

  if c ~= dm(1)
    fprintf(2,'Warning: Picture is truncated.');
  else
   im(r,:)=row';
  end

end  

fclose(fid);





