function [] = pgmsave(m,name,cmt)
%function b = pgmsave(m,name,cmt)
%pgmsave - save a matrix in PGM format.
%------------------------------------------------------------------------------
%SYNOPSIS       pgmsave(M,'STR')
%                 Floor values in [0,256[ and save them as a PGM P5
%                 image. STR is the wanted filename,
%
%               pgmsave(M,'STR','COMMENT')
%                 As above but also writes a comment in the file.
%
%SEE ALSO       pgmload, the pbmplus package.
%
%------------------------------------------------------------------------------
%Jonas Svanberg                        (C) 1994 Image Coding Group. LiU, SWEDEN

%RCSID          $Id: pgmsave.m,v 1.1 1997/12/11 15:33:25 svan Exp harna $



if nargin < 2
  error('Not enough input arguments.')
end

if ~isstr(name)
  error('Second arg. must be a string.')
end

fid = fopen(name,'w');

if fid == -1
  error(['Can''t open ' name '.'])
end

% print the header

if nargin < 3
  fprintf(fid, 'P5\n%d %d\n255\n', fliplr(size(m)));
else
  fprintf(fid, 'P5\n# %s\n%d %d\n255\n', cmt, fliplr(size(m)));
end

% now restrict and save the data

m(m>255)=255;
m(m<0)=0;

c=fwrite(fid, m', 'uint8');

if c ~= prod(size(m))
  error('Error writing to file.')
end

fclose(fid);
