function varargout = MVcode(arg1, arg2, arg3)
% MVcode    Entropy encoding (or decoding) in a JPEG like manner
% Coding is very similar to sequential Huffman coding as described in
%   William B. Pennebaker, Joan L. Mitchell.
%   JPEG: Still Image Data Compression Standard, chapter 11.1
%   (Van Nostrand Reinhold, New York, USA, 1992, ISBN: 0442012721)
% The number of input arguments decide whether it is encoding or decoding,
% three input argument for decoding, two for encoding.
%
% [y,Res] = MVcode(Speed, W);              % encoding  
% y = MVcode(Speed, W);                    % encoding  
% W = MVcode(y, N, L);                     % decoding  
%----------------------------------------------------------------------
% arguments:
%  W        The motion vector matrix (x or y).
%  N        Number of rows in W
%  L        Number of columns in W, [N,L]=size(W)
%  Speed    For complete coding set Speed to 0. Set Speed to 1 to cheat 
%           during encoding, y will then be a sequence of zeros only,
%           but it will be of correct length and the other output 
%           arguments will be correct. 
%  y        a column vector of non-negative integers (bytes) representing 
%           the code, 0 <= y(i) <= 255. 
%  Res      The results (encoding only)        
%           Number og symbols:                 Res(1)            
%           Bits to store Huffman table (HL):  Res(2)            
%           Bits to store Huffman symbols:     Res(3)            
%           Bits to store additional bits:     Res(4)            
%           The total of bits is then:  sum(Res(2:4))
%----------------------------------------------------------------------

% Function needs following m-files: HuffLen, HuffCode, HuffTree

%----------------------------------------------------------------------
% Copyright (c) 1999.  Karl Skretting.  All rights reserved.
% Hogskolen in Stavanger (Stavanger University), Signal Processing Group
% Mail:  karl.skretting@tn.his.no   Homepage:  http://www.ux.his.no/~karlsk/
% 
% HISTORY:
% Ver. 1.0  26.07.99  Karl Skretting, Signal Processing Project 1999
%----------------------------------------------------------------------

global y Byte BitPos Speed 

Mfile='MVcode';
Debug=0;

% check input and output arguments, and assign values to arguments
if (nargin < 2); 
   error([Mfile,': function must have input arguments, see help.']); 
end
if (nargout < 1); 
   error([Mfile,': function must have output arguments, see help.']); 
end

if (nargin == 2)
   Encode=1;Decode=0;
   Speed=arg1;
   W=arg2(:)';
   clear arg1 arg2
   [N,L]=size(W);
   if ((length(Speed(:))~=1)); 
      error([Mfile,': Speed argument is not scalar, see help.']); 
   end
   if Speed; Speed=1; end;
elseif (nargin == 3)
   Encode=0;Decode=1;
   y=arg1(:);         % first argument is y
   N=arg2;
   L=arg3;
   clear arg1 arg2 arg3
   if ((length(N(:))~=1)); 
      error([Mfile,': N argument is not scalar, see help.']); 
   end
   if ((length(L(:))~=1)); 
      error([Mfile,': L argument is not scalar, see help.']); 
   end
else
   error([Mfile,': wrong number of arguments, see help.']); 
end

% Encode
if Encode
   Byte=0;BitPos=1;  % ready to write into first position
   Res=zeros(8,1);
   % take the DC component
   DC=[W(1,1),W(1,2:L)-W(1,1:(L-1))];
   DCsym=abs(DC);
   I=find(DCsym);
   DCsym(I)=floor(log2(DCsym(I)))+1;
   Hi=hist(DCsym,0:15);
   HL=HuffLen(Hi);
   % save HL in 8 byte, then DCsym and extra bits in bits
   bits=8*8+sum(HL.*Hi)+sum(DCsym);
   Res(1)=L;Res(2)=8*8;Res(3)=sum(HL.*Hi);Res(4)=sum(DCsym);
   y=zeros(ceil(bits/8),1);
   if Speed
      % advance Byte and BitPos without writing to y
      Byte=Byte+floor(bits/8);
      BitPos=BitPos-mod(bits,8);
      if (BitPos<=0); BitPos=BitPos+8; Byte=Byte+1; end;
   else
      % save HL 
      for j=1:16
         for (i=4:-1:1); PutBit(bitget(HL(j),i)); end;
      end
      HK=HuffCode(HL);
      % save DCsym and extra bits
      for l=1:L
         if Debug
            disp(['Encode:  DCsym(',int2str(l),')=',int2str(DCsym(l)),...
                  '   DC(',int2str(l),')=',int2str(DC(l))]);
         end
         n=DCsym(l)+1;    % symbol number (value 0 is first symbol, symbol 1)
         for k=1:HL(n)
            PutBit(HK(n,k));
         end
         if DCsym(l)        % extra bits
            if (DC(l)>0)    % the sign '1'='+', '0'='-'
               PutBit(1); 
               n=DC(l)-2^(DCsym(l)-1);
            else; 
               PutBit(0); 
               n=DC(l)+2^DCsym(l)-1;
            end;  
            % if Debug; disp(['   n=',int2str(n)]); end;
            for k=(DCsym(l)-1):-1:1
               PutBit(bitget(n,k));
            end
         end
      end
   end
  
   y=y(1:Byte);   
   varargout(1) = {y};
   if (nargout >= 2); varargout(2) = {Res}; end;

end

% Decode
if Decode
   W=zeros(N,L);
   Byte=0;BitPos=1;  % ready to read from first position
   % first read the HL tab
   HL=zeros(1,16);
   for j=1:16
      for (i=1:4); HL(j)=HL(j)*2+GetBit; end;
   end
   % then the symbols and extra bits
   DCsym=zeros(1,L);
   DC=DCsym;
   Htree=HuffTree(HL);
   root=1;pos=root;     
   l=0;  % number of symbols decoded so far
   while l<L
      if GetBit
         pos=Htree(pos,3);
      else
         pos=Htree(pos,2);
      end
      if Htree(pos,1)           % we have arrived at a leaf
         l=l+1;
         DCsym(l)=Htree(pos,2)-1;   % value is one less than symbol number
         pos=root;              % start at root again
         % is there extra bits
         if DCsym(l)
            DC(l)=2^(DCsym(l)-1);
            if (~(GetBit)); DC(l)=1-2*DC(l); end;
            if (DCsym(l)>1)
               n=0;
               for i=1:(DCsym(l)-1); n=n*2+GetBit; end;
               DC(l)=DC(l)+n;
            end
         end
         if Debug
            disp(['Decode:  DCsym(',int2str(l),')=',int2str(DCsym(l)),...
                  '   DC(',int2str(l),')=',int2str(DC(l))]);
         end
      end
   end
   W(1,1)=DC(1);
   for l=2:L; W(1,l)=W(1,l-1)+DC(l); end;
      
   varargout(1) = {W};
end

return;


% Functions to write and read a Bit
function PutBit(Bit)
global y Byte BitPos
BitPos=BitPos-1;
if (~BitPos); Byte=Byte+1; BitPos=8; end; 
y(Byte) = bitset(y(Byte),BitPos,Bit);
return
   
function Bit=GetBit
global y Byte BitPos
BitPos=BitPos-1;
if (~BitPos); Byte=Byte+1; BitPos=8; end; 
Bit=bitget(y(Byte),BitPos);
return;

% test of function
W=[-7,4,400,-2,0,1;1,0,1,0,-1,3;-1,1,0,0,-1,4;0,0,0,0,-7,8;-6,2,0,0,5,7]
[y,Res] = JPEGlike(0, W);
W2 = JPEGlike(y, 5, 6)


   