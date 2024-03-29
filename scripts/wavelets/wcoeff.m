function [h,g,hh,gg,desp]= wcoeff(code)

%   Filter coefficients for wavelet transforms
%   University of Santiago, http://web.usc.es/~elusive

%   code = 'haar', 'daub4', 'daub6', 'daub8', .., 'daub20',
%          'filter9-7', 'adelson', 'haar2-6', 'haar5-3', 'haar9-3',
%          'villa13-11', 'villa6-10', 'villa10-18', 'odegar'
%          'brislawn'

%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2, or (at your option)
%   any later version.
 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
 
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software Foundation,
%   Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  */

%------------------------------------------------------------*/
%   References:

%   Geoff Davis, Baseline Wavelet Transform Coder Construction Kit
%   gdavis@cs.dartmouth.edu, http://www.cs.dartmouth.edu/~gdavis 
  
%   The Bath Wavelet Warehouse
%   http://dmsun4.bath.ac.uk/wavelets/warehouse.html
   
%   Bank of Wavelet Filters
%   http://www.isds.duke.edu/~brani/filters.html */
%------------------------------------------------------------*/

% Haar and Daubechies filters  
%   I. Daubechies, "Ten Lectures on Wavelets" */ 

hh=[];
desp=0;

if strcmp(code, 'haar')

   h=[ 0.7071067811865475,      0.7071067811865475 ];

elseif strcmp(code, 'daub4')

   h=[ 0.4829629131445341,      0.8365163037378077, ...
       0.2241438680420134,     -0.1294095225512603 ];

elseif strcmp(code, 'daub6')

   h=[ 0.3326705529500825,      0.8068915093110924, ...
       0.4598775021184914,     -0.1350110200102546, ...
      -0.0854412738820267,      0.0352262918857095 ];

elseif strcmp(code, 'daub8')

   h=[ 0.2303778133088964,      0.7148465705529154, ...
       0.6308807679398587,     -0.0279837694168599, ...
      -0.1870348117190931,      0.0308413818355607, ...
       0.0328830116668852,     -0.0105974017850690 ];

elseif strcmp(code, 'daub10')
   
   h=[ 0.1601023979741929,      0.6038292697971895, ...  
       0.7243085284377726,      0.1384281459013203, ...  
      -0.2422948870663823,     -0.0322448695846381, ...   
       0.0775714938400459,     -0.0062414902127983, ...  
      -0.0125807519990820,      0.0033357252854738 ];

elseif strcmp(code, 'daub12')
   
   h=[ 0.1115407433501095,      0.4946238903984533, ...
       0.7511339080210959,      0.3152503517091982, ...
      -0.2262646939654400,     -0.1297668675672625, ...
       0.0975016055873225,      0.0275228655303053, ...
      -0.0315820393174862,      0.0005538422011614, ... 
       0.0047772575109455,     -0.0010773010853085 ]; 
  
elseif strcmp(code, 'daub14')
   
   h=[ 0.0778520540850037,      0.3965393194818912, ...
       0.7291320908461957,      0.4697822874051889, ...
      -0.1439060039285212,     -0.2240361849938412, ...
       0.0713092192668272,      0.0806126091510774, ...
      -0.0380299369350104,     -0.0165745416306655, ...
       0.0125509985560986,      0.0004295779729214, ...
      -0.0018016407040473,      0.0003537137999745 ];
 
elseif strcmp(code, 'daub16')
   
   h=[ 0.0544158422431072,      0.3128715909143166, ...
       0.6756307362973195,      0.5853546836542159, ...
      -0.0158291052563823,     -0.2840155429615824, ...
       0.0004724845739124,      0.1287474266204893, ...
      -0.0173693010018090,     -0.0440882539307971, ...
       0.0139810279174001,      0.0087460940474065, ...
      -0.0048703529934520,     -0.0003917403733770, ...
       0.0006754494064506,     -0.0001174767841248 ];
 
elseif strcmp(code, 'daub18')
 
   h=[ 0.0380779473638778,      0.2438346746125858, ...
       0.6048231236900955,      0.6572880780512736, ...
       0.1331973858249883,     -0.2932737832791663, ...
      -0.0968407832229492,      0.1485407493381256, ...
       0.0307256814793385,     -0.0676328290613279, ...
       0.0002509471148340,      0.0223616621236798, ...
      -0.0047232047577518,     -0.0042815036824635, ...
       0.0018476468830563,      0.0002303857635232, ...
      -0.0002519631889427,      0.0000393473203163 ];
 
elseif strcmp(code, 'daub20')
 
   h=[ 0.0266700579005473,      0.1881768000776347, ...
       0.5272011889315757,      0.6884590394534363, ...
       0.2811723436605715,     -0.2498464243271598, ...
      -0.1959462743772862,      0.1273693403357541, ...
       0.0930573646035547,     -0.0713941471663501, ...
      -0.0294575368218399,      0.0332126740593612, ...
       0.0036065535669870,     -0.0107331754833007, ...
       0.0013953517470688,      0.0019924052951925, ...
      -0.0006858566949564,     -0.0001164668551285, ...
       0.0000935886703202,     -0.0000132642028945 ];

%------------------------------------------------------------*/

%   9/7 Filter from M. Antonini, M. Barlaud, P. Mathieu, and
%   I. Daubechies, "Image coding using wavelet transform", IEEE
%   Transactions on Image Processing", Vol. pp. 205-220, 1992. */

%   J.N. Bradley, C.M. Brislawn, and T. Hopper, "The FBI Wavelet/Scalar 
%   Quantization Standard for Gray-scale Fingerprint Image Compression", 
%   SPIE Proceedings vol 1961, Visual Information Processing II,
%   Orlando, Florida, pp. 293-304, Apr 1993. */

elseif strcmp(code, 'filter9-7')
  
   h=[ 3.782845550699535e-02,  -2.384946501937986e-02, ...
      -1.106244044184226e-01,   3.774028556126536e-01, ...
       8.526986790094022e-01,   3.774028556126537e-01, ...
      -1.106244044184226e-01,  -2.384946501937986e-02, ...
       3.782845550699535e-02 ];

  hh=[-6.453888262893856e-02,  -4.068941760955867e-02, ...
       4.180922732222124e-01,   7.884856164056651e-01, ...
       4.180922732222124e-01,  -4.068941760955867e-02, ...
      -6.453888262893856e-02 ];

%------------------------------------------------------------*/

%   Filter from Eero Simoncelli's PhD thesis -- 
%   used in Edward Adelson's EPIC wavelet coder
%   These are probably the filter coefficients used in Shapiro's EZW paper */

elseif strcmp(code, 'adelson')

   h=[ 0.028220367, -0.060394127, -0.07388188, ...
       0.41394752,   0.7984298,    0.41394752, ... 
      -0.07388188,  -0.060394127,  0.028220367, 0 ];

  hh=[ 0,                                      ...
       0.028220367, -0.060394127, -0.07388188, ...
       0.41394752,   0.7984298,    0.41394752, ... 
      -0.07388188,  -0.060394127,  0.028220367 ];

%------------------------------------------------------------*/

%   Filters evaluated in: J. Villasenor, B. Belzer, J. Liao, 
%   "Wavelet Filter Evaluation for Image Compression." IEEE Transactions 
%   on Image Processing, Vol. 2, pp. 1053-1060, August 1995. */

%   In Haar-like filters, a single normalization factor can be chosen 
%   to make all coefficients integers. */

elseif strcmp(code, 'haar2-6')

   h=[ 7.071067811865476e-01,   7.071067811865476e-01 ];

  hh=[-8.838834764831845e-02,   8.838834764831845e-02, ...
       7.071067811865476e-01,   7.071067811865476e-01, ...
       8.838834764831845e-02,  -8.838834764831845e-02 ];
       
  desp=-2;

elseif strcmp(code, 'haar5-3')

   h=[-1.767766952966369e-01,   3.535533905932738e-01, ...
       1.060660171779821e+00,   3.535533905932738e-01, ...
      -1.767766952966369e-01 ];

  hh=[ 3.535533905932738e-01,   7.071067811865476e-01, ...
       3.535533905932738e-01 ];

elseif strcmp(code, 'haar9-3')

   h=[ 3.314563036811943e-02,  -6.629126073623885e-02, ...
      -1.767766952966369e-01,   4.198446513295127e-01, ...
       9.943689110435828e-01,   4.198446513295127e-01, ...
      -1.767766952966369e-01,  -6.629126073623885e-02, ...
       3.314563036811943e-02 ];

  hh=[ 3.535533905932738e-01,   7.071067811865476e-01, ...
       3.535533905932738e-01 ];
       
  desp=2;

elseif strcmp(code, 'villa13-11')

   h=[-8.472827741318157e-03,   3.759210316686883e-03, ...
       4.728175282882753e-02,  -3.347508104780150e-02, ...
      -6.887811419061032e-02,   3.832692613243884e-01, ...
       7.672451593927493e-01,   3.832692613243889e-01, ...
      -6.887811419061045e-02,  -3.347508104780156e-02, ...
       4.728175282882753e-02,   3.759210316686883e-03, ... 
      -8.472827741318157e-03 ];

  hh=[ 1.418215589126359e-02,   6.292315666859828e-03, ...
      -1.087373652243805e-01,  -6.916271012030040e-02, ...
       4.481085999263908e-01,   8.328475700934288e-01, ...
       4.481085999263908e-01,  -6.916271012030040e-02, ...
      -1.087373652243805e-01,   6.292315666859828e-03, ...
       1.418215589126359e-02 ];

elseif strcmp(code, 'villa6-10')

   h=[-1.290777652578771e-01,   4.769893003875977e-02, ...
       7.884856164056651e-01,   7.884856164056651e-01, ...
       4.769893003875977e-02,  -1.290777652578771e-01 ];

  hh=[ 1.891422775349768e-02,   6.989495243807747e-03, ... 
      -6.723693471890128e-02,   1.333892255971154e-01, ... 
       6.150507673110278e-01,   6.150507673110278e-01, ...
       1.333892255971154e-01,  -6.723693471890128e-02, ...
       6.989495243807747e-03,   1.891422775349768e-02 ];

  desp=-2;

% Unpublished 10/18 filter from Villasenor's group */

elseif strcmp(code, 'villa10-18')

   h=[ 2.885256501123136e-02,   8.244478227504624e-05, ...
      -1.575264469076351e-01,   7.679048884691438e-02, ...
       7.589077294537618e-01,   7.589077294537619e-01, ...
       7.679048884691436e-02,  -1.575264469076351e-01, ...
       8.244478227504624e-05,   2.885256501123136e-02];

  hh=[ 9.544158682436510e-04,  -2.727196296995984e-06, ...
      -9.452462998353147e-03,  -2.528037293949898e-03, ...
       3.083373438534281e-02,  -1.376513483818621e-02, ...
      -8.566118833165798e-02,   1.633685405569902e-01, ...
       6.233596410344172e-01,   6.233596410344158e-01, ...
       1.633685405569888e-01,  -8.566118833165885e-02, ...
      -1.376513483818652e-02,   3.083373438534267e-02, ...
      -2.528037293949898e-03,  -9.452462998353147e-03, ...
      -2.727196296995984e-06,   9.544158682436510e-04];

  desp= 4;

%------------------------------------------------------------*/

% Odegar Filter */

elseif strcmp(code, 'odegar')

   h=[ 5.2865768532960523e-02, -3.3418473279346828e-02, ...
      -9.3069263703582719e-02,  3.8697186387262039e-01, ...
       7.8751377152779212e-01,  3.8697186387262039e-01, ...
      -9.3069263703582719e-02, -3.3418473279346828e-02, ...
       5.2865768532960523e-02 ];
  
  hh=[-8.6748316131711606e-02, -5.4836926902779436e-02, ...
       4.4030170672498536e-01,  8.1678063499210640e-01, ...
       4.4030170672498536e-01, -5.4836926902779436e-02, ...
      -8.6748316131711606e-02 ];

%------------------------------------------------------------*/

% Filters from Chris Brislawn's tutorial code */

elseif strcmp(code, 'brislawn')

   h=[ 0.026913419, -0.032303352, ...
      -0.241109818,  0.054100420, ...
       0.899506092,  0.899506092, ... 
       0.054100420, -0.241109818, ... 
      -0.032303352,  0.026913419];

  hh=[ 0.019843545,  0.023817599, ... 
      -0.023257840,  0.145570740, ...  
       0.541132748,  0.541132748, ... 
       0.145570740, -0.023257840, ... 
       0.023817599,  0.019843545 ];

%------------------------------------------------------------*/

else error('NOT DEFINED FILTER');
end

if(size(hh)==[0 0]) 
   M1=size(h,2);
   hh(: ,[1:M1])=h([M1:-1:1]);
end

M1=size(h,2);M2=size(hh,2); 
if(rem(M1,2)==1) h(M1+1)=0; M1=M1+1; end
if(rem(M2,2)==1) hh(M2+1)=0;M2=M2+1; end
g(:,[1:2:M2])=hh([1:2:M2]);   g(:,[2:2:M2])=-hh([2:2:M2]);
gg(:,[1:2:M1])= -h([1:2:M1]); gg(:,[2:2:M1])=  h([2:2:M1]);
