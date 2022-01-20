clear all
close all

mot = 'abracadabra';

entropy (mot)

x(1) = length(findstr( 'a', mot )); % Nombre d'occurences de 'a'
x(2) = length(findstr( 'b', mot )); % Nombre d'occurences de 'b'
x(3) = length(findstr( 'c', mot )); % Nombre d'occurences de 'c'
x(4) = length(findstr( 'd', mot )); % Nombre d'occurences de 'c'
x(5) = length(findstr( 'r', mot )); % Nombre d'occurences de 'c'

f = x / length(mot);
[code,codelength,lmoy] = huffman(f);

mot_code = ''; % initialisation
for i=1:length(mot)
  switch mot(i)
   case 'a'
     mot_code = strcat( mot_code,code(1,:) );
   case 'b'
     mot_code = strcat( mot_code,code(2,:) );
   case 'c'
     mot_code = strcat( mot_code,code(3,:) );
   case 'd'
     mot_code = strcat( mot_code,code(4,:) );
   case 'r'
     mot_code = strcat( mot_code,code(5,:) );
  end
end
mot_code % affichage du mot code

mot_decode = ''; % initialisation
cur_pos = 1;
while cur_pos <= length(mot_code)
  if strncmp( mot_code(cur_pos:end), code(1,:), codelength(1));
    mot_decode = strcat( mot_decode,'a' ); cur_pos = cur_pos + codelength(1);
  elseif strncmp( mot_code(cur_pos:end), code(2,:), codelength(2));
    mot_decode = strcat( mot_decode,'b' ); cur_pos = cur_pos + codelength(2);
  elseif strncmp( mot_code(cur_pos:end), code(3,:), codelength(3));
    mot_decode = strcat( mot_decode,'c' ); cur_pos = cur_pos + codelength(3);
  elseif strncmp( mot_code(cur_pos:end), code(4,:), codelength(4));
    mot_decode = strcat( mot_decode,'d' ); cur_pos = cur_pos + codelength(4);
  elseif strncmp( mot_code(cur_pos:end), code(5,:), codelength(5));
    mot_decode = strcat( mot_decode,'r' ); cur_pos = cur_pos + codelength(5);
  else
    disp('erreur de decodage'); break;
  end
end
mot_decode % affichage du mot decode
