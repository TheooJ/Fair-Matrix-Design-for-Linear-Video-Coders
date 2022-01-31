clear all
close all

% x = [1 1 1 2 3 4 4 4 5 6 7 7];
% 
% a=x(x<inf);
%  
%  for i = 1:4
%    ind = randperm(numel(x), 1); % select one element out of numel(x) elements, with the probability of occurrence of the element in x
%    r(i) = x(ind);
%    x(ind) = []; % delete this element from the sample, such that the picked elements are unique
%  end
%  

 load N1
 load N2
% 
idx = N_canal1(N_canal1<inf);


 for i = 1:10
   ind = randperm(numel(idx), 1); % select one element out of numel(x) elements, with the probability of occurrence of the element in x
   N_vect_1(i) = idx(ind);
   idx(ind) = []; % delete this element from the sample, such that the picked elements are unique
 end
 
N_1 = diag(N_vect_1);
P_1=find_perm_inv(N_1);
N_1P = P_1*N_1*P_1';
