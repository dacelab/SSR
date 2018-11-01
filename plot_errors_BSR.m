%==========================================================% 
% Monitor curves in cross-validation (block-SSR algorithm) %
%==========================================================%
clc;
clear all;
close all;
format long;

% Modify the name of directory for each test-case:
%  dirname = './Data/UCI_datasets/airfoil/'; 
dirname = './Data/synthetic_datasets/Friedman1/'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1 (CV with coarse grid and NO warm starts) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('*** STEP 1 (coarse lambda-grid, no warm starts) *** \n')

kfold = 5;

filename = strcat([dirname, 'lbd_grid_step1_bsr.dat']);
[fd,err] = fopen(filename,'r');
lbdgrid_c = fscanf(fd,'%f');
fclose(fd);

sz_c = length(lbdgrid_c);

filename = strcat([dirname, 'errorsAver_cv_step1_bsr.dat']);
[fd,err] = fopen(filename,'r');
errorAver_step1 = fscanf(fd,'%f');
fclose(fd);

filename = strcat([dirname, 'errors_cv_step1_bsr.dat']);
[fd,err] = fopen(filename,'r');
errors_step1 = fscanf(fd,'%f');
errors_step1 = reshape(errors_step1,sz_c,kfold);
fclose(fd);

indmin = find(errorAver_step1==min(errorAver_step1));
lbdopt = lbdgrid_c(sz_c-indmin+1)

figure(1)
if (kfold == 5)

 subplot(2,3,1)
 plot(log10(lbdgrid_c), fliplr(reshape(errors_step1(:,1),1,sz_c)), 'k--s')
 xlabel('log_{10}(\lambda)')
 title('split #1')
 
 subplot(2,3,2)
 plot(log10(lbdgrid_c), fliplr(reshape(errors_step1(:,2),1,sz_c)), 'k--s')
 xlabel('log_{10}(\lambda)')
 title('split #2')
 
 subplot(2,3,3)
 plot(log10(lbdgrid_c), fliplr(reshape(errors_step1(:,3),1,sz_c)), 'k--s')
 xlabel('log_{10}(\lambda)')
 title('split #3')
 
 subplot(2,3,4)
 plot(log10(lbdgrid_c), fliplr(reshape(errors_step1(:,4),1,sz_c)), 'k--s')
 xlabel('log_{10}(\lambda)')
 title('split #4')
 
 subplot(2,3,5)
 plot(log10(lbdgrid_c), fliplr(reshape(errors_step1(:,5),1,sz_c)), 'k--s')
 xlabel('log_{10}(\lambda)')
 title('split #5')
 
 subplot(2,3,6)
 plot(log10(lbdgrid_c), fliplr(reshape(errorAver_step1,1,sz_c)), '-r')
 xlabel('log_{10}(\lambda)')
 title('averaged')
 
else

 fprintf('\n Wrong value for kfold');
 fprintf('\n STOP.')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%STEP 2 (CV with fine grid and warm starts):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n')
fprintf('*** STEP 2 (fine lambda-grid, warm starts) *** \n')

filename = strcat([dirname, 'lbd_grid_step2_bsr.dat']);
[fd,err] = fopen(filename,'r');
lbdgrid_f = fscanf(fd,'%f');
fclose(fd);

sz_f = length(lbdgrid_f);

filename = strcat([dirname, 'errorsAver_cv_step2_bsr.dat']);
[fd,err] = fopen(filename,'r');
errorAver_step2 = fscanf(fd,'%f');
fclose(fd);

filename = strcat([dirname, 'errors_cv_step2_bsr.dat']);
[fd,err] = fopen(filename,'r');
errors_step2 = fscanf(fd,'%f');
errors_step2 = reshape(errors_step2,sz_f,kfold);
fclose(fd);

indmin = find(errorAver_step2==min(errorAver_step2));
lbdopt = lbdgrid_f(sz_f-indmin+1)

figure(2)
if (kfold == 5)

 subplot(2,3,1)
 plot(log10(lbdgrid_f), fliplr(reshape(errors_step2(:,1),1,sz_f)), 'k--s')
 xlabel('log_{10}(\lambda)')
 title('split #1')
 
 subplot(2,3,2)
 plot(log10(lbdgrid_f), fliplr(reshape(errors_step2(:,2),1,sz_f)), 'k--s')
 xlabel('log_{10}(\lambda)')
 title('split #2')
 
 subplot(2,3,3)
 plot(log10(lbdgrid_f), fliplr(reshape(errors_step2(:,3),1,sz_f)), 'k--s')
 xlabel('log_{10}(\lambda)')
 title('split #3')
 
 subplot(2,3,4)
 plot(log10(lbdgrid_f), fliplr(reshape(errors_step2(:,4),1,sz_f)), 'k--s')
 xlabel('log_{10}(\lambda)')
 title('split #4')
 
 subplot(2,3,5)
 plot(log10(lbdgrid_f), fliplr(reshape(errors_step2(:,5),1,sz_f)), 'k--s')
 xlabel('log_{10}(\lambda)')
 title('split #5')
 
 subplot(2,3,6)
 plot(log10(lbdgrid_f), fliplr(reshape(errorAver_step2,1,sz_f)), '-r')
 xlabel('log_{10}(\lambda)')
 title('averaged')
 
else

 fprintf('\n Wrong value for kfold');
 fprintf('\n STOP.')

end