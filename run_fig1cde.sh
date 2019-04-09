#
# Figures 1C (diffusion-limited), 1D (excess nutrition) and 1E (no surface tension).
#

cd src
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/Fig1CDE/synthetic_dlm.m', 'fig1cde_dlm'); exit"
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/Fig1CDE/synthetic_excess.m', 'fig1cde_excess'); exit"
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/Fig1CDE/synthetic_notension.m', 'fig1cde_notension'); exit"
cd ..
