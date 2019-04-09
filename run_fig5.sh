#
# Figure 5 incremental growth lines.
#

cd src
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/Fig5/esrf079.m', 'fig5'); exit"
cd ..
