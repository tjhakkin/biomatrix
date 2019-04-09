#
# Figure 7 human, orangutan simulations.
#

cd src
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/Fig7/human.m', 'fig7_human'); exit"
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/Fig7/orangutan.m', 'fig7_orangutan'); exit"
cd ..
