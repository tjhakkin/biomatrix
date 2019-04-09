#
# Figure 2B pig molar, triconid and talonid separately.
#

cd src
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/Fig2B/triconid.m', 'fig2b_triconid'); exit"
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/Fig2B/talonid.m', 'fig2b_talonid'); exit"
cd ..
