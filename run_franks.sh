#
# Frank sphere tests
#

cd src
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/frank_spheres/ref1.m', 'franks'); exit"
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/frank_spheres/ref2.m', 'franks'); exit"
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/frank_spheres/ref4.m', 'franks'); exit"
cd ..

