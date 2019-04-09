#
# Figure 6B human diffusion limited with and without braces, and excess nutrition.
# Figure 6C human molar 3D slices.
#

cd src
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/Fig6BC/hum_dlm.m', 'fig6b_dlm'); exit"
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/Fig6BC/hum_dlm_braces.m', 'fig6b_dlm_braces'); exit"
nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('../parameters/2019/Fig6BC/hum_excess.m', 'fig6b_excess'); exit"
cd ..

# Generates the parameter files for each slice on the fly by appending the
# slice data file name to the beginning of the base parameter file.i
data='./data/2019/Fig6C'
for fname in $data/*.mat; do
    ffull=$(basename $fname)
    fbody="${ffull%.*}"
    par=tmp_$fbody.m
    echo "P.template = '".$data"/"$ffull"';" >> $par
    cat ./parameters/2019/Fig6BC/hum_dlm_braces_base.m >> $par 
    cd src
    nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('"../$par"', 'fig6c_"$fbody"'); exit"
    cd ..
    rm $par
done

