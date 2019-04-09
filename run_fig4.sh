#
# Figure 4 pig cusp 3D slices.
#

# Generates the parameter files for each slice on the fly by appending the
# slice data file name to the beginning of the base parameter file.
data='./data/2019/Fig4'
for fname in $data/*.mat; do
    ffull=$(basename $fname)
    fbody="${ffull%.*}"
    par=tmp_$fbody.m
    echo "P.template = '".$data"/"$ffull"';" >> $par
    cat ./parameters/2019/Fig4/pig_entoconid_cusp_base.m >> $par
    cd src 
    nice -n 19 matlab -nosplash -nodesktop -r "main_matlab('"../$par"', 'fig4_"$fbody"'); exit"
    cd ..
    rm $par
done

