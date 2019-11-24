#this is the procedure I use to update to newer versions of boost in gatb-core
#pretty simple but gets the job done
#to be run within thirdparty/
#-Rayan

newdir=boost_1_71_0/boost/
olddir=boost

for file in `ls $olddir`
do
    echo $file
    cp -R $newdir/$file $olddir/
done
