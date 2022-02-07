for i in *.f90
do
  echo $i
  mkdir dir_$i
  cp $i dir_$i/
  pushd dir_$i
  sh $i
  rm $i
  sh temp.c 
  rm temp.c
  popd
done

for i in dir_*
do 
  mv $i `echo $i | sed -e 's/dir_//g' | sed -e 's/\.f90//g'`
done
