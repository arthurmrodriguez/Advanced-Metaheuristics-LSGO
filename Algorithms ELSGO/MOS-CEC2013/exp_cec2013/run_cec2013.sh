MOSDIR=$PWD/../build/
SOCODIR=$PWD/../build/problems/cec2013_lsgo/
CONFDIR=$PWD/../

for i in `seq 8 15`; do
    for j in `seq 1 25`; do
        LOGFILE=./res/logfile-f${i}-run${j}
        LOGOUT=./res/logfile-f${i}-run${j}.out
        if [ $i -eq 13 -o $i -eq 14 ]; then
            env LD_LIBRARY_PATH=$MOSDIR/gaeda:$HOME/local/packages/hdf5-1.8.10-patch1/lib/:$LD_LIBRARY_PATH $MOSDIR/gaedaexec/gaedaexec --techs-repository=$CONFDIR/techs_cec2013.cfg --cfg $CONFDIR/cec2013.cfg --problem-size=905 --additional-data $SOCODIR/datafiles/ --log-file $LOGFILE $SOCODIR/libcec2013_f${i}.so &> $LOGOUT
        else
            env LD_LIBRARY_PATH=$MOSDIR/gaeda:$HOME/local/packages/hdf5-1.8.10-patch1/lib/:$LD_LIBRARY_PATH $MOSDIR/gaedaexec/gaedaexec --techs-repository=$CONFDIR/techs_cec2013.cfg --cfg $CONFDIR/cec2013.cfg --problem-size=1000 --additional-data $SOCODIR/datafiles/ --log-file $LOGFILE $SOCODIR/libcec2013_f${i}.so &> $LOGOUT
        fi
    done
done
