MOSDIR=$PWD/../build/
SOCODIR=$PWD/../build/cec2012/
CONFDIR=$PWD/../

for i in `seq 1 20`; do
    for j in `seq 1 25`; do
        LOGFILE=./res/logfile-f${i}-run${j}
        LOGOUT=./res/logfile-f${i}-run${j}.out
        env LD_LIBRARY_PATH=$MOSDIR/gaeda:$LD_LIBRARY_PATH $MOSDIR/gaedaexec/gaedaexec --techs-repository=$CONFDIR/techs_soco.cfg --cfg $CONFDIR/soco2010.cfg --problem-size=1000 --log-file $LOGFILE $SOCODIR/libcec2012_f${i}.so &> $LOGOUT
    done
done
