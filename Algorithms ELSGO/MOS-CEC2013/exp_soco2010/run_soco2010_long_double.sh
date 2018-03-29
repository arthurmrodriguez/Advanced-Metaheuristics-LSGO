MOSDIR=$PWD/../build2/
SOCODIR=$PWD/../build2/problems/soco2010/
CONFDIR=$PWD/../

for i in `seq 7 7`; do
    for j in `seq 1 25`; do
        LOGFILE=./res.long_double/logfile-f${i}-run${j}
        LOGOUT=./res.long_double/logfile-f${i}-run${j}.out
        env LD_LIBRARY_PATH=$MOSDIR/gaeda:$LD_LIBRARY_PATH $MOSDIR/gaedaexec/gaedaexec --techs-repository=$CONFDIR/techs_cec2013.cfg --cfg $CONFDIR/cec2013.cfg --problem-size 1000 --log-file $LOGFILE $SOCODIR/libf${i}.so &> $LOGOUT
    done
done

for i in `seq 15 15`; do
    for j in `seq 1 25`; do
        LOGFILE=./res.long_double/logfile-f${i}-run${j}
        LOGOUT=./res.long_double/logfile-f${i}-run${j}.out
        env LD_LIBRARY_PATH=$MOSDIR/gaeda:$LD_LIBRARY_PATH $MOSDIR/gaedaexec/gaedaexec --techs-repository=$CONFDIR/techs_cec2013.cfg --cfg $CONFDIR/cec2013.cfg --problem-size 1000 --log-file $LOGFILE $SOCODIR/libf${i}.so &> $LOGOUT
    done
done
