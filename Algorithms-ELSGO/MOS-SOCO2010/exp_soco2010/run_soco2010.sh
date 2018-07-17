MOSDIR=$PWD/../
SOCODIR=$PWD/../soco2010/
CONFDIR=$PWD/..

#for i in `seq 9 9`; do
#    for j in `seq 1 1`; do
#        LOGFILE=./res/logfile-f${i}-run${j}
#        LOGOUT=./res/logfile-f${i}-run${j}.out
#        env LD_LIBRARY_PATH=$MOSDIR/gaeda:$LD_LIBRARY_PATH $MOSDIR/gaedaexec/gaedaexec --techs-repository=$CONFDIR/techs_soco.cfg --cfg $CONFDIR/soco2010.cfg --problem-size=1000 --log-file $LOGFILE $SOCODIR/libf${i}.so &> $LOGOUT
#    done
#done
LOGFILE=./res/logfile-egg_problem
LOGOUT=./res/logfile-egg_problem.out
env LD_LIBRARY_PATH=$MOSDIR/gaeda:$LD_LIBRARY_PATH $MOSDIR/gaedaexec/gaedaexec --techs-repository=$CONFDIR/techs_soco.cfg --cfg $CONFDIR/soco2010.cfg --problem-size=1000 --log-file $LOGFILE $SOCODIR/libeeg_problem.so &> $LOGOUT