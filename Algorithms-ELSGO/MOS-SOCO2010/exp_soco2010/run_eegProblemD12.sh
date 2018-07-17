
MOSDIR=$PWD/../
SOCODIR=$PWD/../soco2010/
CONFDIR=$PWD/..

LOGFILE=./res/logfile-egg_problemD12
LOGOUT=./res/logfile-egg_problemD12.out
env LD_LIBRARY_PATH=$MOSDIR/gaeda:$LD_LIBRARY_PATH $MOSDIR/gaedaexec/gaedaexec --techs-repository=$CONFDIR/techs_soco.cfg --cfg $CONFDIR/soco2010.cfg --problem-size=3072 --log-file $LOGFILE $SOCODIR/libeeg_problem.so &> $LOGOUT