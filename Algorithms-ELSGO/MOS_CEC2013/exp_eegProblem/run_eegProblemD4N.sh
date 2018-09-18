MOSDIR=$PWD/../src/
SOCODIR=$PWD/../src/problems/EEG_BigOpt/
CONFDIR=$PWD/..

LOGFILE=./res/EEGProblemD4N
LOGOUT=./res/EEGProblemD4N.out
env LD_LIBRARY_PATH=$MOSDIR/gaeda:$LD_LIBRARY_PATH $MOSDIR/gaedaexec/gaedaexec --techs-repository=$CONFDIR/techs_cec2013.cfg --cfg $CONFDIR/cec2013.cfg --problem-size=1025 --log-file $LOGFILE $SOCODIR/libeeg_problem.so &> $LOGOUT
