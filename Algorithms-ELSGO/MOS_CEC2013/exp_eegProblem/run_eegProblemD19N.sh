MOSDIR=$PWD/../src/
SOCODIR=$PWD/../src/problems/EEG_BigOpt/
CONFDIR=$PWD/..

LOGFILE=./res/EEGProblemD19N
LOGOUT=./res/EEGProblemD19N.out
env LD_LIBRARY_PATH=$MOSDIR/gaeda:$LD_LIBRARY_PATH $MOSDIR/gaedaexec/gaedaexec --techs-repository=$CONFDIR/techs_cec2013.cfg --cfg $CONFDIR/cec2013.cfg --problem-size=4865 --log-file $LOGFILE $SOCODIR/libeeg_problem.so &> $LOGOUT
