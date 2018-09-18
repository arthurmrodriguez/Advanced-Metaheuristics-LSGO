# Create the virtualenv
DIR=venv

if [ ! -d $DIR ]; then
	echo "Creating virtualenv..."
	python -m venv $DIR

	if [ $? ]; then
    		python -m venv $DIR --without-pip
    		source $DIR/bin/activate
    		export PATH=$PWD/$DIR/python:$PATH
    		curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
    		python get-pip.py
	else
    		source $DIR/bin/activate
	fi

	# Install dependencies
	pip install cython
	pip install -r requirements.txt
	cd cec2013lsgo
	python setup.py install
	cd ..
	# Install locally package
	cd ea
	python setup.py install
	cd ..
else
	source $DIR/bin/activate
	# This allows the script to be re-run in case
	# recompile is needed due to EEG Problem adaptation
	cd cec2013lsgo
	python setup.py install
	cd ..
fi

sed -i -E "s?/home/dmolina/shadeils/shadeils-env?${PWD}/${DIR}?" shadeils.py
