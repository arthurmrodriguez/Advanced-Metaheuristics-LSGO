CC=g++

all:
	cd cec2013 && $(MAKE)
	$(CC) -o dg2 mainEEG.cpp DG2.cpp Function.cpp CEC2013.cpp -I./cec2013 -L./cec2013 -lcec2013
	cp -r cec2013/cdatafiles cec2013/libcec2013.so .

test:
	g++ -g DG2.cpp Function.cpp test.cpp TestFun.cpp

clean:
	rm -fr dg2 libcec2013.so cdatafiles
	cd cec2013 && $(MAKE) clean
