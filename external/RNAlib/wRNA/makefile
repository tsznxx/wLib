CC = gcc
CPP = g++
PY = ${HOME}/local/include/python2.7
O = wRNA.so

wRNA:
	${CPP} -Wall -shared -fPIC -p -I${PY} -I./ -I../plot -I../fold wRNA.cpp -lm ../plot/libplot.a ../fold/libfold.a  -o ${O}

clean:
	rm -rf *so

