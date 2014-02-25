CC = gcc
PY = ${HOME}/local/include/python2.7
O = wfold.so
wWigIO:
	${CC}  -Wall -shared -fPIC -p -I ${PY} -I../H -I./ wRNA.c ../lib/libRNA.a -lm -o ${O}

clean:
	if [ -f ${O} ]; then rm ${O}; fi

