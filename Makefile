CPP_SRCS = tsp-ga.cc

BIN = tsp-ga

all:
	$(CXX) -g tsp-ga.cc -o ${BIN}

clean:
	rm -f ${BIN}



