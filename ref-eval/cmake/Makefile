wd := $(shell pwd)

install/bin/cmake finished: cmake-2.8.12.1.tar.gz
	tar xvf cmake-2.8.12.1.tar.gz && cd cmake-2.8.12.1 && ./bootstrap --prefix=${wd}/install && make && make install && cd .. && touch finished

clean:
	rm -f finished && rm -rf install && rm -rf cmake-2.8.12.1
