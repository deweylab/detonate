wd := $(shell pwd)

install/lib/libcityhash.a finished: cityhash-1.1.1.tar.gz
	tar xvf cityhash-1.1.1.tar.gz && cd cityhash-1.1.1 && ./configure --prefix=${wd}/install && make && make install && cd .. && touch finished

clean:
	rm -f finished && rm -rf install && rm -rf cityhash-1.1.1
