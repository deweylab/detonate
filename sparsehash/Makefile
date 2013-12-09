wd := $(shell pwd)

include/sparsehash/internal/sparseconfig.h finished: sparsehash-2.0.2.tar.gz
	tar xvf sparsehash-2.0.2.tar.gz
	cd sparsehash-2.0.2 && ./configure && make src/sparsehash/internal/sparseconfig.h
	mkdir include
	mkdir include/sparsehash
	cp -rp sparsehash-2.0.2/src/sparsehash/* include/sparsehash/
	touch finished

clean:
	rm -f finished && rm -rf include && rm -rf sparsehash-2.0.2
