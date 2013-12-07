build/lemon/libemon.a: tip.tar.gz
	tar xvf tip.tar.gz && mkdir build && cd build && cmake ../lemon-main-473c71baff72 && make lemon && touch finished

clean: tip.tar.gz
	rm -rf finished && rm -rf lemon-main-473c71baff72 && rm -rf build
