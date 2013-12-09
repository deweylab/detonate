build/lemon/libemon.a: tip.tar.gz
	tar xvf tip.tar.gz && mkdir -p build && cd build && cmake ../lemon-main-473c71baff72 && make lemon && cd .. && touch finished

clean: tip.tar.gz
	rm -f finished && rm -rf lemon-main-473c71baff72 && rm -rf build
