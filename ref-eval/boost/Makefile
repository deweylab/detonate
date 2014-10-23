wd := $(shell pwd)

finished:
	./bootstrap.sh && ./b2 && touch finished

clean:
	rm -f finished && ./b2 --clean
