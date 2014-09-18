all:
	$(MAKE) -C src
	mkdir -p bin
	mv src/qvz bin/qvz

clean:
	$(MAKE) -C src clean
	rm -f bin/qvz
