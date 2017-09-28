install:
	mkdir -p mafft/binaries mafft/scripts

	cd mafft/src && $(MAKE) $(MAKECMDGOALS)

	mkdir -p cgi-bin/mafft temp
	cp mafft/binaries/* mafft/scripts/*  cgi-bin/mafft/

clean:
	cd mafft/src && $(MAKE) $(MAKECMDGOALS)
	rm cgi-bin/mafft/*

