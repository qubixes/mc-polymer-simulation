.PHONY: clean tar all gpupol analysis cpupol util gpupol2 denspol secstruct

RELEASE=0.9.1
PROJECT=conring
TAR_DIR=../tar
TAR_FILE="$(TAR_DIR)/$(PROJECT)-v$(RELEASE).tar.gz"
DIRS=cpupol util analysis denspol gpupol2 secstruct
export RELEASE
all: $(DIRS)

cpupol:
	$(MAKE) -C cpupol

analysis:
	$(MAKE) -C analysis

util:
	$(MAKE) -C util

denspol: 
	$(MAKE) -C denspol

gpupol2: 
	$(MAKE) -C gpupol2
	
secstruct:
	$(MAKE) -C secstruct


tar: clean
	if [ ! -d $(TAR_DIR) ]; then mkdir -p $(TAR_DIR); fi
	if [ ! -d $(TAR_DIR)/.old ]; then mkdir -p $(TAR_DIR)/.old; fi
	if [ -f $(TAR_DIR)/*.tar.gz ]; then mv $(TAR_DIR)/*.tar.gz $(TAR_DIR)/.old; fi
	tar --exclude data --exclude raw_data --exclude denspol/ee_topo.dat --exclude denspol/ee_topo.dat6 -czf $(TAR_FILE) .

install: all
	for dir in $(DIRS); do $(MAKE) -C $$dir install; done

clean: 
	rm -f $(EXECS) *.o *~ bin/*
	for dir in $(DIRS); do $(MAKE) -C $$dir clean; done
	make -C scripts clean
