.PHONY: clean tar all gpupol analysis cpupol util gpupol2

RELEASE=0.6.0
PROJECT=conring
TAR_DIR=../tar
TAR_FILE="$(TAR_DIR)/$(PROJECT)-v$(RELEASE).tar.gz"
DIRS=cpupol util analysis gpupol2

all: $(DIRS)

cpupol:
	$(MAKE) -C cpupol

analysis:
	$(MAKE) -C analysis

util:
	$(MAKE) -C util

gpupol2: 
	$(MAKE) -C gpupol2

tar: clean
	if [ ! -d $(TAR_DIR) ]; then mkdir -p $(TAR_DIR); fi
	if [ ! -d $(TAR_DIR)/.old ]; then mkdir -p $(TAR_DIR)/.old; fi
	if [ -f $(TAR_DIR)/*.tar.gz ]; then mv $(TAR_DIR)/*.tar.gz $(TAR_DIR)/.old; fi
	tar --exclude data --exclude raw_data -czf $(TAR_FILE) .

install: all
	for dir in $(DIRS); do $(MAKE) -C $$dir install; done

clean: 
	rm -f $(EXECS) *.o *~ bin/*
	for dir in $(DIRS); do $(MAKE) -C $$dir clean; done
