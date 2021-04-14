HOMEL=$(PWD)
TOP=../bin
BIN=$(HOMEL)/../bin

override LDLIBS += -lgsl -lopenblas -lm -lnick
# Some Linux distributions require separate lapacke library
# override LDLIBS += -llapacke
override LDFLAGS += -g  -L./nicksrc
override CFLAGS += -c -g -p -pg -Wimplicit -I./ -I./nicksrc
# Mac additions using homebrew installations
#override CFLAGS += -I/usr/local/opt/openblas/include -I/usr/local/opt/gsl/include
#override LDFLAGS += -L/usr/local/opt/openblas/lib -L/usr/local/opt/gsl/lib
# Harvard Medical School O2 cluster additions
ifdef SLURM_CONF
override CFLAGS += -I/n/app/openblas/0.2.19/include -I/n/app/gsl/2.3/include
override LDFLAGS += -L/n/app/openblas/0.2.19/lib -L/n/app/gsl/2.3/lib/
TOP=../bin
BIN=$(HOMEL)/../bin
endif

ND = nicksrc
NLIB = $(ND)/libnick.a

PROGS= gsmooth refinedate  

CC = cc

all: $(NLIB) $(PROGS) 


$(NLIB):
	$(MAKE) -C $(ND)

## an ubuntu user found he needed: make all LDLIBS="-llapack" 

statsubs.o:     tables nicksrc/statsubs.c
	$(CC)  $(CFLAGS) -o statsubs.o nicksrc/statsubs.c

install:	all 
	mkdir -p $(TOP)
	cp $(PROGS) $(TOP)

uninstall:
	rm -f $(TOP)/*

dirs:	
	mkdir -p  $(HOMEL)/admixtables

clean: 
	rm -f *.o 
	rm -f *junk*
	rm -f core
	$(MAKE) -C $(ND) clean
	rm -f $(PROGS)

clobber: clean uninstall



## gsl-config may help here


gsmooth: $(NLIB) gsmooth.o 

refinedate: $(NLIB) refinedate.o 


