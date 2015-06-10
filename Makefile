INC=-I/usr/include/hdf
LIBS=-lmfhdf -ldf -lm
CXX=g++
LD=g++
CXXFLAGS=-g -O2 -Wall $(INC)
LDFLAGS=$(LIBS) -lopencv_core
TARG=modisresam
OFILES=\
	utils.o\
	resample.o\
	main.o\
	readwrite.o\
	allocate_2d.o\

HFILES=\
	modisresam.h\

all: $(TARG)

$(TARG): $(OFILES)
	$(LD) -o $(TARG) $(OFILES) $(LDFLAGS)

%.o: %.cc $(HFILES)
	$(CXX) $(CXXFLAGS) -c $<

sort.h: scripts/gensortindices.py
	python scripts/gensortindices.py > sort.h

install: $(TARG)
	cp $(TARG) /usr/local/bin/

clean:
	rm -f $(OFILES) $(TARG)
