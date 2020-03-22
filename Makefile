PFLAGS = `pythia8-config --cxxflags --libs`

RFLAGS = `root-config --cflags --libs`

CFLAGS = -Wall

PROGRAMS = lyz_pythia recursion_pythia

all: $(PROGRAMS)

clean:
	rm -f $(PROGRAMS)

debug: CFLAGS += -g
debug: $(PROGRAMS)

phenix: CFLAGS += -std=c++11 -m32
phenix: $(PROGRAMS)

lyz_pythia: lyz_pythia.C
	g++ -o lyz_pythia lyz_pythia.C $(CFLAGS) $(PFLAGS) $(RFLAGS)

recursion_pythia: recursion_pythia.C
	g++ -o recursion_pythia recursion_pythia.C $(CFLAGS) $(PFLAGS) $(RFLAGS)

