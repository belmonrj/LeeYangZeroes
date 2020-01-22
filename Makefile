PFLAGS = `pythia8-config --cxxflags --libs`

RFLAGS = `root-config --cflags --libs`

PROGRAMS = lyz_pythia recursion_pythia

all: $(PROGRAMS)

clean:
	rm -f $(PROGRAMS)

lyz_pythia: lyz_pythia.C
	g++ -o lyz_pythia lyz_pythia.C -Wall $(PFLAGS) $(RFLAGS)

recursion_pythia: recursion_pythia.C
	g++ -o recursion_pythia recursion_pythia.C -Wall $(PFLAGS) $(RFLAGS)

