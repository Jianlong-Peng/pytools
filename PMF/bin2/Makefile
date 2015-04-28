CXX = g++
OB_DIR = /home/xmluo/program/openbabel-2.3.2
CFLAGS = -Wall -O2 -I $(OB_DIR)/include/openbabel-2.0 -fPIC
LDFLAGS = -L $(OB_DIR)/lib -lopenbabel

all: pmf_atom_pairs libatom_type.so

atom_type.o: atom_type.cpp
	$(CXX) $(CFLAGS) $< -c

pmf_atom_pairs.o: pmf_atom_pairs.cpp
	$(CXX) $(CFLAGS) $< -c

pmf_atom_pairs: pmf_atom_pairs.o atom_type.o
	$(CXX) $(LDFLAGS) pmf_atom_pairs.o atom_type.o -o pmf_atom_pairs

libatom_type.so: atom_type.o
	$(CXX) $(LDFLAGS) atom_type.o -shared -o libatom_type.so

