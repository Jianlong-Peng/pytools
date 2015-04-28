CXX = g++
OBROOT = /home/xmluo/program/openbabel-2.3.2
CFLAGS = -Wall -O2 -I $(OBROOT)/include/openbabel-2.0
LDFLAGS = -L $(OBROOT)/lib -lopenbabel
OBJS = atom_type.o pmf.o

all: pmf_step1

%.o: %.cpp
	$(CXX) $(CFLAGS) -c $<

pmf_step1: pmf_step1.o $(OBJS)
	g++ $(LDFLAGS) pmf_step1.o $(OBJS) -o pmf_step1

#pmf_step2: pmf_step2.o $(OBJS)
#	g++ $(LDFLAGS) pmf_step2.o $(OBJS) -o pmf_step2

