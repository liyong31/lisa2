all: T

T: 
	g++ swift.cc  spotutil.cc ltlf2fol.cc spotsynt.cc mona.cc -o swift -lspot -lbddx -lcudd -lmonadfa -lmonamem -lmonabdd -O3


#------------------------------------------------------
clean: #clean
	rm -f *.o main *.cc~ *.h~ Makefile~
#------------------------------------------------------

