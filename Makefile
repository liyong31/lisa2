all: T

T: 
	g++ --std=c++17 lisa.cc  spotutil.cc ltlf2fol.cc mona.cc -o lisa2 -lspot -lbddx -lmonadfa -lmonamem -lmonabdd -O3

#	g++ --std=c++17 swift.cc  spotutil.cc ltlf2fol.cc spotsynt.cc mona.cc -o swift -lspot -lbddx -lmonadfa -lmonamem -lmonabdd -O3


#------------------------------------------------------
clean: #clean
	rm -f *.o main *.cc~ *.h~ Makefile~
#------------------------------------------------------

