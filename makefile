#CXX = icpc
CXX = clang++
#CXX = g++
#CXXFLAGS = -O2 -Wno-unused-result
#CXXFLAGS = -c -g -O2 -Wno-c++11-extensions
CXXFLAGS = -c -O2 -Wno-c++11-extensions -fno-cxx-exceptions


DEPENDENCIES = tricubic_utils.o tree_utils.o kernel_utils.o utilities.o


all:	direct_sum Tricubic 

#-----------------------------------------------------

direct_sum: direct_sum.o kernel_utils.o utilities.o
	$(CXX) direct_sum.o kernel_utils.o utilities.o -o direct_sum

direct_sum.o: direct_sum.cpp
	$(CXX) $(CXXFLAGS) direct_sum.cpp

#-----------------------------------------------------
Tricubic: Tricubic.o $(DEPENDENCIES)
	$(CXX) Tricubic.o $(DEPENDENCIES) -o Tricubic

Tricubic.o: Tricubic.cpp
	$(CXX) $(CXXFLAGS) Tricubic.cpp

#-----------------------------------------------------

tricubic_utils.o: tricubic_utils.cpp tricubic_utils.h tree_utils.o
	$(CXX) $(CXXFLAGS) tricubic_utils.cpp

tree_utils.o: tree_utils.cpp tree_utils.h kernel_utils.o
	$(CXX) $(CXXFLAGS) tree_utils.cpp

kernel_utils.o: kernel_utils.cpp kernel_utils.h utilities.o
	$(CXX) $(CXXFLAGS) kernel_utils.cpp

utilities.o: utilities.cpp utilities.h
	$(CXX) $(CXXFLAGS) utilities.cpp


#-----------------------------------------------------

clean:
	rm *.o *~

