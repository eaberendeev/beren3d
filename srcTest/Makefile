CXX = g++
source_dirs := .
EIGEN_INCLUDE_PATH=-isystem ../../../eigen/

search_wildcards       := $(addsuffix /*.cpp,$(source_dirs)) 
override compile_flags += -Wall -pipe -std=c++11 -O3 

test: $(notdir $(patsubst %.cpp,%.o,$(wildcard $(search_wildcards))))
	$(CXX) $^ $(compile_flags) -o $@

VPATH := $(source_dirs)
     
%.o: %.cpp
	$(CXX) -c -MD $(compile_flags) $(EIGEN_INCLUDE_PATH) $(addprefix -I,$(source_dirs))  $^

clean:
	rm -f test *.o *.d
