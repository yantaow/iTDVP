LIBRARY_DIR=~/itensor3
APP=real_time
CCFILES=$(APP).cc util.cc imps.cc

#################################################################
#################################################################
#################################################################
#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/core.h
#Mappings --------------
OBJECTS=$(patsubst %.cc,%.o, $(CCFILES))
GOBJECTS=$(patsubst %,.debug_objs/%, $(OBJECTS))

#Rules ------------------
%.o: %.cc $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -Wno-reorder -o $@ $<

.debug_objs/%.o: %.cc  $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -Wno-reorder -o $@ $<

#Targets -----------------

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(APP) $(LIBFLAGS) -lgsl -larmadillo -g -fopenmp -Wno-reorder 

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS) 
clean:
	rm -fr .debug_objs *.o $(APP) $(APP)-g

mkdebugdir:
	mkdir -p .debug_objs

