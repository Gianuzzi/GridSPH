FC = gfortran -g
MYFFLAGS = -ffinite-math-only -funsafe-math-optimizations -ffast-math -funroll-loops -fopenmp #-fdefault-real-8
FFLAGS = -Wall -Wextra -march=native -fcheck=all $(MYFFLAGS)
MYLDFLAGS = -O3
LDFLAGS = $(MYLDFLAGS)

TARGETS = main
FIXED_SOURCES = $(wildcard *.f)
FREE_SOURCES = $(wildcard *.f90)
FIXED_OBJECTS = $(patsubst %.f,%.o,${FIXED_SOURCES})
FREE_OBJECTS = $(patsubst %.f90,%.o,${FREE_SOURCES})

MAKEDEP = deps
DEP_FILE_NAME = main.dep
DEP_FILE = $(DEP_FILE_NAME)

all: $(TARGETS) $(DEP_FILE_NAME)

${FIXED_OBJECTS} : %.o : %.f
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

${FREE_OBJECTS} : %.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(DEP_FILE_NAME) $(LDFLAGS)
	
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< $(LDFLAGS)

main: ${FIXED_OBJECTS}
	$(FC) $(FFLAGS) -o main $^ $(LDFLAGS)

deps:
	@echo "Creating dependencies file: $(DEP_FILE_NAME)"
	fortdepend -w -o $(DEP_FILE_NAME) -f *.f90 -i omp_lib

install: # Install fortdepend
	@echo "Installing fortdepend"
	@python -m pip install fortdepend
	make deps

clean:
	@echo "rm -rf *.o *.mod main"
	@rm -rf *.o *.mod main

# Include (or not) and create (if necessary) main.dep
ifeq ($(MAKECMDGOALS), clean)
else ifeq ($(MAKECMDGOALS), install)
else ifeq ($(MAKECMDGOALS), deps)
else
ifneq ("$(wildcard $(DEP_FILE_NAME))","")
-include $(DEP_FILE_NAME)
else
-include $(DEP_FILE)
endif
endif

## The creation of main.dep is made with the python package: fortdepend
## It can be installed via pip: 
### $ python -m pip install fortdepend
${DEP_FILE}:
	fortdepend -w -o $(DEP_FILE_NAME) -f *.f90 -i omp_lib

.PHONY: clean all deps install
