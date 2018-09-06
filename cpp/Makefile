# Makefile for SNAP compilation
# Author: Gustavo Peixoto de Oliveira

# compiler selection according to OS

UNAME = $(shell uname -s)

# Linux
ifeq ($(findstring Linux,${UNAME})),Linux)
 CXX = g++
endif

# Mac OS
ifeq ($(findstring Darwin,${UNAME})),Darwin)
 CXX = clang++
 #CXX = clang-omp++
endif

# flags
CXXFLAGS = -g -std=c++98 -Wall
LFLAGS = 

# directories
BASE_DIR = .

# check SNAP_DIR 
: "$SNAP_DIR:?Environmental variable SNAP_DIR is not set! Set it!"

# inclusions
INCLUDES += -I${BASE_DIR} 
INCLUDES += -I${SNAP_DIR}/snap-core
INCLUDES += -I${SNAP_DIR}/snap-exp
INCLUDES += -I${SNAP_DIR}/snap-adv
INCLUDES += -I${SNAP_DIR}/glib-core
INCLUDES += -I${SNAP_DIR}/glib-adv

# get list of all .cpp files in the current dir
SRC = $(wildcard ${BASE_DIR}/*.cpp) 
SRCH = $(wildcard ${BASE_DIR}/*.h) 

# substitution reference (create .o objects): 
OBJ = $(SRC:%.cpp=%.o)
OBJ2 = $(SNAP_DIR)/snap-core/Snap.o

#=== rules 

all: graphMetrics 
	@echo \
	 Executable compiled.

# === right
%.o: %.cpp ${SRCH}
	${CXX} ${CXXFLAGS} ${INCLUDES} -c $< -o $@ 

# executable
graphMetrics: ./main/mainGraphMetrics.o ${OBJ} ${OBJ2}
	${CXX} ${OBJ} ${OBJ2} ${CXXFLAGS} ${LFLAGS} $< -o $@

#=== phony variables

.PHONY: clean
#clean: 

# clean everything 
.PHONY: deepclean
deepclean:
	@echo
	@rm -f all
	@rm -f ${OBJ} 
	@find ${BASE_DIR} -name "?*.o" -exec rm {} \; 


# clean post-processing dir
.PHONY: cleanPP
cleanPP:
	@find ${PP_DIR} -type f -name "*.*" -exec rm -f {} \;
	@echo Post-processing dir cleaned.

# clean documentation dir
.PHONY: cleandoc
cleandoc:
	@rm -fr ${DOCS_DIR} 
	@echo Documentation dir cleaned.


#=== help 
# automatic variables
# $@ - the file name of the target of the rule
# $< - the name of the first prerequisite (here, source files)
#

# Makefile debugging tools
# - warnings:
# $(warning This is a message that can be placed everywhere in makefile) 
# - prints:
#   'make --just-print' or 'make -n'  
#   'make -p' (data-base)  
#   'make --warn-undefined-variables'
#   'make --debug'
