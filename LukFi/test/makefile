##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of a test main() for LukFiOrcl                                  #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# module name
NAME = test

# basic directory
DIR = .

# basic libraries
LIB = -lm

# debug switches
SW = -g -std=c++14 -ferror-limit=1
# production switches
#SW = -O3 -std=c++11 -DNDEBUG

# compiler
CC = clang++

# default target- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

default: $(NAME)

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(DIR)/*.o $(DIR)/*~ $(NAME)

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object(s) / library
# *LIB is the external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# libNDO
libNDOSDR = ../..
include $(libNDOSDR)/lib/makefile-inc

LukFiSDR = ..
include $(LukFiSDR)/makefile

# main module (linking phase) - - - - - - - - - - - - - - - - - - - - - - - -

# object files
MOBJ =  $(libNDOOBJ) $(LukFiOBJ)

# libraries
MLIB =  $(libNDOLIB) $(LukFiLIB)

$(NAME): $(MOBJ) $(DIR)/test.o
	$(CC) -o $(NAME) $(DIR)/test.o $(MOBJ) $(MLIB) $(SW)

# dependencies: every .o from its .cpp + every recursively included .h- - - -

# include dirs
MINC =  $(libNDOINC) $(LukFiINC)

# includes
MH =    $(libNDOH) $(LukFiH)

$(DIR)/test.o: $(DIR)/test.cpp $(MH)
	$(CC) -c $*.cpp -o $@ $(MINC) $(SW)

############################ End of makefile #################################
