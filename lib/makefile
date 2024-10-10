##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of libNDO                                                       #
#                                                                            #
#   all modules are given in input                                           #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(OPTUxH)    = the OPT****.h files to include                    #
#           $(OPTUxINC)  = the -I$( OPT****.h directory )                    #
#           $(libOSIINC) = the -I$(include directories) for libOSI           #
#           $(libNDOSDR) = the basic directory for the library               #
#                                                                            #
#   Output: $(libNDOOBJ) = the libNDO library (that must be built)           #
#           $(libNDOLIB) = the external libraries + -L< libdirs >            #
#           $(libNDOH)   = all the .h files to include                       #
#           $(libNDOINC) = all the -I$(include directories)                  #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# define & include the necessary modules- - - - - - - - - - - - - - - - - - -
# if a module is not used in the current configuration, just comment out the
# corresponding include line
# each module outputs some macros to be used here:
# *OBJ is the final object/library
# *LIB external libraries + -L< libdirs >
# *H   is the list of all include files
# *INC is the -I< include directories >

# NDOSlver
NDOSFSDR = $(libNDOSDR)/NDOSlver
include $(NDOSFSDR)/makefile

# Bundle
BndleSDR = $(libNDOSDR)/Bundle
include $(BndleSDR)/makefile

# CutPlane
#CutPlSDR = $(libNDOSDR)/CutPlane
#include $(CutPlSDR)/makefile

# MPTester
#MPTstSDR = $(libNDOSDR)/MPTester
#include $(MPTstSDR)/makefile

# OSIMPSolver
OSIMPSDR = $(libNDOSDR)/OSIMPSolver
include $(OSIMPSDR)/makefile

# QPPnltMP
QPPMPSDR = $(libNDOSDR)/QPPnltMP
include $(QPPMPSDR)/makefile

# SubGrad
# if SubGrad_HANDLES_CONSTRAINTS != 0, rather include makefile-c so as to
# have the libCQK stuff added
SubGrSDR = $(libNDOSDR)/SubGrad
include $(SubGrSDR)/makefile

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

# object files
libNDOOBJ =	$(BndleOBJ)\
		$(CutPlOBJ)\
		$(OSIMPOBJ)\
		$(QPPMPOBJ)\
		$(SubGrOBJ)

# libraries
libNDOLIB =	$(BndleLIB)\
		$(CutPlLIB)\
		$(OSIMPLIB)\
		$(QPPMPLIB)\
		$(SubGrLIB)

# includes
libNDOH   =	$(NDOSFH)\
		$(BndleH)\
		$(CutPlH)\
		$(MPTstH)\
		$(OSIMPH)\
		$(QPPMPH)\
		$(SubGrH)\
		$(SubGrDSH)

# include dirs
libNDOINC =	$(NDOSFINC)\
		$(BndleINC)\
		$(CutPlINC)\
		$(MPTstINC)\
		$(OSIMPINC)\
		$(QPPMPINC)\
		$(SubGrINC)

############################ End of makefile #################################
