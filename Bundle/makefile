##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of Bundle                                                       #
#                                                                            #
#   Input:  $(CC)       = compiler command                                   #
#           $(SW)       = compiler options                                   #
#           $(NDOSFH)   = the .h files to include for NDOSolver/FiOracle     #
#           $(NDOSFINC) = the -I$( NDOSolver/FiOracle directories )          #
#           $(BndleSDR) = the directory where the Bundle/MPSolver source is  #
#                                                                            #
#   Output: $(BndleOBJ) = the final object(s) / library                      #
#           $(BndleLIB) = external libraries + -L< libdirs >                 #
#           $(BndleH)   = the .h files to include for Bundle/MPSolver        #
#           $(MPSlvH)   = the .h files to include for MPSolver (only)        #
#           $(BndleINC) = the -I$( Bundle/MPSolver directory )               #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(BndleSDR)/*.o $(BndleSDR)/*~

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

BndleLIB =
BndleOBJ = $(BndleSDR)/Bundle.o
BndleINC = -I$(BndleSDR)
MPSlvH   = $(BndleSDR)/MPSolver.h
BndleH   = $(BndleSDR)/Bundle.h $(MPSlvH)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(BndleOBJ): $(BndleSDR)/Bundle.C $(MPSlvH) $(NDOSFH)
	$(CC) -c $*.C -o $@ $(NDOSFINC) $(SW)

########################## End of makefile ###################################
