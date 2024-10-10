##############################################################################
################################# makefile ###################################
##############################################################################
#                                                                            #
#   makefile of QPPnltMP                                                     #
#                                                                            #
#   Input:  $(CC)       = compiler command                                   #
#           $(SW)       = compiler options                                   #
#           $(OPTUxH)   = the OPT****.h files to include                     #
#           $(OPTUxINC) = the -I$( OPT****.h directory )                     #
#   	    $(NDOSFH)   = the .h files to include for NDOSolver/FiOracle     #
#           $(NDOSFINC) = the -I$(NDOSolver/FiOracle directories)            #
#           $(MPSlvH)   = the .h files to include for MPSolver (only)	     #
#           $(BndleINC) = the -I$( Bundle/MPSolver directory )               #
#           $(QPPMPSDR) = the directory where the QPPnltMP source is         #
#                                                                            #
#   Output: $(QPPMPOBJ) = the final object(s) / library                      #
#           $(QPPMPLIB) = external libraries + -L< libdirs >                 #
#           $(QPPMPH)   = the .h files to include                            #
#           $(QPPMPINC) = the -I$( QPPnltMP directory )                      #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(QPPMPSDR)/*.o $(QPPMPSDR)/*~

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

# comment out BMinQuad.o to build the unconstrained version
QPPMPOBJ = $(QPPMPSDR)/QPPnltMP.o $(QPPMPSDR)/MinQuad.o \
	$(QPPMPSDR)/BMinQuad.o
QPPMPLIB =
# comment out BMinQuad.h to build the unconstrained version
QPPMPH   = $(QPPMPSDR)/QPPnltMP.h $(QPPMPSDR)/MinQuad.h \
	$(QPPMPSDR)/BMinQuad.h
QPPMPINC = -I$(QPPMPSDR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(QPPMPSDR)/MinQuad.o: $(QPPMPSDR)/MinQuad.C $(QPPMPSDR)/MinQuad.h $(OPTUxH)
	$(CC) -c $*.C -o $@ $(OPTUxINC) $(SW)

# comment the following three lines to build the unconstrained version
$(QPPMPSDR)/BMinQuad.o: $(QPPMPSDR)/BMinQuad.C $(QPPMPSDR)/BMinQuad.h \
	$(QPPMPSDR)/MinQuad.h $(OPTUxH)
	$(CC) -c $*.C -o $@ $(OPTUxINC) $(SW)

$(QPPMPSDR)/QPPnltMP.o: $(QPPMPSDR)/QPPnltMP.C $(QPPMPH) $(MPSlvH) $(NDOSFH)
	$(CC) -c $*.C -o $@ $(NDOSFINC) $(BndleINC) $(SW)

########################### End of makefile ##################################
