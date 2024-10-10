##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of the NDOSolver/FiOracle interface pair                        #
#                                                                            #
#   Input:  $(NDOSFSDR) = the directory where the NDOSolver/FiOracle source  #
#                         is                                                 #
#           $(OPTUxH)   = the OPT****.h files to include                     #
#           $(OPTUxINC) = the -I$( OPT****.h directory )                     #
#                                                                            #
#   Output: $(NDOSFH)   = the .h files to include for NDOSolver/FiOracle     #
#           $(NDOFiH)   = the .h files to include for FiOracle (only)        #
#           $(NDOSFINC) = the -I$( NDOSolver/FiOracle directory )            #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(NDOSFSDR)/*~

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

NDOFiH   = $(NDOSFSDR)/FiOracle.h $(OPTUxH)
NDOSFH   = $(NDOSFSDR)/NDOSlver.h $(NDOFiH)
NDOSFINC = -I$(NDOSFSDR) $(OPTUxINC)

# dependencies: every .o from its .C + every recursively included .h- - - - -

########################## End of makefile ###################################
