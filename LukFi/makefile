##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of LukFi                                                        #
#                                                                            #
#   Input:  $(CC)       = compiler command                                   #
#           $(SW)       = compiler options                                   #
#           $(NDOFiH)   = the include files for FiOracle                     #
#           $(NDOSFINC) = the -I$( NDOSolver/FiOracle directories )          #
#           $(LukFiSDR) = the directory where the LukFiOrcl source is        #
#                                                                            #
#   Output: $(LukFiOBJ) = the final object(s) - library                      # 
#           $(LukFiLIB) = external libraries + -L< libdirs >                 #
#           $(LukFiH)   = the .h files to include                            #
#           $(LukFiINC) = the -I$( LukFiOrcl directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(LukFiSDR)/*.o $(LukFiSDR)/*~

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

LukFiOBJ = $(LukFiSDR)/LukFiOrcl.o
LukFiLIB =
LukFiH   = $(LukFiSDR)/LukFiOrcl.h 
LukFiINC = -I$(LukFiSDR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(LukFiOBJ): $(LukFiSDR)/LukFiOrcl.C $(LukFiH) $(NDOFiH)
	$(CC) -c $*.C -o $@ $(NDOSFINC) $(SW)

########################## End of makefile ###################################
