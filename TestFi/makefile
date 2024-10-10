##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of TestFi                                                       #
#                                                                            #
#   Input:  $(CC)       = compiler command                                   #
#           $(SW)       = compiler options                                   #
#           $(NDOFiH)   = the include files for FiOracle                     #
#           $(NDOSFINC) = the -I$( NDOSolver/FiOracle directories )          #
#           $(TstFiSDR) = the directory where the TestFi source is           #
#                                                                            #
#   Output: $(TstFiOBJ) = the final object(s) - library                      #
#           $(TstFiLIB) = external libraries + -L< libdirs >                 #
#           $(TstFiH)   = the .h files to include                            #
#           $(TstFiINC) = the -I$( TestFi directory )                        #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(TstFiSDR)/*.o $(TstFiSDR)/*~

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

TstFiOBJ = $(TstFiSDR)/TestFi.o
TstFiLIB =
TstFiH   = $(TstFiSDR)/TestFi.h 
TstFiINC = -I$(TstFiSDR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(TstFiOBJ): $(TstFiSDR)/TestFi.C $(TstFiH) $(NDOFiH)
	$(CC) -c $*.C -o $@ $(NDOSFINC) $(SW)

########################## End of makefile ###################################
