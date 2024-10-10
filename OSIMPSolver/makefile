##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of OSIMPSolver                                                  #
#                                                                            #
# Note: OSIMPSolver.C has a switch WHICH_OSI_MP that allows to set the       #
#       type of the OsiSolverInterface to either Cplex or GuRoBi. This is    #
#       an hack necessary to implement the quadratic stabilization, since    #
#       OsiSolverInterface does not allow to manage the quadratic objective  #
#       function necessary for proximal stabilization, thus methods of the   #
#       underlying solver's API have to be directly accessed. The downside   #
#       of this hack is that OSIMPSolver.C may then have a need top access   #
#       the corresponding include files (cplex.h or gurobi_c++.h), and       #
#       therefore the corresponding -I directives need be added to the       #
#       compile command. This is why $(libCPXINC) and $(libGRBINC) are in    #
#       the input list: only at most one of them is required, depending on   #
#       the value of WHICH_OSI_MP.                                           #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(NDOSFH)    = the .h files to include for NDOSolver/FiOracle    #
#           $(NDOSFINC)  = the -I$(NDOSolver/FiOracle directories)           #
#           $(MPSlvH)    = the .h files to include for MPSolver (only)       #
#           $(BndleINC)  = the -I$(Bundle/MPSolver directory)                #
#           $(libOSIINC) = the -I$( include directories ) for libOSI         #
#           $(libCPXINC) = the -I$( include directories ) for libCPX         #
#           $(libGRBINC) = the -I$( include directories ) for libGRB         #
#           $(OSIMPSDR)  = the directory where the OSIMPSolver source is     #
#                                                                            #
#   Output: $(OSIMPOBJ)  = the final object(s) / library                     #
#           $(OSIMPLIB)  = external libraries + -L< libdirs >                #
#           $(OSIMPH)    = the .h files to include                           #
#           $(OSIMPINC)  = the -I$( OSIMPSolver directory )                  #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(OSIMPSDR)/*.o $(OSIMPSDR)/*~

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

OSIMPOBJ = $(OSIMPSDR)/OSIMPSolver.o
OSIMPLIB = 
OSIMPH   = $(OSIMPSDR)/OSIMPSolver.h
OSIMPINC = -I$(OSIMPSDR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(OSIMPSDR)/OSIMPSolver.o: $(OSIMPSDR)/OSIMPSolver.C $(OSIMPH) $(MPSlvH) \
	$(NDOSFH)
	$(CC) -c $*.C -o $@ $(BndleINC) $(libOSIINC) $(libCPXINC) \
	$(libGRBINC) $(NDOSFINC) $(SW)

########################## End of makefile ###################################
