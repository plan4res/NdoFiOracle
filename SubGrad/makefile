##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of SubGrad                                                      #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(NDOSFH)    = the .h files to include for NDOSolver/FiOracle    #
#           $(NDOSFINC)  = the -I$( NDOSolver/FiOracle directories )         #
#           $(libCQKH)   = all the .h files of the CQK library               #
#           $(libCQKINC) = all the -I$( include dirs ) of the CQK library    #
#                          note: libCQK* are only needed if the macro        #
#                          SubGrad_HANDLES_CONSTRAINTS in SubGrad.h is != 0  #
#           $(SubGrSDR)  = the directory where the SubGrad source is         #
#                                                                            #
#   Output: $(SubGrOBJ)  = the final object(s) / library                     #
#           $(SubGrLIB)  = external libraries + -L< libdirs >                #
#           $(SubGrH)    = the .h files to include for SubGrad               #
#           $(SubGrDSH)  = the .h files to include for the Deflection and    #
#                          Stepsize components of SubGrad (which are not     #
#                          needed by the SubGrad object file)                #
#           $(SubGrINC)  = the -I( SubGrad directory )                       #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(SubGrSDR)/*.o $(SubGrSDR)/*~

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

SubGrOBJ = $(SubGrSDR)/SubGrad.o
SubGrLIB =
SubGrINC = -I$(SubGrSDR) -I$(SubGrSDR)/Deflection -I$(SubGrSDR)/Stepsize \
	$(libCQKINC)
SubGrH   = $(SubGrSDR)/SubGrad.h $(SubGrSDR)/Stepsize.h \
	$(SubGrSDR)/Deflection.h $(libCQKH)
SubGrSDH = $(SubGrSDR)/Deflection/PrimalDual.h \
	$(SubGrSDR)/Deflection/Volume.h \
	$(SubGrSDR)/Stepsize/ColorTV.h \
	$(SubGrSDR)/Stepsize/FumeroTV.h \
	$(SubGrSDR)/Stepsize/Polyak.h

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(SubGrOBJ): $(SubGrSDR)/SubGrad.C $(SubGrH) $(NDOSFH)
	$(CC) -c $*.C -o $@ $(NDOSFINC) $(libCQKINC) $(SW)

########################## End of makefile ###################################
