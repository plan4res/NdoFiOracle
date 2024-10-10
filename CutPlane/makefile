##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of CutPlane                                                     #
#                                                                            #
#   Input:  $(CC)        = compiler command                                  #
#           $(SW)        = compiler options                                  #
#           $(NDOSFH)    = the .h files to include for NDOSolver/FiOracle    #
#           $(NDOSFINC)  = the -I$( NDOSolver/FiOracle directories )         #
#           $(libOSIINC) = the -I$( include directories ) for libOSI         #
#           $(CutPlSDR)  = the directory where the source is                 #
#                                                                            #
#   Output: $(CutPlOBJ)  = the final object(s) / library                     #
#           $(CutPlLIB)  = external libraries + -L< libdirs >                #
#           $(CutPlH)    = the .h files to include for CutPlane              #
#           $(CutPlINC)  = the -I$( CutPlane directory )                     #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(CutPlSDR)/*.o $(CutPlSDR)/*~

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

CutPlOBJ = $(CutPlSDR)/CutPlane.o
CutPlLIB = 
CutPlH   = $(CutPlSDR)/CutPlane.h
CutPlINC = -I$(CutPlSDR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

$(CutPlOBJ): $(CutPlSDR)/CutPlane.C $(CutPlH) $(NDOSFH)
	$(CC) -c $*.C -o $@ $(NDOSFINC) $(libOSIINC) $(SW)

########################## End of makefile ###################################
