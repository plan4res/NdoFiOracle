##############################################################################
################################ makefile ####################################
##############################################################################
#                                                                            #
#   makefile of MPTester                                                     #
#                                                                            #
#   Input:  $(MPTstSDR) = the directory where the MPTester source is         #
#                                                                            #
#   Output: $(MPTstH)   = the .h files to include for MPTester               #
#           $(MPTstINC) = the -I$(MPTester directory)                        #
#                                                                            #
#                              Antonio Frangioni                             #
#                         Dipartimento di Informatica                        #
#                             Universita' di Pisa                            #
#                                                                            #
##############################################################################

# clean - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clean::
	rm -f $(MPTstSDR)/*~

# macros to be exported - - - - - - - - - - - - - - - - - - - - - - - - - - -

MPTstH   = $(MPTstSDR)/MPTester.h
MPTstINC = -I$(MPTstSDR)

# dependencies: every .o from its .C + every recursively included .h- - - - -

########################## End of makefile ###################################
