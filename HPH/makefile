#
# Set the variables for the working directory, module directory, compile
# options, and library locations
#
MQCDir       = $(mqcinstall)
GAUOPENDIR   = $(gauopenroot)
MQCMODS      = $(MQCDir)/PGI/mod
MQCLIB       = $(MQCDir)/PGI/lib
LIBS         = -llapack -lblas -L$(MQCLIB)
MQC_CPPFLAGS = -DMATRIXFILE -DPGI
#F03Flags    = -std=f2008
#RunF        = gfortran -fdefault-integer-8 -fdefault-real-8 
F03Flags     = 
RunF         = pgfortran -i8 -r8 -Mallocatable=03
#
#
# The 'all' rule.
#
all: fmo.exe

#
# Generic rules for building module (*.mod) and object (*.o) files.
#
%.f03: %.F03
	$(RunF) -F $(MQC_CPPFLAGS) $*.F03

%.mod: %.f03
	$(RunF) $(F03Flags) -c $*.f03

%.o: %.f03
	$(RunF) $(F03Flags) -I$(MQCMODS) -c $*.f03

#
# Generic rule for building general executable program (*.exe) from a standard
# Fortran 2003 source (*.f03) file.
#
%.exe: %.f03 $(MQCLIB)/libmqc.a
	$(RunF) $(LIBS) $(Prof) -I$(MQCMODS) -o $*.exe $*.f03 $(MQCLIB)/libmqc.a
