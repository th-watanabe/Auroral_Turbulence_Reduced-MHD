# for SX-aurora
FFC = mpinfort -compiler /opt/nec/ve/nfort/3.0.4/bin/nfort
###FFLAGS  = -report-all -fcheck=all -fpp $(EXTRA) #-mparallel
###FFLAGS  = -report-all -O3 -fpp $(EXTRA) #-mparallel
###FFLAGS  = -report-all -O1
###FFLAGS  = -report-all -O2
FFLAGS  = -report-all -O3
FFLAGS  += -fdiag-vector=2 -fdiag-inline=2 -fdiag-parallel=2
FFLAGS  += -proginf -ftrace

NLC_HOME=/opt/nec/ve/nlc/2.1.0
INC = -I$(NLC_HOME)/include
LIB = -L$(NLC_HOME)/lib  -laslfftw3 -lasl_sequential -ftrace
LIB += -static

SRC = ../src_SX/
MYL = ../src_SX/

## for FX100
#
## FFC = mpifrt
#FFC = mpifrtpx
#
#FFLAGS = -fw -Kfast,parallel,openmp -Kocl -Qt -Cpp -Koptmsg=2
#
## FFLAGS = -fw -Kfast,parallel -Kocl -Qt -Cpp -Haefosux

## FFC = mpifrtpx
## FFLAGS = -Kfast -Cpp
## FFLAGS = -Kfast,parallel -Cpp
## FFLAGS = -Kfast,parallel,openmp -Cpp

#LIB = -SSL2 -lfftw3 -lm -SSL2BLAMP

#SRC = ../src/
#MYL = ../src/

CWD = ./
OUT = rmhds.exe

MATH = rmhds_r0.1_math_SSL2

gkvp:	$(CWD)rmhds_r0.1_header.f90\
	$(SRC)rmhds_r0.1_mpienv.f90\
	$(MYL)$(MATH).f90\
	$(SRC)rmhds_r0.1_clock.f90\
	$(SRC)rmhds_r0.4_intgrl.f90\
	$(SRC)rmhds_r0.1_bndry.f90\
	$(SRC)rmhds_r0.5_pssn_fftwmm.f90\
	$(SRC)rmhds_r0.4_fld.f90\
	$(SRC)rmhds_r0.1_feedback.f90\
	$(SRC)rmhds_r0.4_set.f90\
	$(SRC)rmhds_r0.9_diag.f90\
	$(SRC)rmhds_r0.4_advnc.f90\
	$(SRC)rmhds_r0.1_main.f90
	
	$(FFC) $(FFLAGS) -c $(CWD)rmhds_r0.1_header.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhds_r0.1_mpienv.f90
	$(FFC) $(FFLAGS) -c $(MYL)$(MATH).f90 $(INC)
	$(FFC) $(FFLAGS) -c $(SRC)rmhds_r0.1_clock.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhds_r0.4_intgrl.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhds_r0.1_bndry.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhds_r0.5_pssn_fftwmm.f90 $(INC)
	$(FFC) $(FFLAGS) -c $(SRC)rmhds_r0.4_fld.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhds_r0.1_feedback.f90
	$(FFC) $(FFLAGS) -O1 -c $(SRC)rmhds_r0.4_set.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhds_r0.9_diag.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhds_r0.4_advnc.f90
	$(FFC) $(FFLAGS) -c $(SRC)rmhds_r0.1_main.f90
	
	$(FFC) $(FFLAGS) rmhds_r0.1_main.o\
			rmhds_r0.1_header.o\
			rmhds_r0.1_mpienv.o\
			$(MATH).o\
			rmhds_r0.1_clock.o\
			rmhds_r0.4_intgrl.o\
			rmhds_r0.1_bndry.o\
			rmhds_r0.5_pssn_fftwmm.o\
			rmhds_r0.4_fld.o\
			rmhds_r0.1_feedback.o\
			rmhds_r0.4_set.o\
			rmhds_r0.9_diag.o\
			rmhds_r0.4_advnc.o\
			-o $(OUT) $(LIB)

clean:
	rm -f ./*.o ./*.mod ./*.L ./*.exe ./*.err ./*.out

clear:
	rm -f ./*.o ./*.mod ./*.L ./*.err ./*.out
