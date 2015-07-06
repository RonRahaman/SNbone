# makefile directives: make -j4 COMPILER=intel THREAD=yes PAPI=yes IBMHPM=yes OPT=3 USEMETIS=yes USEGREEDY=yes
# examples  : make -j4 COMPILER=gfortran
#           : make -j4 OPT=2 PAPI=no
#           : make -j4 OPT=1 DEBUG="-g -C"
#           : make -j4 OPT=1 DEBUG="-g -C -DLocal_Debug_FGMRES_Driver" THREAD=no
# -------------------
# -j4 specifies that 4 processes are available for parallel make 
# You can select a compiler by using make COMPILER=intel or COMPILER=gfortran or defining a machine name below
COMPILER=nvcc
# You can enable or disable OPENMP threading by using "THREAD=no"
# You can enable or disable the two PAPI options by using PAPIFLOP and PAPISET1. They cannot both be used
# You can change the optimization settings by using "OPT=2"
ARCH=sm_35
OPT=3
DEBUG=-g -G
FASTMATH=--use_fast_math
USEMETIS=no
USEGREEDY=no
METISLIB=

# Initialization of null variables
FF=echo
CC=echo
LD=echo
PAPIFLAGS=
PAPIINC=
PAPILIB=
IBMFLAGS=
IBMLIB=
# ------------------------------
# NVCC NVCC NVCC
# ------------------------------
CC = nvcc
CFLAGS := -c -arch ${ARCH} -O${OPT} ${FASTMATH} ${DEBUG}
LD = ${CC}
LDFLAGS = -lm -lcurand

# RR: not sure how these will be integrated
# ifeq (${USEGREEDY},yes)
#    EXTFLAGS=-DUSEGREEDY
#    METISLIB=
# endif
# ifeq (${USEMETIS},yes)
#    EXTFLAGS=-DUSEMETIS
#    METISLIB=../metis/libmetis.a
# endif

TARGET = cuSNaCFE.x

OBJS =  CommonBlock.o GenerateXb.o PrintSummary.o Verify.o Zero.o \
        Import_PMesh.o BuildAngleCubature.o AssembleNZmatrix.o StencilNZmatrix.o \
        SolveWGS.o FGMRES_Threaded.o \
        SolveWGS_PassThrough_AVE1.o \
        SolveWGS_PassThrough_AVE2.o \
        SolveWGS_PassThrough_AVE3.o \
        SolveWGS_PassThrough_PC.o \
        ApplyA_AVE1_Tet_SUPG.o ApplyA_AVE2_Tet_SUPG.o 
all: $(TARGET)

print_vars: 
	echo "host: " ${COMPILER}
	echo "comp: " ${CC} ${CFLAGS}
	echo "link: " ${LD} ${LDFLAGS}

BuildAngleCubature.o: BuildAngleCubature.cu CommonBlock.o
	${CC} ${CFLAGS} -o BuildAngleCubature.o BuildAngleCubature.cu
StencilNZmatrix.o: StencilNZmatrix.cu CommonBlock.o
	${CC} ${CFLAGS} -o StencilNZmatrix.o StencilNZmatrix.cu
AssembleNZmatrix.o: AssembleNZmatrix.cu CommonBlock.o
	${CC} ${CFLAGS} -o AssembleNZmatrix.o AssembleNZmatrix.cu

SolveWGS.o: SolveWGS.cu CommonBlock.o
	${CC} ${CFLAGS} -o SolveWGS.o SolveWGS.cu
FGMRES_Threaded.o: FGMRES_Threaded.cu CommonBlock.o 
	${CC} ${CFLAGS} -o FGMRES_Threaded.o FGMRES_Threaded.cu
Import_PMesh.o: Import_PMesh.cu CommonBlock.o
	${CC} ${CFLAGS} -o Import_PMesh.o Import_PMesh.cu

CommonBlock.o: CommonBlock.cu PROTEUS_Preprocess.h ApplyA_common.h
	${CC} ${CFLAGS} -o CommonBlock.o CommonBlock.cu

Verify.o: Verify.cu CommonBlock.o
	${CC} ${CFLAGS} -o Verify.o Verify.cu
Zero.o: Zero.cu
	${CC} ${CFLAGS} -o Zero.o Zero.cu
GenerateXb.o: GenerateXb.cu CommonBlock.o
	${CC} ${CFLAGS} -o GenerateXb.o GenerateXb.cu
PrintSummary.o: PrintSummary.cu CommonBlock.o
	${CC} ${CFLAGS} -o PrintSummary.o PrintSummary.cu

SolveWGS_PassThrough_AVE1.o: SolveWGS_PassThrough_AVE1.cu CommonBlock.o
	${CC} ${CFLAGS} -o SolveWGS_PassThrough_AVE1.o SolveWGS_PassThrough_AVE1.cu
SolveWGS_PassThrough_AVE2.o: SolveWGS_PassThrough_AVE2.cu CommonBlock.o
	${CC} ${CFLAGS} -o SolveWGS_PassThrough_AVE2.o SolveWGS_PassThrough_AVE2.cu

ApplyA_AVE1_Tet_SUPG.o: ApplyA_AVE1_Tet_SUPG.cu CommonBlock.o
	${CC} ${CFLAGS} -o ApplyA_AVE1_Tet_SUPG.o ApplyA_AVE1_Tet_SUPG.cu
ApplyA_AVE2_Tet_SUPG.o: ApplyA_AVE2_Tet_SUPG.cu CommonBlock.o
	${CC} ${CFLAGS} -o ApplyA_AVE2_Tet_SUPG.o ApplyA_AVE2_Tet_SUPG.cu
SolveWGS_PassThrough_AVE3.o: SolveWGS_PassThrough_AVE3.cu CommonBlock.o
	${CC} ${CFLAGS} -o SolveWGS_PassThrough_AVE3.o SolveWGS_PassThrough_AVE3.cu
SolveWGS_PassThrough_PC.o: SolveWGS_PassThrough_PC.cu CommonBlock.o
	${CC} ${CFLAGS} -o SolveWGS_PassThrough_PC.o SolveWGS_PassThrough_PC.cu

SetCommonFromFortran.o: SetCommonFromFortran.cu ApplyA_common.h
	${CC} ${CFLAGS} -o SetCommonFromFortran.o SetCommonFromFortran.cu

${TARGET}: ${OBJS} Driver.cu
	$(LD) ${LDFLAGS} -o ${TARGET} ${OBJS} Driver.cu ${METISLIB}

clean:
	rm -f ${OBJS} *.mod
