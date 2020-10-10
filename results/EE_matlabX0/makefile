defaults: EE_main
CFLAGS                          =
FFLAGS                          =
CPPFLAGS                        =
FPPFLAGS                        =

BIN_DIR=.

OBJS = \
${BIN_DIR}/EE_main.o\
${BIN_DIR}/ReadCsvData.o\
${BIN_DIR}/ea.o\
${BIN_DIR}/NormalizeData.o\
${BIN_DIR}/ee.o\
${BIN_DIR}/other_fun.o\

${BIN_DIR}/%.o : %.c myhead.h

EE_main: ${OBJS}
	-${CC} -o ${BIN_DIR}/EE_main ${OBJS}
	${RM} ${OBJS}

clear:
	-${RM} ${BIN_DIR}/EE_main

run:
	${BIN_DIR}/EE_main hspcs_pca.csv hspcs_matlab_X0.csv hspcs_matlabX0_ee
