CXX       := "/home/bartunix/MulticoreBSP-for-C/tools/bspcc"
CXX_FLAGS := -O3

BIN     := bin
SRC     := src
INCLUDE := include
INCLUDE_BSP := "/home/bartunix/MulticoreBSP-for-C/"

LIBRARIES   := 
EXECUTABLE  := main

all: ${BIN}/${EXECUTABLE}

run: clean all
	clear
	./${BIN}/${EXECUTABLE}

${BIN}/${EXECUTABLE}: ${SRC}/*.c
	${CXX} ${CXX_FLAGS} -I${INCLUDE} -I${INCLUDE_BSP} $^ -o $@ ${LIBRARIES} -lm

clean:
	-rm ${BIN}/*
