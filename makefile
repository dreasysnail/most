TARGET=../most
OBJ=bwt.o common.o motif.o most.o FFT.o
CPP=g++
$(TARGET):$(OBJ)
	$(CPP) -g -o $(TARGET) $(OBJ)
	rm *.o
	@echo "success!"
motif.o:motif.cpp motif.h
	$(CPP) -g -c motif.cpp
common.o:common.h common.cpp
	$(CPP) -g -c common.cpp
most.o: bwt.o common.o motif.o
	$(CPP) -g -c most.cpp
bwt.o:bwt.h bwt.cpp
	$(CPP) -g -c bwt.cpp
FFT.o:FFT.h FFT.cpp
	$(CPP) -g -c FFT.cpp

