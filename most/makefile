TARGET=../most
OBJ=bwt.o common.o motif.o most.o FFT.o
CPP=g++
$(TARGET):$(OBJ)
	$(CPP) -o $(TARGET) $(OBJ)
	rm *.o
	@echo "success!"
motif.o:motif.cpp motif.h
	$(CPP) -c motif.cpp
common.o:common.h common.cpp
	$(CPP) -c common.cpp
most.o: bwt.o common.o motif.o
	$(CPP) -c most.cpp
bwt.o:bwt.h bwt.cpp
	$(CPP) -c bwt.cpp
FFT.o:FFT.h FFT.cpp
	$(CPP) -c FFT.cpp

