LIBS=-larmadillo
CFLAGS=-O3

testTrans_oo: testTrans_oo.o simon_util_cpp/combinatorics.o toyjets/gaus.o toyjets/gen.o jetinfo.o maxDR.o compositions.o
	g++ $^ -o $@ $(LIBS) $(CFLAGS)

testBin: eec.o testBin.o simon_util_cpp/combinatorics.o toyjets/gaus.o toyjets/gen.o
	g++ $^ -o $@ $(LIBS) $(CFLAGS)

testCov: eec.o testCov.o simon_util_cpp/combinatorics.o toyjets/gaus.o toyjets/gen.o
	g++ $^ -o $@ $(LIBS) $(CFLAGS)

testTrans: eec.o testTrans.o simon_util_cpp/combinatorics.o toyjets/gaus.o toyjets/gen.o
	g++ $^ -o $@ $(LIBS) $(CFLAGS)

testWt: eec.o testWt.o simon_util_cpp/combinatorics.o toyjets/gaus.o
	g++ $^ -o $@ $(LIBS) $(CFLAGS)

testDR: eec.o testDR.o simon_util_cpp/combinatorics.o
	g++ $^ -o $@ $(LIBS) $(CFLAGS)

testComp: eec.o testComp.o simon_util_cpp/combinatorics.o
	g++ $^ -o $@ $(LIBS) $(CFLAGS)

testTrans_oo.o : testTrans_oo.cc eec_oo.h jetinfo.h compositions.h maxDR.h adj.h
	g++ -c -o $@ testTrans_oo.cc -I/usr/local/include/Minuit2 $(CFLAGS)

%.o: %.cc
	g++ -c -o $@ $^ -I/usr/local/include/Minuit2 $(CFLAGS)

clean: 
	rm -f *.o
	find . -maxdepth 1 -type f -executable -exec rm {} +
