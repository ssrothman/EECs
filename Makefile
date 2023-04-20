LIBS="-larmadillo"
CFLAGS="-pg"

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

%.o: %.cc
	g++ -c -o $@ $^ -I/usr/local/include/Minuit2 $(CFLAGS)

clean: 
	rm -f *.o
	find . -maxdepth 1 -type f -executable -exec rm {} +
