testBin: eec.o testBin.o simon_util_cpp/combinatorics.o toyjets/gaus.o toyjets/gen.o
	g++ $^ -o $@ -larmadillo

testCov: eec.o testCov.o simon_util_cpp/combinatorics.o toyjets/gaus.o toyjets/gen.o
	g++ $^ -o $@ -larmadillo

testTrans: eec.o testTrans.o simon_util_cpp/combinatorics.o toyjets/gaus.o toyjets/gen.o
	g++ $^ -o $@ -larmadillo

testWt: eec.o testWt.o simon_util_cpp/combinatorics.o toyjets/gaus.o
	g++ $^ -o $@ -larmadillo

testDR: eec.o testDR.o simon_util_cpp/combinatorics.o
	g++ $^ -o $@ -larmadillo

testComp: eec.o testComp.o simon_util_cpp/combinatorics.o
	g++ $^ -o $@ -larmadillo

%.o: %.cc
	g++ -c -o $@ $^ -I/usr/local/include/Minuit2 

clean: 
	rm -f *.o
	find . -maxdepth 1 -type f -executable -exec rm {} +
