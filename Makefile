testCov: eec.cc testCov.cc simon_util_cpp/combinatorics.cc toyjets/gaus.cc toyjets/gen.cc
	g++ $^ -o $@ -larmadillo

testTrans: eec.cc testTrans.cc simon_util_cpp/combinatorics.cc toyjets/gaus.cc toyjets/gen.cc
	g++ $^ -o $@ -larmadillo

testWt: eec.cc testWt.cc simon_util_cpp/combinatorics.cc toyjets/gaus.cc
	g++ $^ -o $@

testDR: eec.cc testDR.cc simon_util_cpp/combinatorics.cc
	g++ $^ -o $@

testComp: eec.cc testComp.cc simon_util_cpp/combinatorics.cc
	g++ $^ -o $@

clean: 
	rm -f *.o
	find . -maxdepth 1 -type f -executable -exec rm {} +
