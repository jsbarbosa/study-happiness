plot.png : plotter.py functions.dat energies.dat
	python plotter.py
#	rm *.dat a.out

functions.dat energies.dat : a.out
	./a.out

a.out : seeker.c
	gcc -fopenmp seeker.c -lm
