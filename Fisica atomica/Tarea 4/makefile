Tarea\ 4.pdf : Tarea\ 4.tex plot.pdf energies.pdf
	pdflatex -interaction=nonstopmode Tarea\ 4.tex

plot.pdf energies.pdf : plotter.py functions.dat energies.dat omegas.dat
	python plotter.py

functions.dat energies.dat omegas.dat : a.out
	./a.out

a.out : seeker.c
	gcc -fopenmp seeker.c -lm
clean : 
	rm *.dat a.out *.log *.out *.aux
