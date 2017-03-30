all : plotter.py data.dat
	python plotter.py

data.dat : a.out
	./a.out
a.out : atom.c
	gcc -lm atom.c
	
