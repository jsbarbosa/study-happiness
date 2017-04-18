all : plotter.py 1_data.dat
	python plotter.py

1_data.dat : a.out
	./a.out

a.out : atom.c
	gcc -lm atom.c

clean :
	rm *.dat a.out
	
