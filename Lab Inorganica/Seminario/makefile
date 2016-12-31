single.png: plotter.py simulated.dat
	python plotter.py

simulated.dat: magnet_field
	./magnet_field > simulated.dat

magnetic_field: magnet_field.c
	gcc -o magnet_field magnet_field.c
