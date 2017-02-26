temp.gif : boxes.py core.cpython-34m.so
	time python boxes.py

core.cpython-34m.so core.c : core.pyx setup.py
#	rm -r build core.c core.cpython-34m.so
	python setup.py build_ext --inplace
