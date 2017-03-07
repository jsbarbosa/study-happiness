temp : test.py temp.cpython-34m.so
	python test.py
temp.cpython-34m.so temp.c : temp.pyx setup.py
#	rm -r build core.c core.cpython-34m.so
	python setup.py build_ext --inplace
