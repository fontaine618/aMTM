all: libguiledistribLOH.so distribLOH_wrap.c

%_wrap.c: %.i
	swig -guile -Linkage passive $<

libguile%.so: %.c %_wrap.c
	gcc --std=gnu99 -fpic $$CFLAGS `gsl-config --cflags` -shared $^ $$LDFLAGS `gsl-config --libs` -o $@
