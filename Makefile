help:
	@echo "Usage: make {all rall clean} "
	@echo "make all    -- compile/run all C examples for feast_sparse interface "
	@echo "make rall   -- compile and run all C examples for feast_sparse interface "
	@echo "make clean  -- clean all C examples for feast_sparse interface "
	@echo "!!!!Please correct accordingly compiler and libraries paths, change compiler options " 
	@echo " in file ../../make.inc !!!!"
	@echo
CC=icc
CFLAGS=-O3
INCL=-I/home/quantum/FEAST/3.0/include
LIB=-L/home/quantum/FEAST/3.0/lib/x64 -lfeast_sparse -lfeast -mkl

#==============================================================
# Include the LIB (feast and  pardiso- lapack - blas)  
#==============================================================
#LIB = $(LOCLIBS) $(FEAST_SPARSE) $(FEAST) $(CLIBS)
#==============================================================
# List of codes to be compiled 
#==============================================================
EXAMPLES = main
#==============================================================
# Compile Link Execute
#==============================================================
all: examples 


examples: 
	@echo $(EXAMPLES)
	@for file in $(EXAMPLES); \
	do \
		echo $(CC)  $(CFLAGS) $(INCL) -c $$file.c;\
		$(CC)  $(CFLAGS) $(INCL) -c $$file.c;\
		echo $(CC)   -o $$file $$file.o $(LIB) ;\
		$(CC)   -o $$file $$file.o $(LIB) ;\
		echo rm $$file.o;\
		rm $$file.o ;\
	done


rall: 	examples
	@for file in $(EXAMPLES); \
	do \
                ./$$file; \
	done
#==========================================
# Clean up directory: delete object files 
#==========================================
clean: 
	-@rm  $(EXAMPLES) *.o
