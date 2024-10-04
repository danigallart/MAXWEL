# Variables
COMPILER = ifx
FLAGS = -g
EXEC = maxwel
OBJ_DIR = ./
F90_FILES := $(wildcard *.f90)
OBJFILES = $(patsubst %.f90, %.o, $(F90_FILES))

#Compile the program
$(EXEC) : $(OBJFILES)
	$(COMPILER) $(FLAGS) -o ./$(EXEC) $(OBJFILES)

#Dependencies

%.o: %.f90
	$(COMPILER) -c -o $@ $<

2D_harmonic_main.o: def_io.o def_variables.o def_vectors.o

make clean:
	rm $(OBJ_DIR)/*.o
	rm $(OBJ_DIR)/*.mod
	rm $(EXEC)