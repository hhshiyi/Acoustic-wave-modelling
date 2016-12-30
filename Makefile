modelling : modelling.o
	gcc -o modelling modelling.o -lm -lGL -lglut

modelling.o : modelling.c
	gcc -c modelling.c

clean: 
	rm *.o modelling
