.SUFFIXES: .c .cpp .o .ex


  CC =mpic++  -fopenmp  -ggdb  -Wall  -std=gnu++17  -Wfatal-errors 



INCLUDE= -I /public1/soft/eigen/3.4.0/include/eigen3

.cpp.o:
	$(CC) $(INCLUDE) -c $< 

.c.o:
	gcc $(INCLUDE) -c $< 

.o.ex:
	@echo g++ ... -o $@ $< ... $(OBJS) ... $(LIBS) 
	@$(CC) -o $@ $< $(OBJS) $(LIBS)

######################################################################
# targets
######################################################################
OBJS  = 

objs: $(OBJS)

clean:
	rm -f *.o *.ex *~

