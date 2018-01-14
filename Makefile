CC = gcc-7

LIBS  = gsl gslcblas m 
CCFLAGS += -g -Wall -std=gnu99 -fmax-errors=5 #-Werror 
#CCFLAGS += -ffast-math -O2 -ftree-vectorize 

OBJS = IO_ass2.o

all : IO_ass2 stud_t pet rj # $(OBJS)  

# Ecc_SPA.o : Ecc_SPA.c Ecc_SPA.h Ecc_Binary.h Constants.h Ecc_Math.h
# 	$(CC) $(CCFLAGS) -c Ecc_SPA.c

IO_ass2 : IO_ass2.c Constants.h
	$(CC) $(CCFLAGS) -c IO_ass2.c 
	
	
	
	
stud_t : $(OBJS) stud_t_sampler.c Constants.h
	$(CC) $(CCFLAGS) -o stud_t stud_t_sampler.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)
	
pet : $(OBJS) PE_stud_t.c Constants.h
	$(CC) $(CCFLAGS) -o pet PE_stud_t.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)
	
rj : $(OBJS) rj_norm_t.c Constants.h
	$(CC) $(CCFLAGS) -o rj rj_norm_t.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)
	
clean: 
	rm stud_t pet rj *.o  
