CC = gcc-7

LIBS  = gsl gslcblas m 
CCFLAGS += -g -Wall -std=gnu99 -fmax-errors=5 #-Werror 
#CCFLAGS += -g -ffast-math -Wall -O2 -ftree-vectorize -std=gnu99 -fmax-errors=5  

# OBJS = Ecc_SPA.o 

all : stud_t # $(OBJS)  

# Ecc_SPA.o : Ecc_SPA.c Ecc_SPA.h Ecc_Binary.h Constants.h Ecc_Math.h
# 	$(CC) $(CCFLAGS) -c Ecc_SPA.c

	
stud_t : $(OBJS) stud_t_sampler.c
	$(CC) $(CCFLAGS) -o stud_t stud_t_sampler.c $(OBJS) $(INCDIR:%=-I%) $(LIBDIR:%=-L%) $(LIBS:%=-l%)
	
clean: 
	rm stud_t # *.o  
