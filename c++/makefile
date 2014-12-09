CC		= g++
DEBUG		= -g
OPTIMIZE	= -O
CFLAGS		= $(OPTIMIZE) -fpermissive -w
LDFLAGS		= $(OPTIMIZE)
BOOM		= BOOM
OBJ		= obj
GZ		= genezilla
LIBS		= -LBOOM -lBOOM -lgsl -lgslcblas
#---------------------------------------------------------
all: \
	BOOM \
	obj \
	compress-junctions \
	compress-pileup \
	dump-junctions \
	dump-pileup \
	extract-pileup-interval \
	genezilla

BOOM:
	ln -s ../BOOM

obj: 
	mkdir obj
.PHONY : clean
clean:
	@rm -f $(OBJ)/*.o
#---------------------------------------------------------
genezilla:
	cd genezilla ; make genezilla
#---------------------------------------------------------
$(OBJ)/compress-pileup.o:\
		compress-pileup.C
	$(CC) $(CFLAGS) -o $(OBJ)/compress-pileup.o -c \
		compress-pileup.C
#---------------------------------------------------------
compress-pileup: \
		$(OBJ)/compress-pileup.o
	$(CC) $(LDFLAGS) -o compress-pileup \
		$(OBJ)/compress-pileup.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/dump-pileup.o:\
		dump-pileup.C
	$(CC) $(CFLAGS) -o $(OBJ)/dump-pileup.o -c \
		dump-pileup.C
#---------------------------------------------------------
dump-pileup: \
		$(OBJ)/dump-pileup.o
	$(CC) $(LDFLAGS) -o dump-pileup \
		$(OBJ)/dump-pileup.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/extract-pileup-interval.o:\
		extract-pileup-interval.C
	$(CC) $(CFLAGS) -o $(OBJ)/extract-pileup-interval.o -c \
		extract-pileup-interval.C
#---------------------------------------------------------
extract-pileup-interval: \
		$(OBJ)/extract-pileup-interval.o \
		$(OBJ)/RnaJunction.o
	$(CC) $(LDFLAGS) -o extract-pileup-interval \
		$(OBJ)/extract-pileup-interval.o \
		$(OBJ)/RnaJunction.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/train-segmenter.o:\
		train-segmenter.C
	$(CC) $(CFLAGS) -o $(OBJ)/train-segmenter.o -c \
		train-segmenter.C
#---------------------------------------------------------
train-segmenter: \
		$(OBJ)/train-segmenter.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/RnaJunction.o
	$(CC) $(LDFLAGS) -o train-segmenter \
		$(OBJ)/train-segmenter.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/RnaJunction.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/RnaJunction.o:\
		$(GZ)/RnaJunction.H \
		$(GZ)/RnaJunction.C
	$(CC) $(CFLAGS) -o $(OBJ)/RnaJunction.o -c \
		$(GZ)/RnaJunction.C
#---------------------------------------------------------
$(OBJ)/RnaJunctions.o:\
		$(GZ)/RnaJunctions.H \
		$(GZ)/RnaJunctions.C
	$(CC) $(CFLAGS) -o $(OBJ)/RnaJunctions.o -c \
		$(GZ)/RnaJunctions.C
#---------------------------------------------------------
$(OBJ)/compress-junctions.o:\
		compress-junctions.C
	$(CC) $(CFLAGS) -o $(OBJ)/compress-junctions.o -c \
		compress-junctions.C
#---------------------------------------------------------
compress-junctions: \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/compress-junctions.o
	$(CC) $(LDFLAGS) -o compress-junctions \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/compress-junctions.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/get-intergenic-distances.o:\
		get-intergenic-distances.C
	$(CC) $(CFLAGS) -o $(OBJ)/get-intergenic-distances.o -c \
		get-intergenic-distances.C
#---------------------------------------------------------
get-intergenic-distances: \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/get-intergenic-distances.o
	$(CC) $(LDFLAGS) -o get-intergenic-distances \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/RnaJunctions.o \
		$(OBJ)/get-intergenic-distances.o \
		$(LIBS)
#---------------------------------------------------------
$(OBJ)/dump-junctions.o:\
		dump-junctions.C
	$(CC) $(CFLAGS) -o $(OBJ)/dump-junctions.o -c \
		dump-junctions.C
#---------------------------------------------------------
dump-junctions: \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/dump-junctions.o
	$(CC) $(LDFLAGS) -o dump-junctions \
		$(OBJ)/RnaJunction.o \
		$(OBJ)/dump-junctions.o \
		$(LIBS)
#---------------------------------------------------------

#--------------------------------------------------------
$(OBJ)/normalize-wig.o:\
		normalize-wig.C
	$(CC) $(CFLAGS) -o $(OBJ)/normalize-wig.o -c \
		normalize-wig.C
#---------------------------------------------------------
normalize-wig: \
		$(OBJ)/normalize-wig.o
	$(CC) $(LDFLAGS) -o normalize-wig \
		$(OBJ)/normalize-wig.o \
		$(LIBS)
#---------------------------------------------
#--------------------------------------------------------
$(OBJ)/evaluate.o:\
		evaluate.C
	$(CC) $(CFLAGS) -o $(OBJ)/evaluate.o -c \
		evaluate.C
#---------------------------------------------------------
evaluate: \
		$(OBJ)/evaluate.o
	$(CC) $(LDFLAGS) -o evaluate \
		$(OBJ)/evaluate.o \
		$(LIBS)
#---------------------------------------------


#--------------------------------------------------------
$(OBJ)/extend-segments-by-reads.o:\
		extend-segments-by-reads.C
	$(CC) $(CFLAGS) -o $(OBJ)/extend-segments-by-reads.o -c \
		extend-segments-by-reads.C
#---------------------------------------------------------
extend-segments-by-reads: \
		$(OBJ)/extend-segments-by-reads.o
	$(CC) $(LDFLAGS) -o extend-segments-by-reads \
		$(OBJ)/extend-segments-by-reads.o \
		$(LIBS)
#---------------------------------------------
