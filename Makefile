all:bam_to_mods


bam_to_mods: obj/bam_mods_to_text.o
	$(CC) -g -o $@ obj/bam_mods_to_text.o -lhts $(CFLAGS) $(LDFLAGS)  $(LIBS) -lpthread


obj/bam_mods_to_text.o: src/bam_mods_to_text.c
	$(CC) -c -g -o $@ src/bam_mods_to_text.c -lhts ${CFLAGS} $(LIBS) -lpthread