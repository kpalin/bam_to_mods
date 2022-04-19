all:bam_to_mods


bam_to_mods: obj/bam_mods_to_text.o
	$(CXX) -g -o $@ obj/bam_mods_to_text.o  -lhts $(CFLAGS) $(LDFLAGS)  $(LIBS) -lpthread


obj/bam_mods_to_text.o: src/bam_mods_to_text.cpp
	$(CXX) -c -g -o $@ $^  -lhts ${CFLAGS} $(LIBS) -lpthread 

clean: 
	rm obj/bam_mods_to_text.o bam_to_mods