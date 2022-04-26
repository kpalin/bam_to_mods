all:bam_to_mods


bam_to_mods: obj/bam_mods_to_text.o
	$(CXX) -O0  -g -o $@ obj/bam_mods_to_text.o  -lhts $(LDFLAGS)  $(LIBS) -lpthread


obj/bam_mods_to_text.o: src/bam_mods_to_text.cpp
	$(CXX)  -c  -o $@ $^  -lhts ${CFLAGS} $(LIBS) -lpthread -O0  -fno-tree-vectorize -Wall -g
	#$(CXX) -O0  -c -g -o $@ $^  -lhts ${CFLAGS} $(LIBS) -lpthread 

clean: 
	rm obj/bam_mods_to_text.o bam_to_mods