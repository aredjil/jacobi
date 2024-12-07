CC=mpic++
FLAGS= -I../include -DNON_BLOCKING -O3 -march=native -Wall -Wextra 
NP=1
main.x: main.cpp
	$(CC) $(FLAGS) main.cpp -o $@
run: 
	mpirun --mca -n $(NP) ./main.x 12000 10 
clean:
	rm -rf *.x
flush:
	rm -rf ../data/*.dat
	rm -rf ../data*.csv
	rm -rf ../data/*.txt
	rm -rf ../gif/*.gif
