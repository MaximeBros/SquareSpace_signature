main:
	mkdir -p ./build
	gcc -g src/main.c src/flint_utils.c src/square_space.c src/sha2.c -lflint -lmpfr -lgmp -o build/main

clean:
	rm -rf ./build
