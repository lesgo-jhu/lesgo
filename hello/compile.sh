mkdir -p obj
g++ -c src/hello.cpp -o obj/hello.o
g++ -c src/main.cpp -o obj/main.o
g++ -o hello obj/main.o obj/hello.o
