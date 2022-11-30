qcc -O3 -Wall main.c -o a.out -lm
./a.out >out.ppm 2>log
rm a.out
