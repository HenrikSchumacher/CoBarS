#This is how I compile it on my Apple M1.

clang++ -Wall -Wextra -mmacos-version-min=13.0 -std=c++20 -Ofast -fno-math-errno -pthread -flto=auto -fenable-matrix -march=native -mtune=native -DNDEBUG -I.. main.cpp -oExample_ConfidenceSample_clang
