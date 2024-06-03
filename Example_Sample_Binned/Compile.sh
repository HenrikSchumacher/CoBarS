clang++ -Wall -Wextra -mmacos-version-min=13.0 -std=c++20 -Ofast -flto -fno-math-errno -pthread -fenable-matrix -march=native -mtune=native -DNDEBUG -I.. main.cpp -oExample_Sample_Binned
