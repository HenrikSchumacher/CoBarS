#This is how I compile it on my Apple M1.

g++-14 -Wall -Wextra -Wno-ignored-qualifiers -mmacosx-version-min=13.0 -std=c++20 -m64 -Ofast -flto=auto -fno-math-errno -pthread -mcpu=apple-m1 -mtune=native -DNDEBUG -I.. main.cpp -oExample_ConfidenceSample_gcc
