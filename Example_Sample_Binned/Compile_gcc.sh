# This is how I compile it on my Apple M1.
# It works with gcc 14 installed via homebrew. (brew install gcc@14.)

g++-14 -Wall -Wextra -Wno-ignored-qualifiers -std=c++20 -m64 -Ofast -flto=auto -fno-math-errno -pthread -mcpu=apple-m1 -mtune=native -DNDEBUG -I.. main.cpp -oExample_Sample_Binned_gcc

# On other systems you might want to replace -mcpu=apple-m1 by -march=native. (-march=native does not work on Apple Silicon due to a bug(?) in gcc.
