rm ./obj
clang -std=c++1z -stdlib=libc++ -Wall -Ofast -o obj ../system.cpp obj.cpp -lc++
./obj -i ae86.obj -o ae86.73
