rm game
bob/bob -t base.ninja -o build.ninja
export NINJA_STATUS='[%t/%f] %es '
ninja
rm build.ninja
./game
