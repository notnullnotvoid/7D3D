#ifndef FILEBUFFER_HPP
#define FILEBUFFER_HPP

#include <stdlib.h>

//TODO: add method to re-align "head" pointer
class FileBuffer {
public:
    uint8_t * block;
    uint8_t * head;
    uint8_t * end;
public:
    FileBuffer(size_t bytes = 1024) {
        block = (uint8_t *) malloc(bytes);
        head = block;
        end = block + bytes;
    }

    void check(size_t bytes) {
        size_t byte = head - block;
        size_t size = end - block;

        //done in a loop so adding extremely large payloads will work
        while (byte + bytes >= size) {
            size *= 2;
        }

        block = (uint8_t *) realloc(block, size);
        head = block + byte;
        end = block + size;
    }

    template<typename TYPE>
    size_t write(TYPE t) {
        check(sizeof(t));
        *((TYPE *)head) = t;
        size_t pos = head - block;
        head += sizeof(t);
        return pos;
    }

    template<typename TYPE>
    size_t update(size_t pos, TYPE t) {
        uint8_t * ptr = block + pos;
        *((TYPE *)ptr) = t;
        return pos;
    }

    size_t size() {
        return head - block;
    }
};

#endif
