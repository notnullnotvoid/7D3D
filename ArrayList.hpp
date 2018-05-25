#ifndef ARRAYLIST_HPP
#define ARRAYLIST_HPP

#include <stdlib.h>

template <typename TYPE>
struct ArrayList {
    TYPE * data;
    size_t len;
    size_t max;

    void init(size_t reserve = 1024 / sizeof(TYPE) + 1) {
        data = (TYPE *) malloc(reserve * sizeof(TYPE));
        max = reserve;
        len = 0;
    }

    void add(TYPE t) {
        if (len == max) {
            max *= 2;
            data = (TYPE *) realloc(data, max * sizeof(TYPE));
        }

        data[len] = t;
        ++len;
    }

    void remove(size_t index) {
        --len;
        for (int i = index; i < len; ++i) {
            data[i] = data[i + 1];
        }
    }

    //removes elements in range [first, last)
    void remove(size_t first, size_t last) {
        size_t range = last - first;
        len -= range;
        for (int i = first; i < len; ++i) {
            data[i] = data[i + range];
        }
    }

    void shrink_to_fit() {
        data = (TYPE *) realloc(data, len * sizeof(TYPE));
    }

    void finalize() {
        free(data);
    }

    TYPE & operator[](size_t index) {
        return data[index];
    }
};

//convenience function for same-line declaration+initialization
template<typename TYPE>
ArrayList<TYPE> create_array_list(size_t reserve = 1024 / sizeof(TYPE) + 1) {
    ArrayList<TYPE> list;
    list.init(reserve);
    return list;
}

#endif
