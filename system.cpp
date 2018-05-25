#include "system.hpp"

#ifdef __APPLE__

#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>

//TODO: windows version of this
char * read_entire_file(const char * filepath) {
	int filedesc = open(filepath, O_RDONLY);
	if (filedesc == -1) {
		fprintf(stderr, "ERROR: Could not open file %s\n", filepath);
		return nullptr;
	}

	off_t size = lseek(filedesc, 0, SEEK_END);
	if (size == -1) {
		fprintf(stderr, "ERROR: Could not seek to end of file %s\n", filepath);
		return nullptr;
	}

	off_t seekErr = lseek(filedesc, 0, SEEK_SET);
	if (seekErr == -1) {
		fprintf(stderr, "ERROR: Could not seek to beginning of file %s\n", filepath);
		return nullptr;
	}

	char * text = (char *) malloc(size + 1);
	//TODO: error check malloc

	ssize_t bytesRead = read(filedesc, text, size);
	if (bytesRead == -1) {
		free(text);
		fprintf(stderr, "ERROR: Could not read file %s\n", filepath);
		return nullptr;
	} else if (bytesRead != size) {
		free(text);
		fprintf(stderr, "ERROR: Could only read part of file %s\n", filepath);
		return nullptr;
	}

	int closeErr = close(filedesc);
	if (closeErr == -1) {
		free(text);
		fprintf(stderr, "ERROR: Could not close file %s\n", filepath);
		return nullptr;
	}

	text[size] = '\0';

	return text;
}

#endif //__APPLE__
