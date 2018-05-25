//this program walks the root and lib/ directories of the project to generate a .ninja file
//yes this is a pun on "bob the builder" don't @ me

//bob/build.sh builds bob.
//you only need to run it once before you can build the rest of the project

//TODO: walk folder structures recursively
//TODO: generate build files for Windows
//TODO: generate .app bundles for macOS
//TODO: test and make sure this works with symlinks

//wishful thinking section
//TODO: build for Android
//TODO: build for iOS

#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>
#include <string.h>

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

enum ArgType {
    NONE, ROOT, TEMP, OUTPUT,
};

//poor man's ad hoc linked list
struct Node {
    Node * next;
    char * name;
};

Node * add_files(
        const char * root,
        const char * inDir,
        const char * outDir,
        const char * inExtension,
        const char * outExtension,
        const char * rule,
        Node * list,
        FILE * out)
{
    //open directory
    char * dir = (char *) malloc(strlen(root) + strlen(inDir) + 2);
    strcpy(dir, root);
    strcat(dir, "/");
    strcat(dir, inDir);
    DIR * dp = opendir(dir);
    if (dp == nullptr) {
        fprintf(stderr, "WARNING: Could not open directory %s, assuming empty\n", dir);
        return nullptr;
    }
    free(dir);

    //scan for source files
    dirent * ep = readdir(dp);
    while (ep != nullptr) {
        char * ext = strrchr(ep->d_name, '.');
        // if (ext != nullptr && (strcmp(ext + 1, "cpp") == 0 ||
        //                        strcmp(ext + 1, "cc" ) == 0 ||
        //                        strcmp(ext + 1, "c"  ) == 0)) {
        if (ext != nullptr && strcmp(ext + 1, inExtension) == 0) {
            //add file to head of linked list
            int len = ext - ep->d_name;
            Node * head = (Node *) malloc(sizeof(Node));
            head->next = list;
            head->name = (char *) malloc(len + 1);
            strncpy(head->name, ep->d_name, len);
            head->name[len] = '\0';
            list = head;

            //write build command to output file
            fprintf(out, "build %s/%s.%s: %s %s/%s\n",
                outDir, head->name, outExtension, rule, inDir, ep->d_name);
        }
        ep = readdir(dp);
    }

    closedir(dp);

    return list;
}

int main(int argc, char ** argv) {
    const char * root = ".";
    const char * temp = "template.ninja";
    const char * output = "build.ninja";

    //parse command line arguments
    ArgType type = NONE;
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-r") == 0) {
            type = ROOT;
        } else if (strcmp(argv[i], "-t") == 0) {
            type = TEMP;
        } else if (strcmp(argv[i], "-o") == 0) {
            type = OUTPUT;
        } else if (type == ROOT) {
            root = argv[i];
            type = NONE;
        } else if (type == TEMP) {
            temp = argv[i];
            type = NONE;
        } else if (type == OUTPUT) {
            output = argv[i];
            type = NONE;
        }
    }

    //read template ninja file
    char * part = read_entire_file(temp);
    if (part == nullptr) {
        printf("Error while reading file %s\n", temp);
        exit(-1);
    }

    //paste template into output file
    FILE * out = fopen(output, "w");
    fprintf(out, "%s\n", part);
    free(part);



    //gather asset files into a linked list
    Node * assetList = add_files(root, "assets", "res", "obj", "73", "obj", nullptr, out);
    //free immediately because we're not doing anything with these
    while (assetList != nullptr) {
        Node * next = assetList->next;
        free(assetList->name);
        free(assetList);
        assetList = next;
    }



    //gather cpp files into a linked list
    Node * node = add_files(root, "lib", "$builddir", "cpp", "o", "lib", nullptr, out);
    fprintf(out, "\n");
    node = add_files(root, ".", "$builddir", "cpp", "o", "cc", node, out);

    //write link command to output file and free linked list
    fprintf(out, "\nbuild game: link");
    while (node != nullptr) {
        fprintf(out, " $builddir/%s.o", node->name);
        Node * next = node->next;
        free(node->name);
        free(node);
        node = next;
    }
    fprintf(out, "\n");

    fclose(out);

    return 0;
}
