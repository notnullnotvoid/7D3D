#undef DEBUG_PRINT

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "../types.hpp"
#include "../system.hpp"
#include "../FileBuffer.hpp"
#include "../ArrayList.hpp"

struct Vec3 {
	float x, y, z;
};

struct VertPair {
	u16 p, n; //position index, normal index
};

struct Triangle {
	VertPair pairs[3];
	u16 mat; //material index
};

struct Model {
	size_t vertexCount;
	Vec3 * vertices;
	void * verts; //allocated within the game

	size_t normalCount;
	Vec3 * normals;
	void * norms; //allocated within the game

	size_t triangleCount;
	Triangle * triangles;
};

enum ArgType {
	ARG_NONE, ARG_INFILE, ARG_OUTFILE,
};

int main(int argc, char ** argv) {
	const char * infile = "";
	const char * outfile = "";

	//parse command line arguments
	ArgType type = ARG_NONE;
	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-i") == 0) {
			type = ARG_INFILE;
		} else if (strcmp(argv[i], "-o") == 0) {
			type = ARG_OUTFILE;
		} else if (type == ARG_INFILE) {
			infile = argv[i];
			type = ARG_NONE;
		} else if (type == ARG_OUTFILE) {
			outfile = argv[i];
			type = ARG_NONE;
		}
	}



	char * file = read_entire_file(infile);

	//split file into lines
	ArrayList<char *> lines = create_array_list<char *>();
	char * line = strtok(file, "\r\n");
	while (line != nullptr) {
		lines.add(line);
		line = strtok(nullptr, "\r\n");
	}

	ArrayList<Vec3> positions = create_array_list<Vec3>();
	ArrayList<Vec3> normals = create_array_list<Vec3>();
	ArrayList<Triangle> faces = create_array_list<Triangle>();
	ArrayList<const char *> materialNames = create_array_list<const char *>();
	materialNames.add(""); //slot 0 = default (unassigned) material

	//parse lines
	u16 matIndex = 0;
	for (int i = 0; i < lines.len; ++i) {
		char * token = strtok(lines[i], " \t");

		if (strcmp(token, "v") == 0) {
			float x = atof(strtok(nullptr, " \t"));
			float y = atof(strtok(nullptr, " \t"));
			float z = atof(strtok(nullptr, " \t"));
			positions.add({ x, y, z });
		} else if (strcmp(token, "vn") == 0) {
			float x = atof(strtok(nullptr, " \t"));
			float y = atof(strtok(nullptr, " \t"));
			float z = atof(strtok(nullptr, " \t"));
			normals.add({ x, y, z });
		} else if (strcmp(token, "f") == 0) {
			Triangle tri = {};
			for (int j = 0; j < 3; ++j) {
				char * vert = strtok(nullptr, " \t");
				//assume the format uses two forward slashes
				char * slash = strchr(vert, '/');
				*slash = '\0';
				tri.pairs[j] = { (u16) (atoi(vert) - 1), (u16) (atoi(slash + 2) - 1) };
			}
			tri.mat = matIndex;
			faces.add(tri);
		} else if (strcmp(token, "usemtl") == 0) {
			char * name = strtok(nullptr, " \t");
			//find material name in list, or append if it doesn't exist
			bool found = false;
			for (int j = 0; j < materialNames.len; ++j) {
				if (strcmp(name, materialNames[j]) == 0) {
					found = true;
					matIndex = j;
					break;
				}
			}
			if (!found) {
				matIndex = materialNames.len;
				materialNames.add(name);
			}
		}

		//if the line doesn't begin with a known identifier, do nothing
		//this way, comments and properties we don't care about are ignored
	}



	FileBuffer data;

	//write model
	Model model = { positions.len, nullptr, nullptr,
					normals.len, nullptr, nullptr,
					faces.len, nullptr };
	size_t modelPos = data.write(model);

	//write positions
	data.update(modelPos + offsetof(Model, vertices), data.size());
	for (int i = 0; i < positions.len; ++i) {
		data.write(positions[i]);
	}

	//write normals
	data.update(modelPos + offsetof(Model, normals), data.size());
	for (int i = 0; i < normals.len; ++i) {
		data.write(normals[i]);
	}

	//write faces
	data.update(modelPos + offsetof(Model, triangles), data.size());
	for (int i = 0; i < faces.len; ++i) {
		data.write(faces[i]);
	}

	FILE * fp = fopen(outfile, "wb");
	fwrite(data.block, data.head - data.block, 1, fp);
	fclose(fp);



	#ifdef DEBUG_PRINT
		printf("\n\n\npositions:\n");
		for (int i = 0; i < positions.len; ++i) {
			printf("    %f, %f, %f\n", positions[i].x, positions[i].y, positions[i].z);
		}

		printf("\n\n\nnormals:\n");
		for (int i = 0; i < normals.len; ++i) {
			printf("    %f, %f, %f\n", normals[i].x, normals[i].y, normals[i].z);
		}

		printf("\n\n\nfaces:\n");
		for (int i = 0; i < faces.len; ++i) {
			printf("    %d//%d %d//%d %d//%d : %d : %s\n",
				faces[i].pairs[0].p, faces[i].pairs[0].n,
				faces[i].pairs[1].p, faces[i].pairs[1].n,
				faces[i].pairs[2].p, faces[i].pairs[2].n,
				faces[i].mat, materialNames[faces[i].mat]);
		}
	#endif



	faces.finalize();
	normals.finalize();
	positions.finalize();
	lines.finalize();
	free(file);

	return 0;
}
