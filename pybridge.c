#include <stdio.h>

FILE* py_popen(const char* cmd, const char* mode) {
    return popen(cmd, mode);
}

int py_pclose(FILE* f) {
    return pclose(f);
}

int py_write_int(FILE* f, int v) {
    return fprintf(f, "%d\n", v);
}

int py_write_double(FILE* f, double v) {
    return fprintf(f, "%.17e\n", v);
}

int py_flush(FILE* f) {
    return fflush(f);
}

int py_read_int(FILE* f, int* v) {
    return fscanf(f, "%d", v);
}

int py_read_double(FILE* f, double* v) {
    return fscanf(f, "%lf", v);
}
