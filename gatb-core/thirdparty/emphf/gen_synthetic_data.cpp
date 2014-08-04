#include <stdio.h>

int main(int argc, char** argv)
{
    if (argc != 3) {
        printf("Usage: %s <filename> <number>\n", argv[0]);
    }

    const char* filename = argv[1];
    unsigned long n = 0;
    sscanf(argv[2], "%lu", &n);

    printf("Generating %lu strings\n", n);

    auto f = fopen(filename, "w");
    for (unsigned long i = 0; i < n; ++i) {
        fprintf(f, "%lu\n", i);
    }
    fclose(f);
}
