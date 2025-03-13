#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void convert_csv_to_header(const char *csv_filename, const char *header_filename, const char *array_name) {
    FILE *csv = fopen(csv_filename, "r");
    FILE *header = fopen(header_filename, "w");

    if (!csv || !header) {
        printf("Error opening file: %s or %s\n", csv_filename, header_filename);
        return;
    }

    fprintf(header, "#ifndef %s_H\n#define %s_H\n\n", array_name, array_name);
    fprintf(header, "#define %s_LENGTH /* fill this in after counting elements */\n\n", array_name);
    fprintf(header, "double %s[] = {\n", array_name);

    char line[1024];
    int count = 0;
    while (fgets(line, sizeof(line), csv)) {
        char *token = strtok(line, ",\n");
        while (token != NULL) {
            fprintf(header, "    %s,", token);
            count++;
            token = strtok(NULL, ",\n");
        }
    }

    fprintf(header, "\n};\n\n#define %s_LENGTH %d\n\n", array_name, count);
    fprintf(header, "#endif // %s_H\n", array_name);

    fclose(csv);
    fclose(header);

    printf("Converted %s → %s (Total elements: %d)\n", csv_filename, header_filename, count);
}

int main() {
    // List of CSV files and matching header/array names
    const char *csv_files[] = {
        "navic_prn_navic_prn.csv"
    };
    const char *array_names[] = {
        "navic_input"
    };

    for (int i = 0; i < 5; i++) {
        char header_filename[256];
        sprintf(header_filename, "%s.h", array_names[i]);
        convert_csv_to_header(csv_files[i], header_filename, array_names[i]);
    }

    return 0;
}
