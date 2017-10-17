#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include <omp.h>
#include <time.h>

int compare(const void* x, const void* y) {
    return (*(int*)x - *(int*)y);
}

void chipping(int* array, int n, int m) {
#pragma omp parallel for
    for (int i = 0; i < n; i += m) {
        qsort(array + i, (size_t)(i + m <= n ? m : n - i), sizeof(int), compare);
    }
}

void swap(int *a, int *b) {
    int temp = *b;
    *b = *a;
    *a = temp;
}

int binary_search(int* array, int key, int left, int right) {
    right++;
    while (left < right) {
        int mid = (left + right) / 2;
        if (key > array[mid]) {
            left = mid + 1;
        }
        else {
            right = mid;
        }
    }
    return right;
}

void single_flow_merge(int* array, int left_1, int right_1, int left_2, int right_2, int* result, int left_3) {
    int it_1 = 0;
    int it_2 = 0;

    while(left_1 + it_1 <= right_1 && left_2 + it_2 <= right_2) {
        if (array[left_1 + it_1] < array[left_2 + it_2]) {
            result[left_3 + it_1 + it_2] = array[left_1 + it_1];
            it_1++;
        } else {
            result[left_3 + it_1 + it_2] = array[left_2 + it_2];
            it_2++;
        }
    }

    while (left_1 + it_1 <= right_1) {
        result[left_3 + it_1 + it_2] = array[left_1 + it_1];
        it_1++;
    }

    while (left_2 + it_2 <= right_2) {
        result[left_3 + it_1 + it_2] = array[left_2 + it_2];
        it_2++;
    }
}

void merge(int* array, int left_1, int right_1, int left_2, int right_2, int* result, int left_3) {
    int size_1 = right_1 - left_1 + 1;
    int size_2 = right_2 - left_2 + 1;

    if (size_1 < size_2) {
        swap(&left_1, &left_2);
        swap(&right_1, &right_2);
        swap(&size_1, &size_2);
    }
    if (size_1 != 0) {
        int mid_1 = (right_1 + left_1) / 2;
        int mid_2 = binary_search(array, array[mid_1], left_2, right_2);
        int mid_3 = left_3 + (mid_1 - left_1) + (mid_2 - left_2);
        result[mid_3] = array[mid_1];
#pragma omp task
        {
            single_flow_merge(array, left_1, mid_1 - 1, left_2, mid_2 - 1, result, left_3);
        }
#pragma omp task
        {
            single_flow_merge(array, mid_1 + 1, right_1, mid_2, right_2, result, mid_3 + 1);
        }
    }
}

void merge_sort(int* array, int n, int m, int* result_array) {
    chipping(array, n, m);
    memcpy(result_array, array, sizeof(int) * n);

    for(int i = 2 * m; i <= 2 * n; i *= 2) {
#pragma omp parallel for shared(array, result_array)
        for (int j = 0; j < n; j += i) {
            int left_1 = j;
            int right_1 = j + i / 2 - 1;
            int left_2 = j + i / 2;
            int right_2 = j + i - 1;
            if (right_2 < n) {
                merge(array, left_1, right_1, left_2, right_2, result_array, left_1);
                continue;
            }
            if (left_2 < n) {
                merge(array, left_1, right_1, left_2, n - 1, result_array, left_1);
            }
        }
        memcpy(array, result_array, sizeof(int) * n);
    }
}

int main(int argc, char **argv) {
    srand(time(NULL));

    if (argc != 4) {
        printf("Not enough arguments!\n");
        return 0;
    }

    int n = atoi(argv[1]);
    int m = atoi(argv[2]);
    int p = atoi(argv[3]);

    int* array = (int*)malloc(sizeof(int) * n);
    if (array == NULL) {
        printf("Memory error");
    }

    for (int i = 0; i < n; i++) {
        array[i] = rand() % 100;
    }

    FILE* data = fopen("data.txt", "w");

    for (int i = 0; i < n; i++) {
        fprintf(data, "%d ", array[i]);
    }
    fprintf(data, "\n");

    omp_set_num_threads(p);

    int* array_for_qsort = (int*)malloc(sizeof(int) * n);
    if (array_for_qsort == NULL) {
        printf("Memory error");
    }
    memcpy(array_for_qsort, array, sizeof(int) * n);

    int* array_result = (int*)malloc(sizeof(int) * n);
    if (array_result == NULL) {
        printf("Memory error");
    }
    double start_time_merge = omp_get_wtime();
    merge_sort(array, n, m, array_result);
    double end_time_merge = omp_get_wtime();

    double start_time_quick = omp_get_wtime();
    qsort(array_for_qsort, n, sizeof(int), compare);
    double end_time_quick = omp_get_wtime();

    printf("Merge: %lf\n", end_time_merge - start_time_merge);
    printf("Qsort: %lf\n", end_time_quick - start_time_quick);

    for (int i = 0; i < n; i++) {
        fprintf(data, "%d ", array[i]);
    }

    FILE * stats = fopen("stats.txt", "w");
    fprintf(stats, "%lfs %d %d %d", end_time_merge - start_time_merge, n, m, p);

    fclose(data);
    fclose(stats);
    free(array);
    free(array_result);
    free(array_for_qsort);

    return 0;
}
