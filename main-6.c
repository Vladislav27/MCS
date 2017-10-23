#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <time.h>
#include <pthread.h>
#include <omp.h>

struct qsort_args {
    int* array;
    int i;
    int n;
    int m;
};

struct single_flow_merge_args {
    int* array;
    int left_1;
    int right_1;
    int left_2;
    int right_2;
    int left_3;
    int* result;
};

struct merge_args {
    int* array;
    int* result;
    int n;
    pthread_t* threads;
    int i;
    int j;
    int p;
    int thread_id;
    struct single_flow_merge_args* single_flow_merge_args;
};

int compare(const void* x, const void* y) {
    return (*(int*)x - *(int*)y);
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

void* single_flow_chipping(void* args) {
    struct qsort_args* arg = args;
    qsort(arg->array + arg->i, arg->i + arg->m <= arg->n ? (size_t)arg->m : (size_t)(arg->n - arg->i),
          sizeof(int), compare);
}

void chipping(int* array, int n, int m, pthread_t* threads, int p) {
    int counter_threads = 0;

    struct qsort_args* args = (struct qsort_args*)malloc(sizeof(struct qsort_args) * p);
    for (int i = 0; i < n; ) {
        for (int j = 0; j < p; j++) {
            if (i >= n) {
                break;
            }

            struct qsort_args arg;
            arg.array = array;
            arg.i = i;
            arg.n = n;
            arg.m = m;
            args[j] = arg;
            pthread_create(threads + j, NULL, single_flow_chipping, args + j);
            counter_threads++;
            i += m;
        }
        for (int k = 0; k < counter_threads; k++) {
            pthread_join(threads[k], NULL);
        }
        counter_threads = 0;
    }
    free(args);
}

void* single_flow_merge(void* args) {
    struct single_flow_merge_args* arg = args;

    int* array = arg->array;
    int left_1 = arg->left_1;
    int right_1 = arg->right_1;
    int left_2 = arg->left_2;
    int right_2 = arg->right_2;
    int* result = arg->result;
    int left_3 = arg->left_3;

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

void merge(int* array, int left_1, int right_1, int left_2, int right_2, int* result, int left_3,
           pthread_t* threads, int p, int thread_id, struct single_flow_merge_args* single_flow_merge_arg) {
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

        struct single_flow_merge_args arg1;
        struct single_flow_merge_args arg2;

        arg1.left_1 = left_1;
        arg1.right_1 = mid_1 - 1;
        arg1.left_2 = left_2;
        arg1.right_2 = mid_2 - 1;
        arg1.left_3 = left_3;
        arg1.array = array;
        arg1.result = result;

        arg2.left_1 = mid_1 + 1;
        arg2.right_1 = right_1;
        arg2.left_2 = mid_2;
        arg2.right_2 = right_2;
        arg2.left_3 = mid_3 + 1;
        arg2.array = array;
        arg2.result = result;

        if (p > 1) {
            single_flow_merge_arg[thread_id] = arg1;
            single_flow_merge_arg[thread_id + p / 2] = arg2;
            pthread_create(threads + thread_id, NULL, single_flow_merge, single_flow_merge_arg + thread_id);
            single_flow_merge(single_flow_merge_arg + thread_id + p / 2);
            pthread_join(threads[thread_id], NULL);
        } else {
            single_flow_merge(&arg1);
            single_flow_merge(&arg2);
        }
    }
}

void* single_merge(void* args) {
    struct merge_args* arg = args;
    int left_1 = arg->j;
    int right_1 = arg->j + arg->i / 2 - 1;
    int left_2 = arg->j + arg->i / 2;
    int right_2 = arg->j + arg->i - 1;
    if (right_2 < arg->n) {
        merge(arg->array, left_1, right_1, left_2, right_2, arg->result, left_1, arg->threads,
              arg->p, arg->thread_id, arg->single_flow_merge_args);
    } else {
        if (left_2 < arg->n) {
            merge(arg->array, left_1, right_1, left_2, arg->n - 1, arg->result, left_1, arg->threads,
                  arg->p, arg->thread_id, arg->single_flow_merge_args);
        }
    }
}

void merge_sort(int* array, int n, int m, int* result_array, pthread_t* threads, int p) {
    chipping(array, n, m, threads, p);
    memcpy(result_array, array, sizeof(int) * n);

    if (p > 1) {
        int counter_threads = 0;

        struct merge_args* merge_arg = (struct merge_args*)malloc(sizeof(struct merge_args) * p / 2);
        struct single_flow_merge_args* single_flow_merge_arg = (struct single_flow_merge_args*)malloc(sizeof(struct single_flow_merge_args) * p);

        for (int i = 2 * m; i <= 2 * n; i *= 2) {
            for (int j = 0; j < n;) {
                for (int k = 0; k < p / 2; k++) {
                    if (j >= n) {
                        break;
                    }
                    struct merge_args arg;
                    arg.array = array;
                    arg.result = result_array;
                    arg.n = n;
                    arg.threads = threads;
                    arg.p = p;
                    arg.thread_id = k;
                    arg.i = i;
                    arg.j = j;
                    arg.single_flow_merge_args = single_flow_merge_arg;
                    merge_arg[k] = arg;
                    pthread_create(threads + k, NULL, single_merge, merge_arg + k);
                    counter_threads++;
                    j += i;
                }
                for (int l = 0; l < counter_threads; l++) {
                    pthread_join(threads[l], NULL);
                }
                counter_threads = 0;
            }
            memcpy(array, result_array, sizeof(int) * n);
        }
        free(merge_arg);
        free(single_flow_merge_arg);
    } else {
        struct merge_args* merge_arg = (struct merge_args*)malloc(sizeof(struct merge_args));
        for (int i = 2 * m; i <= 2 * n; i *= 2) {
            for (int j = 0; j < n; j += i) {
                struct merge_args arg;
                arg.array = array;
                arg.result = result_array;
                arg.n = n;
                arg.threads = threads;
                arg.p = p;
                arg.thread_id = 0;
                arg.i = i;
                arg.j = j;
                arg.single_flow_merge_args = NULL;
                merge_arg[0] = arg;
                single_merge(merge_arg);
            }
            memcpy(array, result_array, sizeof(int) * n);
        }
        free(merge_arg);
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

    pthread_t* threads = (pthread_t*)malloc((sizeof(pthread_t) * p));

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
    merge_sort(array, n, m, array_result, threads, p);
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
    free(threads);

    return 0;
}