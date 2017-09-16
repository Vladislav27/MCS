#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv) {
    double start_time = omp_get_wtime();
    if(argc != 7) {
        printf("Not enough arguments!\n");
        return 0;
    }
    int a =  atoi(argv[1]);
    int b =  atoi(argv[2]);
    int x =  atoi(argv[3]);
    int n =  atoi(argv[4]);
    double p =  atof(argv[5]);
    int threads_count =  atoi(argv[6]);

    omp_set_num_threads(threads_count);

    srand(time(NULL));
    int *seed = (int*)malloc(threads_count*sizeof(int));
    for (int i = 0; i < threads_count; i++)
        seed[i] = rand();

    double counter_b = 0;
    double quantity_steps = 0;
	int steps;

#pragma omp parallel for reduction(+: counter_b, quantity_steps) private(steps)
    for (int i = 0; i < n; i++) {
        steps = x;
        while (1) {
            if (steps == a){
                break;
            }
            if (steps == b){
                counter_b++;
                break;
            }
            if (rand_r(&seed[omp_get_thread_num()]) % 100 < p * 100) {
                steps++;
            } else {
                steps--;
            }
            quantity_steps++;
        }
    }
    double probability = counter_b / n;
    double average_steps = quantity_steps / n;

    double finish_time = omp_get_wtime();

    FILE *f = fopen("stats.txt", "w");
    if (f == NULL) {
        printf("Dont find the file stats.txt\n");
        return 0;
    } else {
        fprintf(f, "%.2f %.1f %.4fs %d %d %d %d %.2f %d", probability, average_steps, finish_time - start_time, a, b, x,
                n, p, threads_count);
        fclose(f);
    }
    free(seed);
    return 0;
}
