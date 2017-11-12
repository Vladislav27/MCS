#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <sys/time.h>
#include <mpi.h>

/*
 * left = 0
 * right = 1
 * up = 2
 * down = 4
 */

#define MAX(x, y) (((x) >= (y)) ? (x) : (y))

struct point {
    int x;
    int y;
};

struct point creater_point(int x, int y) {
    struct point arg;
    arg.x = x;
    arg.y = y;
    return arg;
}

struct run_args{
    int l;
    int a;
    int b;
    int n;
    int N;
    float p_l;
    float p_r;
    float p_u;
    float p_d;
    int rank;
    int size;
};

struct run_args creater_run_args(int l, int a, int b, int n, int N,
                                 float p_l, float p_r, float p_u, float p_d, int rank, int size) {
    struct run_args arg;
    arg.l = l;
    arg.a = a;
    arg.b = b;
    arg.n = n;
    arg.N = N;
    arg.p_l = p_l;
    arg.p_r = p_r;
    arg.p_u = p_u;
    arg.p_d = p_d;
    arg.rank = rank;
    arg.size = size;
    return arg;
}

struct area{
    struct point points;
    int rand;
    int n;
};

struct area creater_area(struct point points, int rand, int n){
    struct area arg;
    arg.points = points;
    arg.rand = rand;
    arg.n = n;
    return arg;
}

void insert(struct area** array, int* size, int* max_size, struct area* el) {
    if (*size < *max_size) {
        (*array)[*size] = *el;
        (*size)++;
    } else {
        *array = realloc(*array, *max_size * 2 * sizeof(struct area));
        *max_size *= 2;
        (*array)[*size] = *el;
        (*size)++;
    }
}

void delete(struct area** array, int* size, int i) {
    (*array)[i] = (*array)[(*size) - 1];
    (*size)--;
}

int max(float p0, float p1, float p2, float p3){
    float t1 = MAX(p0, p1);
    float t2 = MAX(p2, p3);
    float t3 = MAX(t1, t2);
    if (t3 == p0) {
        return 0;
    } else if (t3 == p1) {
        return 1;
    } else if (t3 == p2) {
        return 2;
    } else {
        return 3;
    }
}

int cheaker(int x, int y, int a, int b, int lrud) {
    if (lrud == 0) {
        if (x < 0) {
            x = a - 1;
        }
    }
    if (lrud == 1) {
        if (x >= a) {
            x = 0;
        }
    }
    if (lrud == 2) {
        if (y < 0) {
            y = b - 1;
        }
    }
    if (lrud == 3) {
        if (y >= b) {
            y = 0;
        }
    }
    return y * a + x;
}

int get_rank(struct point points, int a, int b, int lrud) {
    if (lrud == 0) {
        points.x -= 1;
    }
    if (lrud == 1) {
        points.x += 1;
    }
    if (lrud == 2) {
        points.y -= 1;
    }
    if (lrud == 3) {
        points.y += 1;
    }
    return cheaker(points.x, points.y, a, b, lrud);
}

void* run(void* args) {

    struct timeval start, finish;
    gettimeofday(&start, NULL);

    struct run_args* arg = args;
    int l = arg->l;
    int a = arg->a;
    int b = arg->b;
    int n = arg->n;
    int N = arg->N;
    int rank = arg->rank;
    int size = arg->size;
    float p_l = arg->p_l;
    float p_r = arg->p_r;
    float p_u = arg->p_u;
    float p_d = arg->p_d;
    struct point points = creater_point(rank % a, rank / a);

    int* ranks = (int*)malloc(sizeof(int) * 4);

    for(int i = 0; i < 4; i++) {
        ranks[i] = get_rank(points, a, b, i);
    }

    int* rands = malloc(sizeof(int) * size);
    if (rank == 0) {
        srand(time(NULL));
        for (int i = 0; i < size; i++)
            rands[i] = rand();
    }
    int seed;
    MPI_Scatter(rands, 1, MPI_INT, &seed, 1, MPI_INT, 0, MPI_COMM_WORLD);

    struct area* areas = (struct area*)malloc(sizeof(struct area) * N);
    srand(seed);
    for (int i = 0; i < N; i++) {
        areas[i] = creater_area(creater_point(rand() % l, rand() % l), rand(), n);
    }
    int areas_size = N;
    int areas_capacity = N;

    int* size_lrud = (int*)malloc(sizeof(int) * 5);
    int* real_size_lrud = (int*)malloc(sizeof(int) * 5);
    int* new_real_size_lrud = (int*)malloc(sizeof(int) * 4);
    for (int i = 0; i < 5; i++){
        size_lrud[i] = 0;
        real_size_lrud[i] = N;
    }

    for (int i = 0; i < 4; i++){
        new_real_size_lrud[i] = 0;
    }

    struct area* passing_particles0 = (struct area*)malloc(sizeof(struct area) * real_size_lrud[0]);
    struct area* passing_particles1 = (struct area*)malloc(sizeof(struct area) * real_size_lrud[1]);
    struct area* passing_particles2 = (struct area*)malloc(sizeof(struct area) * real_size_lrud[2]);
    struct area* passing_particles3 = (struct area*)malloc(sizeof(struct area) * real_size_lrud[3]);
    struct area* passing_particles_comp = (struct area*)malloc(sizeof(struct area) * real_size_lrud[4]);

    while (1) {
        int i = 0;
        while (i < areas_size) {
            struct area* areas_i = areas + i;
            int flag2 = 1;
            for (int t = 0; t < 190; t++) {
                if (areas_i->n == 0) {
                    insert(&passing_particles_comp, &size_lrud[4], &real_size_lrud[4], areas_i);
                    delete(&areas, &areas_size, i);
                    flag2 = 0;
                    break;
                }
                int direction = max(rand_r((unsigned int*) &areas_i->rand) * p_l, rand_r((unsigned int*) &areas_i->rand) * p_r,
                              rand_r((unsigned int*) &areas_i->rand) * p_u, rand_r((unsigned int*) &areas_i->rand) * p_d);
                if (direction == 0) {
                    areas_i->points.x -= 1;
                } else if (direction == 1) {
                    areas_i->points.x += 1;
                } else if (direction == 2) {
                    areas_i->points.y -= 1;
                } else {
                    areas_i->points.y += 1;
                }
                areas_i->n -= 1;
                if (areas_i->points.x < 0) {
                    areas_i->points.x = l - 1;
                    insert(&passing_particles0, &size_lrud[0], &real_size_lrud[0], areas_i);
                    delete(&areas, &areas_size, i);
                    flag2 = 0;
                    break;
                }
                if (areas_i->points.x >= l) {
                    areas_i->points.x = 0;
                    insert(&passing_particles1, &size_lrud[1], &real_size_lrud[1], areas_i);
                    delete(&areas, &areas_size, i);
                    flag2 = 0;
                    break;
                }
                if (areas_i->points.y < 0) {
                    areas_i->points.y = l - 1;
                    insert(&passing_particles2, &size_lrud[2], &real_size_lrud[2], areas_i);
                    delete(&areas, &areas_size, i);
                    flag2 = 0;
                    break;
                }
                if (areas_i->points.y >= l) {
                    areas_i->points.y = 0;
                    insert(&passing_particles3, &size_lrud[3], &real_size_lrud[3], areas_i);
                    delete(&areas, &areas_size, i);
                    flag2 = 0;
                    break;
                }
            }
            if (flag2) {
                i += 1;
            }
        }

        MPI_Request* data1 = (MPI_Request*) malloc(sizeof(MPI_Request) * 8);
        for(int r = 0; r < 4; r++){
            MPI_Isend(&size_lrud[r], 1, MPI_INT, ranks[r], r, MPI_COMM_WORLD, data1 + r);
        }
        MPI_Irecv(&new_real_size_lrud[0], 1, MPI_INT, ranks[0], 1, MPI_COMM_WORLD, data1 + 4);
        MPI_Irecv(&new_real_size_lrud[1], 1, MPI_INT, ranks[1], 0, MPI_COMM_WORLD, data1 + 5);
        MPI_Irecv(&new_real_size_lrud[2], 1, MPI_INT, ranks[2], 3, MPI_COMM_WORLD ,data1 + 6);
        MPI_Irecv(&new_real_size_lrud[3], 1, MPI_INT, ranks[3], 2, MPI_COMM_WORLD, data1 + 7);
        MPI_Waitall(8, data1, MPI_STATUS_IGNORE);

        struct area* receive0 = (struct area*) malloc(new_real_size_lrud[0] * sizeof(struct area));
        struct area* receive1 = (struct area*) malloc(new_real_size_lrud[1] * sizeof(struct area));
        struct area* receive2 = (struct area*) malloc(new_real_size_lrud[2] * sizeof(struct area));
        struct area* receive3 = (struct area*) malloc(new_real_size_lrud[3] * sizeof(struct area));

        MPI_Request* data = (MPI_Request*) malloc(sizeof(MPI_Request) * 8);
        MPI_Issend(passing_particles0, sizeof(struct area) * size_lrud[0], MPI_BYTE, ranks[0], 0, MPI_COMM_WORLD, data+0);
        MPI_Issend(passing_particles1, sizeof(struct area) * size_lrud[1], MPI_BYTE, ranks[1], 1, MPI_COMM_WORLD, data+1);
        MPI_Issend(passing_particles2, sizeof(struct area) * size_lrud[2], MPI_BYTE, ranks[2], 2, MPI_COMM_WORLD, data+2);
        MPI_Issend(passing_particles3, sizeof(struct area) * size_lrud[3], MPI_BYTE, ranks[3], 3, MPI_COMM_WORLD, data+3);
        MPI_Irecv(receive0, sizeof(struct area) * new_real_size_lrud[0], MPI_BYTE, ranks[0], 1, MPI_COMM_WORLD, data+4);
        MPI_Irecv(receive1, sizeof(struct area) * new_real_size_lrud[1], MPI_BYTE, ranks[1], 0, MPI_COMM_WORLD, data+5);
        MPI_Irecv(receive2, sizeof(struct area) * new_real_size_lrud[2], MPI_BYTE, ranks[2], 3, MPI_COMM_WORLD, data+6);
        MPI_Irecv(receive3, sizeof(struct area) * new_real_size_lrud[3], MPI_BYTE, ranks[3], 2, MPI_COMM_WORLD, data+7);
        MPI_Waitall(8, data, MPI_STATUS_IGNORE);

        for (int j = 0; j < new_real_size_lrud[0]; j++) {
            insert(&areas, &areas_size, &areas_capacity, receive0 + j);
        }
        for (int j = 0; j < new_real_size_lrud[1]; j++) {
            insert(&areas, &areas_size, &areas_capacity, receive1 + j);
        }
        for (int j = 0; j < new_real_size_lrud[2]; j++) {
            insert(&areas, &areas_size, &areas_capacity, receive2 + j);
        }
        for (int j = 0; j < new_real_size_lrud[3]; j++) {
            insert(&areas, &areas_size, &areas_capacity, receive3 + j);
        }
        for (int k = 0; k < 4; k++){
            size_lrud[k] = 0;
        }

        int* buff = (int*) malloc(sizeof(int));
        MPI_Reduce(&size_lrud[4], buff, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        int flag_tmp = 0;
        if (rank == 0) {
            if (*buff == size * N) {
                flag_tmp = 1;
            }
        }
        MPI_Bcast(&flag_tmp, 1, MPI_INT, 0, MPI_COMM_WORLD);

        free(buff);
        free(data1);
        free(data);
        free(receive0);
        free(receive1);
        free(receive2);
        free(receive3);

        if (flag_tmp) {
            int* result = (int*) malloc(sizeof(int) * size);
            MPI_Gather(&size_lrud[4], 1, MPI_INT, result, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (rank == 0) {
                gettimeofday(&finish, NULL);
                double delta = ((finish.tv_sec - start.tv_sec) * 1000000 + finish.tv_usec - start.tv_usec) / 1000000.0;
                FILE *file;
                file = fopen("stats.txt", "w");
                fprintf(file, "%d %d %d %d %d %f %f %f %f %fs\n", l, a, b, n, N, p_l, p_r, p_u, p_d, delta);
                for (int y = 0; y < size; y++) {
                    fprintf(file, "%d: %d\n", y, result[y]);
                }
                fclose(file);
            }
            free(result);
            break;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }

    free(passing_particles0);
    free(passing_particles1);
    free(passing_particles2);
    free(passing_particles3);
    free(rands);
    free(areas);
    free(passing_particles_comp);

    return NULL;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int l = atoi(argv[1]);
    int a = atoi(argv[2]);
    int b = atoi(argv[3]);
    int n = atoi(argv[4]);
    int N = atoi(argv[5]);
    float p_l = atof(argv[6]);
    float p_r = atof(argv[7]);
    float p_u = atof(argv[8]);
    float p_d = atof(argv[9]);

    struct run_args arg = creater_run_args(l, a, b, n, N, p_l, p_r, p_u, p_d, rank, size);

    pthread_t thread;
    pthread_create(&thread, NULL, run, &arg);
    pthread_join(thread, NULL);

    MPI_Finalize();
    return 0;
}
