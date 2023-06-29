#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <mpi.h>

typedef struct point
{
    int x;
    int y;
    char dir;
} point;


int main(int argc, char** argv) 
{
    int rank, size, n, m, k, time, arr_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    point *arr, *recv_arr;

    MPI_Datatype point_type;
    MPI_Datatype type[3] = {MPI_INT, MPI_INT, MPI_CHAR};
    int blocklen[3] = {1, 1, 1};
    MPI_Aint offset[3];
    offset[0] = offsetof(point, x);
    offset[1] = offsetof(point, y);
    offset[2] = offsetof(point, dir);

    MPI_Type_create_struct(3, blocklen, offset, type, &point_type);
    MPI_Type_commit(&point_type);

    if(rank==0)
    {
        scanf("%d %d %d %d", &n, &m, &k, &time);
        int flag = 0;
        if(k%size != 0)
        {
            flag = 1;
            int remainder = k % size;
            int add_num = size - remainder;
            arr_size = k + add_num;
        }
        else arr_size = k;
        arr = (point*)malloc(arr_size * sizeof(point));
        int i;
        for(i=0; i<k; i++)
        {
            int temp_x, temp_y;
            char temp_dir; 
            scanf("%d %d %c", &temp_x, &temp_y, &temp_dir);
            arr[i].x = temp_x;
            arr[i].y = temp_y;
            arr[i].dir = temp_dir;
        }
        if(flag)
        {
            while(i<arr_size)
            {
                arr[i].x = -1;
                arr[i].y = -1;
                arr[i].dir = 'X';
                i++;
            }
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&time, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&arr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    recv_arr = (point*)malloc((arr_size/size) * sizeof(point));

    for (int t = 0; t < time; t++) 
    {
    
        if(rank==0)
        {
            point *temp_arr; 
            temp_arr = (point*)malloc(arr_size * sizeof(point)); 
            for(int i=0; i<arr_size; i++) temp_arr[i] = arr[i];

            for(int i=0; i<arr_size-1; i++)
            {
                int count=1;
                for(int j=i+1; j<arr_size; j++)
                {
                    if(temp_arr[i].dir != 'X' && temp_arr[j].dir != 'X' && temp_arr[i].x == temp_arr[j].x && temp_arr[i].y == temp_arr[j].y) count++;
                }

                if(count == 2)
                {
                
                    for(int j=i+1; j<arr_size; j++)
                    {
                        if(temp_arr[i].dir != 'X' && temp_arr[j].dir != 'X' && temp_arr[i].x == temp_arr[j].x && temp_arr[i].y == temp_arr[j].y)
                        {
                            if(temp_arr[i].dir == 'R' && temp_arr[j].dir == 'L' && temp_arr[i].y != 0 && temp_arr[i].y != m-1)
                            {
                                arr[j].y += 1;
                                arr[j].dir = 'u';
                                
                                arr[i].y -= 1;
                                arr[i].dir = 'd';

                                temp_arr[i].dir = 'X';
                                temp_arr[j].dir = 'X';
                                count=0;
                            }
                            else if(temp_arr[i].dir == 'L' && temp_arr[j].dir == 'R' && temp_arr[i].y != 0 && temp_arr[i].y != m-1)
                            {
                                arr[i].y += 1;
                                arr[i].dir = 'u';
                                
                                arr[j].y -= 1;
                                arr[j].dir = 'd';

                                temp_arr[i].dir = 'X';
                                temp_arr[j].dir = 'X';
                                count=0;
                            }
                            else if(temp_arr[i].dir == 'D' && temp_arr[j].dir == 'U' && temp_arr[i].x != 0 && temp_arr[i].x != n-1)
                            {
                                arr[j].x += 1;
                                arr[j].dir = 'r';

                                arr[i].x -= 1;
                                arr[i].dir = 'l';

                                temp_arr[i].dir = 'X';
                                temp_arr[j].dir = 'X';
                                count=0;
                            }
                            else if(temp_arr[i].dir == 'U' && temp_arr[j].dir == 'D' && temp_arr[i].x != 0 && temp_arr[i].x != n-1)
                            {
                                arr[i].x += 1;
                                arr[i].dir = 'r';

                                arr[j].x -= 1;
                                arr[j].dir = 'l';

                                temp_arr[i].dir = 'X';
                                temp_arr[j].dir = 'X';
                                count=0;
                            }
                    
                        }

                        if(count==0) break;
                    }

                }
            }

        }
        
   
        MPI_Scatter(arr, arr_size/size, point_type, recv_arr, arr_size/size, point_type, 0, MPI_COMM_WORLD);
    
        for(int i=0; i<arr_size/size; i++) 
        {
            if(recv_arr[i].dir == 'L')
            {
                if(recv_arr[i].x == 0)
                {
                    recv_arr[i].x += 1;
                    recv_arr[i].dir = 'R';
                }

                else recv_arr[i].x -= 1;
                
            }
            else if(recv_arr[i].dir == 'R')
            {
                if(recv_arr[i].x == n-1)
                {
                    recv_arr[i].x -= 1;
                    recv_arr[i].dir = 'L';
                }

                else recv_arr[i].x += 1;
            }
            else if(recv_arr[i].dir == 'U')
            {
                if(recv_arr[i].y == m-1)
                {
                    recv_arr[i].y -= 1;
                    recv_arr[i].dir = 'D';
                }

                else recv_arr[i].y += 1;
                
            }
            else if(recv_arr[i].dir == 'D')
            {
                if(recv_arr[i].y == 0)
                {
                    recv_arr[i].y += 1;
                    recv_arr[i].dir = 'U';
                }

                else recv_arr[i].y -= 1;
            }
            else if(recv_arr[i].dir == 'r') recv_arr[i].dir = 'R';
            else if(recv_arr[i].dir == 'l') recv_arr[i].dir = 'L';
            else if(recv_arr[i].dir == 'd') recv_arr[i].dir = 'D';
            else if(recv_arr[i].dir == 'u') recv_arr[i].dir = 'U';
        }
        
        

        MPI_Gather(recv_arr, arr_size/size, point_type, arr, arr_size/size, point_type, 0, MPI_COMM_WORLD);

    }
    

    if (rank == 0) 
    {
        for (int i = 0; i < k; i++)  
        {
            printf("%d %d %c\n", arr[i].x, arr[i].y, arr[i].dir);
        }
    }
    MPI_Type_free(&point_type);
    free(recv_arr);
    if(rank == 0) free(arr);
    MPI_Finalize();
    return 0;
}

