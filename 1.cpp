#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <mpi.h>


struct Complex {
    double real;
    double imag;
};

double abs_complex(struct Complex z) {

    return sqrt(z.real * z.real + z.imag * z.imag);
}

struct Complex square_complex(struct Complex z) 
{
    struct Complex result;
    result.real = z.real * z.real - z.imag * z.imag;
    result.imag = 2 * z.real * z.imag;
    return result;
}

struct Complex add_complex(struct Complex z1, struct Complex z2) 
{
    struct Complex result;
    result.real = z1.real + z2.real;
    result.imag = z1.imag + z2.imag;
    return result;
}

int mandelbort_set(struct Complex c, int k)
{
    struct Complex z;
    int iter=0;
    z.real = 0;
    z.imag = 0;

    while (abs_complex(z) <= 2 && iter < k) 
    {
        z = add_complex(square_complex(z), c);
        iter++;
    }

    if (iter == k) return 1;
    else return 0;
       
    
}
int main(int argc, char *argv[]) 
{
    MPI_Init(&argc, &argv);
    int i, j, n, m, k, size, rank;
    struct Complex c;

    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    if (rank == 0) scanf("%d %d %d", &n, &m, &k);

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int flag = 0;
    int arr_size;
    if((n*m)%size != 0)
    {
        flag = 1;
        int remainder = (n*m) % size;
        int add_num = size - remainder;
        arr_size = (n*m) + add_num;
    }
    else arr_size = (n*m);

    std::complex<double> arr[arr_size];
    std::complex<double> recv_arr[arr_size/size];

    int sub_result_arr[arr_size/size];
    int result_arr[arr_size];

    if (rank == 0) 
    {
        int count=0;
        for (i = 0; i < m; i++) 
        {
            for (j = 0; j < n; j++) 
            {
                double x = -1.5 + (double)j * (2.5 / (double)(n-1));
                double y = 1.0 - (double)i * (2.0 / (double)(m-1));
                std::complex<double>  c_temp(x, y);
                arr[count] = c_temp;
                count++;
            }
        }
        if(flag)
        {
            for(i=count+1; i<arr_size; i++) 
            {
                std::complex<double> c_temp(10.0, 10.0);
                arr[i] = c_temp;
            }
        }
    }
    MPI_Scatter(arr, arr_size/size, MPI_C_DOUBLE_COMPLEX, recv_arr, arr_size/size,MPI_C_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    for(i=0; i<arr_size/size; i++)
    {
        c.real = real(recv_arr[i]);
        c.imag = imag(recv_arr[i]);
        sub_result_arr[i] = mandelbort_set(c, k);
    
    }

    MPI_Gather(sub_result_arr, arr_size/size, MPI_INT, result_arr, arr_size/size, MPI_INT, 0, MPI_COMM_WORLD);

    if(rank==0)
    {
        for(i=0; i<n*m; i++)
        {
            printf("%d ", result_arr[i]);
            if((i+1)%n==0) printf("\n");
        }
    }
    
    MPI_Finalize();
    return 0;
}