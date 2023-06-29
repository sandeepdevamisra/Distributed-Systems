#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <limits.h>
#include <math.h>
#include <mpi.h>
#include <algorithm>

int cost[360000];
int root[360000];
int parent[600];
int key[600];
int freq[600];

typedef struct node
{
    int key;
    int freq;
 
} node;

bool compare(node a, node b)
{
    return a.key < b.key;
}


void merge_arrays(node* arr1, node* arr2, int n1, int n2, node* arr3)
{
    int i=0, j=0, k=0;
    while (i<n1 && j<n2) 
    {
        if (arr1[i].key < arr2[j].key) arr3[k++] = arr1[i++];
        else arr3[k++] = arr2[j++];
    }

    while (i<n1) arr3[k++] = arr1[i++];
    while (j<n2) arr3[k++] = arr2[j++];
}

void merge_k_arrays(node* temp_arr[], int i, int j, int arr_size, int size, node* output)
{
    
    if (i == j) 
    {
        int k=0;
        while(k < arr_size/size)
        {
            output[k] = temp_arr[i][k];
            k++;
        }
        return;
    }
 
    if (j == i + 1) 
    {
        merge_arrays(temp_arr[i], temp_arr[j], arr_size/size, arr_size/size, output);
        return;
    }

    int n1, n2;
 
    n1 = (arr_size/size) * (((i + j) / 2) - i + 1);
    n2 = (arr_size/size) * (j - ((i + j) / 2));

    node* out1 = (node*)malloc(n1* sizeof(node));
    node* out2 = (node*)malloc(n2* sizeof(node));

    merge_k_arrays(temp_arr, i, (i + j) / 2, arr_size, size, out1);
    merge_k_arrays(temp_arr, (i + j) / 2 + 1, j, arr_size, size, out2);
    merge_arrays(out1, out2, n1, n2, output);
    
    free(out1);
    free(out2);
}

void merge(node* arr, int l, int m, int r)
{
    node* temp = (node*)malloc((r-l+1)* sizeof(node));
    int i=l, j=m+1, k=0;

    while (i<=m && j<=r) 
    {
        if(arr[i].key <= arr[j].key) temp[k++] = arr[i++];
        else temp[k++] = arr[j++];
    }

    while (i<=m) temp[k++] = arr[i++];
    while (j<=r) temp[k++] = arr[j++];
       
    k=0;
    for(i=l; i<=r; i++) arr[i] = temp[k++];
    
}

void mergeSort(node* arr, int l, int r)
{
    if (l < r) 
    {
        int m = l + (r - l) / 2;
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
        merge(arr, l, m, r);
    }
}

void find_parent(int i, int j, int curr_root, int n)
{
    if(i>j) return;
    int index = root[(n+1)*i + j];
    parent[index] = curr_root;
    find_parent(i, index-1, index, n);
    find_parent(index+1, j, index, n);
}

int main(int argc, char** argv) 
{
    int rank, size, n, arr_size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    node *arr, *recv_arr, *key_freq;
    

    MPI_Datatype node_type;
    MPI_Datatype type[2] = {MPI_INT, MPI_INT};
    int blocklen[2] = {1, 1};
    MPI_Aint offset[2];
    offset[0] = offsetof(node, key);
    offset[1] = offsetof(node, freq);
  
    MPI_Type_create_struct(2, blocklen, offset, type, &node_type);
    MPI_Type_commit(&node_type);

    if(rank==0)
    {
        scanf("%d", &n);
        int flag = 0;
        if(n%size != 0)
        {
            flag = 1;
            int remainder = n % size;
            int add_num = size - remainder;
            arr_size = n + add_num;
        }
        else arr_size = n;
        arr = (node*)malloc(arr_size * sizeof(node));
    
        int i;
        for(i=0; i<n; i++) scanf("%d %d", &arr[i].key, &arr[i].freq);
        
        
        if(flag)
        {
            while(i<arr_size)
            {
                arr[i].key = 1000;
                arr[i].freq = -1;
                i++;
            }
        }

    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&arr_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    recv_arr = (node*)malloc((arr_size/size) * sizeof(node));

    MPI_Scatter(arr, arr_size/size, node_type, recv_arr, arr_size/size, node_type, 0, MPI_COMM_WORLD);
    for(int i=0; i<arr_size/size; i++)
    {
        mergeSort(recv_arr, 0, (arr_size/size)-1);
    }
    MPI_Gather(recv_arr, arr_size/size, node_type, arr, arr_size/size, node_type, 0, MPI_COMM_WORLD);

    if(rank == 0)
    {
        
        node* temp_arr[size];
        int i;
        for (i = 0; i < size; i++) temp_arr[i] = (node*)malloc((arr_size/size) * sizeof(node));
        for(i = 0; i<arr_size; i++)
        {
            int x = i/(arr_size/size);
            int y = i%(arr_size/size);
            temp_arr[x][y] = arr[i];
        }

        merge_k_arrays(temp_arr, 0, size-1, arr_size, size, arr);

        for(i = 0; i<size; i++) free(temp_arr[i]);

        
        key[0] = -1;
        freq[0] = -1;
        
    
        for(i=0; i<n; i++)
        {
            key[i+1] = arr[i].key; 
            freq[i+1] = arr[i].freq;
        }
        
    
        int t=1;
        while(t<=2)
        {
            int i=1, j=t;
            while(j<n+1)
            {
                switch(t)
                {
                    case 1:
                        cost[(n+1)*i + j] = freq[i];
                        root[(n+1)*i + j] = key[i];
                        break;

                    case 2:
                        int f1 = freq[i];
                        int f2 = freq[j];
                        if(f1 + 2*f2 < f2 + 2*f1)
                        {
                            cost[(n+1)*i + j] = f1 + 2*f2;
                            root[(n+1)*i + j] = key[i];
                        }
                        else
                        {
                            cost[(n+1)*i + j] = f2 + 2*f1;
                            root[(n+1)*i + j] = key[j];
                        }
                        break;
                }

                i++, j++;
            }

            t++;
        }
        
    }


    MPI_Bcast(cost, ((n+1) * (n+1)), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(root, ((n+1) * (n+1)), MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(key, n+1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(freq, n+1, MPI_INT, 0, MPI_COMM_WORLD);
    
    int t = 3;
    while(t<n+1)
    {
        int temp_flag = 0;
        int i, j;
        int temp_size = n+2-t;
        if(temp_size%size != 0)
        {
            temp_flag = 1;
            int remainder = temp_size % size;
            int add_num = size - remainder;
            temp_size = temp_size + add_num;
        }
    
        int big_arr_cost[temp_size];
        int big_arr_root[temp_size];
        int small_arr_cost[temp_size/size];
        int small_arr_root[temp_size/size];

        if(rank==0)
        {
            big_arr_cost[0] = -1;
            big_arr_root[0] = -1;
            i=1, j=t;
            while(j<n+1)
            {
                int index = (n+1)*i + j;
                big_arr_cost[i] = index;
                big_arr_root[i] = index;
                i++, j++;
            }
            if(temp_flag)
            {
                while(i<temp_size)
                {
                    big_arr_cost[i] = -1;
                    big_arr_root[i] = -1;
                    i++;
                }
            }
        
        }

        MPI_Scatter(big_arr_cost, temp_size/size, MPI_INT, small_arr_cost, temp_size/size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Scatter(big_arr_root, temp_size/size, MPI_INT, small_arr_root, temp_size/size, MPI_INT, 0, MPI_COMM_WORLD);

        int c=0;
        while(c<temp_size/size)
        {   
       
            int index = small_arr_cost[c];

            if(index > -1)
            {
                int x = index/(n+1);
                int y = index%(n+1);
                
                int mini = INT_MAX;
                int freq_sum = 0;
                int f=x;
                while(f<=y)
                {
                    freq_sum += freq[f];
                    f++;
                }
                
                int root_index;
                int r=x;
                while(r<=y)
                {
                    int left=0, right=0;
                    if(r!=x) left = cost[(n+1)*x + (r-1)];
                    if(r!=y) right = cost[(n+1)*(r+1) + y];
            
                    if(left + right + freq_sum < mini)
                    {
                        mini = left + right + freq_sum;
                        root_index = r;
                    }
                    r++;
                }
                small_arr_cost[c] = mini;
                small_arr_root[c] = key[root_index];

            }

            c++;
        }

        MPI_Gather(small_arr_cost, temp_size/size, MPI_INT, big_arr_cost, temp_size/size, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Gather(small_arr_root, temp_size/size, MPI_INT, big_arr_root, temp_size/size, MPI_INT, 0, MPI_COMM_WORLD);

        if(rank == 0)
        {
            i=1, j=t;
            while(j<n+1)
            {
                cost[(n+1)*i + j] = big_arr_cost[i];
                root[(n+1)*i + j] = big_arr_root[i]; 
                i++, j++; 
            }
    
        }
    
        MPI_Bcast(cost, ((n+1) * (n+1)), MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(root, ((n+1) * (n+1)), MPI_INT, 0, MPI_COMM_WORLD);

        t++;
    }

    if(rank==0)
    {
        printf("%d\n", cost[2*n+1]);
        find_parent(1, n, 0, n);
        for(int i=1; i<=n; i++) printf("%d ", parent[i]);
    }

    MPI_Type_free(&node_type);
    free(recv_arr);
    if(rank == 0) free(arr);
    MPI_Finalize();
    return 0;
}

