/*
 * main.cpp
 *
 *  Created on: 13.05.2019
 *      Author: rene
 */
#include <vector>
#include <thread>
#include "mpi.h"
#include <iostream>
#include <cstdlib>

void main_task(int rank, int world_size)
{
	/*     m          n           n
	 *   _____      _____       _____
	 * 	|          |           |
	 * l|       x m|        = l|
	 *  |          |           |
	 *
	 */
	int m = 2;
	int l = 2;
	int n = 2;
	std::vector<std::vector<double> > matA;
	std::vector<std::vector<double> > matB;
	std::vector<std::vector<double> > matC;
	std::vector<std::vector<double> > matC_ref;
	int part_size_l;
	int part_size_l_last;
	int part_size_n;
	int part_size_n_last;

	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&l, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	std::vector<double> tmp_m(m);
	std::vector<double> tmp_l(l);

	for(int i = 0; i < m; i++){
		matA.push_back(tmp_l);
	}
	for(int i = 0; i < n; i++){
		matB.push_back(tmp_m);
		matC.push_back(tmp_l);
		matC_ref.push_back(tmp_l);
	}
	std::cout << "matA" << std::endl;
	for(int i = 0;i<l;i++)
	{
		for(int j = 0;j<m;j++)
		{
			matA[j][i] = rand() % 4;
			std::cout << matA[j][i] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << "matB" << std::endl;
	for(int i = 0;i<m;i++)
	{
		for(int j = 0;j<n;j++)
		{
			matB[j][i] = rand() % 4;
			std::cout << matB[j][i] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	if(l%world_size == 0)
	{
		part_size_l = l / world_size;
		part_size_l_last = part_size_l;
	}
	else
	{
		part_size_l = l / (world_size-1);
		part_size_l_last = l % (world_size-1);
	}
	if(n%world_size == 0)
	{
		part_size_n = n / world_size;
		part_size_n_last = part_size_n;
	}
	else
	{
		part_size_n = n / (world_size-1);
		part_size_n_last = n % (world_size-1);
	}

	std::vector<int> partitions(world_size);
	std::fill(partitions.begin(), partitions.end(),part_size_l);
	partitions[world_size-1] = part_size_l_last;
	std::vector<int> space(world_size);
	for(int i = 1; i<world_size;i++)
	{
		space[i] = space[i-1] + partitions[i-1];
	}
	std::vector<double> tmp(part_size_l);
	for(int i = 0; i<m;i++)
	{
		MPI_Scatterv(&(matA[i][0]),&partitions[0],&space[0],MPI_DOUBLE, NULL,0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	for(int i = 1; i< world_size;i++)
	{
		int start = i*part_size_n;
		int end = i!=world_size-1 ? (i+1) * part_size_n : i*part_size_n + part_size_n_last;
		MPI_Send(&start, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		MPI_Send(&end, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		for(int j = start;j<end;j++)
		{
			MPI_Send(&matB[j][0], m, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
	int b_start = 0;
	int b_end = part_size_n;
	for(int i = world_size;i>0;i--)
	{
		if(i!=world_size)
		{
			b_start = i*part_size_n;
			b_end = i!=world_size-1 ? (i+1) * part_size_n : i*part_size_n + part_size_n_last;
		}
		for(int j = b_start;j<b_end;j++)
		{
			for(int k = 0;k<part_size_l;k++)
			{
				for(int h = 0;h<m;h++)
				{
					matC[j][k] += matA[h][k] * matB[j][h];
				}
			}
		}
		MPI_Send(&b_start, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
		MPI_Send(&b_end, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
		for(int j = b_start;j<b_end;j++)
		{
			MPI_Send(&matB[j][0], m, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
		}
	}
	for(int i = 0; i<m;i++)
	{
		MPI_Gatherv(NULL, 0, MPI_INT, &matC[i][0], &partitions[0], &space[0],  MPI_INT, 0, MPI_COMM_WORLD);
	}

    for(unsigned i = 0; i < matA.size(); i++){
		for(unsigned y = 0; y < matB.size(); y++){
			for(unsigned z = 0; z < matB.size(); z++){
				matC_ref[i][y] += matA[i][z] * matB[z][y];
			}
		}
	}

    std::cout << "matc" << std::endl;
	for(int i = 0;i<l;i++)
	{
		for(int j = 0;j<n;j++)
		{
			std::cout << matC[j][i] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
    std::cout << "matc_ref" << std::endl;
	for(int i = 0;i<l;i++)
	{
		for(int j = 0;j<n;j++)
		{
			std::cout << matC_ref[j][i] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl << "equal" << std::endl;
    bool equal = true;
	for(int i = 0;i<n;i++)
	{
		for(int j = 0;j<l;j++)
		{
			if(matC[i][j] != matC_ref[i][j])
				equal = false;
		}
	}
	std::cout << equal << std::endl;
}

void worker_task(int rank, int world_size)
{
	/*     m          n           n
	 *   _____      _____       _____
	 * 	|          |           |
	 * l|       x m|        = l|
	 *  |          |           |
	 *
	 */
	int m = 0;
	int l = 0;
	int n = 0;
	int part_size_l = 0;
	int part_size_n;
	int part_size_n_last;
	std::vector<std::vector<double> > matA_part;
	std::vector<std::vector<double> > matB_part;
	std::vector<std::vector<double> > matC_part;

	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&l, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(l%world_size == 0)
	{
		part_size_l = l / world_size;
	}
	else
	{
		part_size_l = l / (world_size-1);
		if(rank == world_size -1)
			part_size_l = l % (world_size-1);
	}
	if(n%world_size == 0)
	{
		part_size_n = n / world_size;
		part_size_n_last = part_size_n;
	}
	else
	{
		part_size_n = n / (world_size-1);
		part_size_n_last = n % (world_size-1);
	}


	std::vector<double> tmp_l_part(part_size_l);
	matA_part.reserve(m);
	for(int i = 0; i < m; i++)
	{
		matA_part.push_back(tmp_l_part);
	}
	matA_part.reserve(n);
	for(int i = 0; i< n;i++)
	{
		matC_part.push_back(tmp_l_part);
	}

	for(int i = 0;i<m;i++)
	{
		MPI_Scatterv(NULL,NULL,NULL ,MPI_DOUBLE, &(matA_part[i][0]), matA_part[i].size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	matB_part.reserve(part_size_n);
	std::vector<double> tmp(m);
	for(int i = 0;i<part_size_n;i++)
	{
		matB_part.push_back(tmp);
	}
	int b_start;
	int b_end;
	for(int i = 0;i< world_size;i++)
	{
		std::cout << i << std::endl;
		MPI_Status status;
		MPI_Recv(&b_start, 1,MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&b_end , 1,MPI_INT, MPI_ANY_SOURCE,0, MPI_COMM_WORLD, &status);
		for(int i = 0;i< b_end-b_start;i++)
		{
			MPI_Recv(&tmp[0],m,MPI_DOUBLE, MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
			matB_part[i] = tmp;
		}
		for(int j = b_start;j<b_end;j++)
		{
			for(int k = 0;k<part_size_l;k++)
			{
				for(int h = 0;h<m;h++)
				{
					matC_part[j][k] += matA_part[h][k] * matB_part[j-b_start][h];
				}
			}
		}
		if(rank!=world_size-1)
		{
			MPI_Send(&b_start, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
			MPI_Send(&b_end, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
			for(int j = 0;j<b_end-b_start;j++)
			{
				MPI_Send(&matB_part[j][0], m, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
			}
		}
	}
    std::cout << "matc_part" << std::endl;
	for(int i = 0;i<part_size_l;i++)
	{
		for(int j = 0;j<n;j++)
		{
			std::cout << matC_part[j][i] << " ";
		}
		std::cout << std::endl;
	}
	for(int i = 0;i<n;i++)
	{
		MPI_Gatherv(&matC_part[i][0],part_size_l, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);
	}
}

int main(int argc, char **argv)
{
	int rank;
	int comm_size;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
	//setup
	if(rank == 0)
	{
		std::cout << "Comm Size:" << comm_size << std::endl;
		main_task(rank,comm_size);
	}
	else
	{
		worker_task(rank,comm_size);
	}
    MPI_Finalize();
	return 0;
}
