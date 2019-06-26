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

void main_task(int rank, int world_size)
{
	/*     m          n           n
	 *   _____      _____       _____
	 * 	|          |           |
	 * l|       x m|        = l|
	 *  |          |           |
	 *
	 */
	int m = 20;
	int l = 20;
	int n = 20;
	std::vector<std::vector<double> > matA;
	std::vector<std::vector<double> > matB;
	std::vector<std::vector<double> > matC;
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
	}


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
	std::cout << "start" << std::endl;
	for(int i = 0; i<m;i++)
	{
		MPI_Scatterv(&(matA[i][0]),&(partitions[0]),&space[0],MPI_DOUBLE, NULL,0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	for(int i = 1; i< world_size;i++)
	{
		std::cout << "test: " << i << std::endl;
		int start = i*part_size_n;
		int end = i!=world_size-1 ? (i+1) * part_size_n : i*part_size_n + part_size_n_last;
		for(int j = start;j<end;j++)
		{
			MPI_Send(&matB[j][0], m, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
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
	int part_size_n = 0;
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
	}
	else
	{
		part_size_n = n / (world_size-1);
		if(rank == world_size-1)
			part_size_n = n % (world_size-1);
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
	for(int i = 0;i< part_size_n;i++)
	{
		MPI_Status status;
		MPI_Recv(&tmp[0],m,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
		matB_part.push_back(tmp);
	}
	for(int i = 0;i< world_size;i++)
	{

		for(int j = 0;j<m;j++)
		{
			for(int k = 0;k<part_size_l;k++)
			{
				for(int l = 0;l<m;l++)
				{
					matC_part[j][k] += matA_part[l][k] * matB_part[j][l];
				}
			}
		}
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
