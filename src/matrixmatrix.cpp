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

	std::vector<double> tmp_m(m);
	std::vector<double> tmp_l(l);

	for(int i = 0; i < m; i++){
		matA.push_back(tmp_l);
	}
	for(int i = 0; i < n; i++){
		matB.push_back(tmp_m);
		matC.push_back(tmp_l);
	}

	for(int i = 0;i<l;i++)
	{
		matA[0][i] = i;
	}

	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&l, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int part_size = n / (world_size-1);
	int part_size_last = n % (world_size-1);
	std::vector<int> partitions(world_size);
	std::fill(partitions.begin(), partitions.end(),part_size);
	partitions[world_size-1] = part_size_last;
	std::vector<int> space(world_size);
	for(int i = 1; i<world_size;i++)
	{
		space[i] = space[i-1] + partitions[i-1];
	}
	std::vector<double> tmp(part_size);
	for(int i = 0; i<m;i++)
	{
		MPI_Scatterv(&(matA[i][0]),&(partitions[0]),&space[0],MPI_DOUBLE, &(tmp[0]),part_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&l, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int part_size = n / (world_size-1);
	if(rank == world_size -1)
		part_size = n % (world_size-1);
	std::vector<std::vector<double> > matA_part;
	std::vector<std::vector<double> > matB;
	std::vector<std::vector<double> > matC;
	std::vector<double> tmp_l_part(part_size);
	for(int i = 0; i < m; i++){
		matA_part.push_back(tmp_l_part);
	}
	if(rank == 1)
	{
		std::cout << matA_part.size() << std::endl;
		std::cout << matA_part[0].size() << std::endl;
	}

	for(int i = 0;i<m;i++)
	{
		MPI_Scatterv(NULL,NULL,NULL ,MPI_DOUBLE, &(matA_part[m][0]), part_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	if(rank == 1)
	{
		std::cout << matA_part.size() << std::endl;
		std::cout << matA_part[0].size() << std::endl;
		for(double & x: matA_part[0])
		{
			std::cout << "test" << std::endl;
			std::cout << x << std::endl;
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
