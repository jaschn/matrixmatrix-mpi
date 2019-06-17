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
	int n = 2000;
	std::vector<std::vector<double> > matA;
	std::vector<std::vector<double> > matB;
//	std::vector<std::vector<double> > matC_ref;
	std::vector<std::vector<double> > matC;

//	std::vector<double> tmp_filled(2000);
	std::vector<double> tmp(n);
//	int x = 0;
//	for(unsigned i = 0; i < tmp_filled.size();i++)
//	{
//		tmp_filled.at(i) = x;
//		x++;
//		if(i>10)
//			x = 0;
//	}
//	for(unsigned i = 0; i < tmp.size(); i++){
//		matA.push_back(tmp_filled);
//		matB.push_back(tmp_filled);
//		matC.push_back(tmp);
//	}
	for(unsigned i = 0; i < tmp.size(); i++){
		matA.push_back(tmp);
		matB.push_back(tmp);
		matC.push_back(tmp);
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int part_size = n / world_size;
	int part_size_last = part_size + (n % world_size);
	for(int i = 1; i<world_size;i++)
	{
		int size = part_size;
		if(i == world_size-1)
			size = part_size_last;
		for(int j = part_size*i; j < part_size*i+size;j++)
		{
			MPI_Send(&matA[part_size*i][0], n, MPI_INT, i, 0, MPI_COMM_WORLD);
		}
	}
}

void worker_task(int rank, int world_size)
{
	int n = 0;
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int part_size = n / world_size;
	if(rank == world_size -1)
		part_size += n % world_size;
	std::vector<std::vector<double> > matA_part(part_size);
	std::vector<std::vector<double> > matB;
	std::vector<std::vector<double> > matC;
	for(int i = 0;i<part_size;i++)
	{
		MPI_Recv(&matA_part[i][0],)
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
