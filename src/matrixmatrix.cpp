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
#include "omp.h"
#include <chrono>

std::chrono::time_point<std::chrono::high_resolution_clock> start;
std::chrono::time_point<std::chrono::high_resolution_clock> finish;
std::chrono::milliseconds milliseconds;
void start_timer()
{

	start = std::chrono::high_resolution_clock::now();
}

void stop_timer(std::string name)
{
    finish = std::chrono::high_resolution_clock::now();
    milliseconds = std::chrono::duration_cast<std::chrono::milliseconds>(finish-start);
	std::cout << name << " took " <<  milliseconds.count() << " ms" << std::endl;
}

void check_with_reference(std::vector<std::vector<double>> &matA, std::vector<std::vector<double>> &matB,
						  std::vector<std::vector<double>> &matC, std::vector<std::vector<double>> &matC_ref, int world_size)
{
	std::cout << "Calculate reference with " << world_size << " threads." << std::endl;
	start_timer();
	omp_set_num_threads(world_size);
#pragma omp parallel for
    for(unsigned i = 0; i < matB.size(); i++){
		for(unsigned y = 0; y < matA[0].size(); y++){
			for(unsigned z = 0; z < matA.size(); z++){
				matC_ref[i][y] += matA[z][y] * matB[i][z];
			}
		}
	}
    stop_timer("omp");

    bool equal = true;
	for(unsigned i = 0;i<matC.size();i++)
	{
		for(unsigned j = 0;j<matC[0].size();j++)
		{
			if(matC[i][j] != matC_ref[i][j])
			{
				equal = false;
			}
		}
	}
	if(equal)
		std::cout << "both are equal" << std::endl;
	else
		std::cout << "not equal" << std::endl;
}

void get_dimensions(int &m, int &l, int &n, bool &check_reference)
{
	std::string check_ref;
	std::cout << "     m          n           n" << std::endl
			  << "   _____      _____       _____" << std::endl
			  << "  |          |           |" << std::endl
			  << " l|       x m|        = l|" << std::endl
			  << "  |          |           |" << std::endl << std::endl;
	std::cout << "m: ";
	std::cin >> m;
	std::cout << "l: ";
	std::cin >> l;
	std::cout << "n: ";
	std::cin >> n;
	std::cout << "Calculate Reference? [y/n]";
	std::cin >> check_ref;
	if(check_ref == "y" || check_ref == "Y")
	{
		check_reference = true;
	}
	else
	{
		check_reference = false;
	}
}

void main_task(int rank, int world_size)
{
	/*     m          n           n
	 *   _____      _____       _____
	 * 	|          |           |
	 * l|       x m|        = l|
	 *  |          |           |
	 *
	 */
	int m = 1000;
	int l = 1000;
	int n = 1000;
	bool check_ref;
	get_dimensions(m,l,n,check_ref);
	std::vector<std::vector<double> > matA;
	std::vector<std::vector<double> > matB;
	std::vector<std::vector<double> > matC;
	std::vector<std::vector<double> > matC_ref;
	int part_size_l;
	int part_size_l_last;
	int part_size_n;
	int part_size_n_last;
	int part_size_n_max;
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
	for(int i = 0;i<l;i++)
	{
		for(int j = 0;j<m;j++)
		{
			matA[j][i] = rand() % 100;
		}
	}
	for(int i = 0;i<m;i++)
	{
		for(int j = 0;j<n;j++)
		{
			matB[j][i] = rand() % 100;
		}
	}


	std::cout << "initialized" << std::endl;
	std::cout << "calculate matrix using MPI. World Size: " << world_size << std::endl;
	start_timer();
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
	part_size_n_max = std::max(part_size_n,part_size_n_last);

	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&part_size_l, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&part_size_l_last, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&part_size_n_max, 1, MPI_INT, 0, MPI_COMM_WORLD);

	std::vector<int> partitions(world_size);
	std::fill(partitions.begin(), partitions.end(),part_size_l);
	partitions[world_size-1] = part_size_l_last;
	std::vector<int> space(world_size);
	for(int i = 1; i<world_size;i++)
	{
		space[i] = space[i-1] + partitions[i-1];
	}

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
			MPI_Send(&matB[j][0], matB[j].size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
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
		if(world_size > 1)
		{
			MPI_Send(&b_start, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
			MPI_Send(&b_end, 1, MPI_INT, rank+1, 0, MPI_COMM_WORLD);
			for(int j = b_start;j<b_end;j++)
			{
				MPI_Send(&matB[j][0], matB[j].size(), MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
			}
		}
	}
	for(int i = 0; i<n;i++)
	{
		MPI_Gatherv(NULL, 0, MPI_DOUBLE, &matC[i][0], &partitions[0], &space[0],  MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	stop_timer("MPI");
	std::cout << "finished" << std::endl;
    MPI_Finalize();
    if(check_ref)
    	check_with_reference(matA,matB,matC,matC_ref,world_size);
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
	int n = 0;
	int part_size_l = 0;
	int part_size_n_max = 0;
	int b_start;
	int b_end;
	std::vector<std::vector<double> > matA_part;
	std::vector<std::vector<double> > matB_part;
	std::vector<std::vector<double> > matC_part;
	{
		int p_l;
		int p_l_last;
		MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&p_l, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&p_l_last, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&part_size_n_max, 1, MPI_INT, 0, MPI_COMM_WORLD);

		if(rank == world_size-1)
			part_size_l = p_l_last;
		else
			part_size_l = p_l;
	}
	std::vector<double> tmp_l_part(part_size_l);
	std::vector<double> tmp(m);

	matA_part.reserve(m);
	for(int i = 0; i < m; i++)
	{
		matA_part.push_back(tmp_l_part);
	}
	matC_part.reserve(n);
	for(int i = 0; i< n;i++)
	{
		matC_part.push_back(tmp_l_part);
	}
	matB_part.reserve(part_size_n_max);
	for(int i = 0;i<part_size_n_max;i++)
	{
		matB_part.push_back(tmp);
	}

	for(int i = 0;i<m;i++)
	{
		MPI_Scatterv(NULL,NULL,NULL ,MPI_DOUBLE, &(matA_part[i][0]), matA_part[i].size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	for(int i = 0;i< world_size;i++)
	{
		MPI_Status status;
		MPI_Recv(&b_start, 1,MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(&b_end , 1,MPI_INT, MPI_ANY_SOURCE,0, MPI_COMM_WORLD, &status);
		for(int i = 0;i< b_end-b_start;i++)
		{
			MPI_Recv(&tmp[0],tmp.size(),MPI_DOUBLE, MPI_ANY_SOURCE,0,MPI_COMM_WORLD,&status);
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
				MPI_Send(&matB_part[j][0], matB_part[j].size(), MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
			}
		}
	}
	for(int i = 0;i<n;i++)
	{
		MPI_Gatherv(&matC_part[i][0],part_size_l, MPI_DOUBLE, NULL, NULL, NULL, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	}
    MPI_Finalize();
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
	return 0;
}
