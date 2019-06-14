/*
 * main.cpp
 *
 *  Created on: 13.05.2019
 *      Author: rene
 */
#include <vector>
#include <thread>
#include <omp.h>

void process_portion(std::vector<std::vector<double> > *matA, std::vector<std::vector<double> > *matB,
		std::vector<std::vector<double> > *matC, unsigned start, unsigned tpor){
	for(unsigned i = start; i < start+tpor; i++){
		for(unsigned y = 0; y < matB->size(); y++){
			for(unsigned z = 0; z < matB->size(); z++){
				(*matC)[i][y] += (*matA)[i][z] * (*matB)[z][y];
			}
		}
	}
}


int main(){

	std::vector<std::vector<double> > matA;
	std::vector<std::vector<double> > matB;
	std::vector<std::vector<double> > matC;

	std::vector<double> tmp(2000);

	for(unsigned i = 0; i < tmp.size(); i++){
		matA.push_back(tmp);
		matB.push_back(tmp);
		matC.push_back(tmp);
	}

	// compute tread portions
//    unsigned tpor = tmp.size() / 4;
//	std::thread first (process_portion, &matA, &matB, &matC, tpor*0, tpor);
//	std::thread second (process_portion, &matA, &matB, &matC, tpor*1, tpor);
//	std::thread third (process_portion, &matA, &matB, &matC, tpor*2, tpor);
//	std::thread fourth (process_portion, &matA, &matB, &matC, tpor*3, tpor);
//
//	first.join();
//	second.join();
//	third.join();
//	fourth.join();

	// openmp
omp_set_num_threads(4);
//
#pragma omp parallel for
   // openmp

    for(unsigned i = 0; i < matA.size(); i++){
		for(unsigned y = 0; y < matB.size(); y++){
			for(unsigned z = 0; z < matB.size(); z++){
				matC[i][y] += matA[i][z] * matB[z][y];
			}
		}
	}

	// clean up
	for(unsigned i = 0; i < tmp.size(); i++){
		matA[i].clear();
		matB[i].clear();
		matC[i].clear();

	}
	matA.clear();
	matB.clear();
	matC.clear();
    tmp.clear();

	return 0;
}

