#include <iostream>
#include <mpi.h>

int main()
{
	MPI_Init(NULL, NULL);

	int world_size, world_rank, n, r0, r1;
	double sum = 0, min = 0, max = 0;
	int offset = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	if (world_rank == 0) {
		std::cout << "Podaj n: "; std::cin >> n;
		std::cout << "Podaj r0: "; std::cin >> r0;
		std::cout << "Podaj r1: "; std::cin >> r1;

		int* A = new int[n];

		srand(NULL);

		for (int i = 0; i < n; i++) {
			A[i] = rand() % (r1 - r0) + r0;
			std::cout << A[i] << " ";
		}

		std::cout << std::endl;

		delete[] A;

		int send_count = 0;
		for (int i = 1; i < world_size; i++) {
			
		}
	}

	if (world_rank == 0) {
		double srednia = sum / n;
		std::cout << "Srednia: " << srednia << std::endl;
		std::cout << "Minimum : " << min << std::endl;
		std::cout << "Maximum : " << max << std::endl;
	}

	//delete[] A;

	MPI_Finalize();

	/*
	MPI_Init(NULL, NULL);

	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
	
	std::cout << "Hello world from processor " << processor_name << " rank " <<
		world_rank << " out of " << world_size << " processors.\n";
	
	MPI_Finalize();*/

	return 0;
}