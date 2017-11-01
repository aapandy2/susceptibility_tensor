susceptibility_tensor: susceptibility_tensor.c integrands.c integrator.c  
	gcc -o susceptibility_tensor susceptibility_tensor.c integrands.c integrator.c -lm -lgsl -lgslcblas -fopenmp
