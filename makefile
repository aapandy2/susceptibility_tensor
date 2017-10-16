susceptibility_tensor: susceptibility_tensor.c integrands.c integrator.c spline_calc.c 
	gcc -o susceptibility_tensor susceptibility_tensor.c integrands.c integrator.c spline_calc.c -lm -lgsl -lgslcblas -fopenmp
