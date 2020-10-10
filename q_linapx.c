/*
*********************************************************************
*	About: The code implements an algorithm to find Linear Approx.	*
*	Usage: Run in command prompt "q_linapx_repeat.exe <times>"		*
*	Replace <times> by a positive integer to repeat <times> times	*
*	Author: Ashwini Kumar Malviya									*
*	Email: ashwinixar@gmail.com										*
*********************************************************************
*/

//Find the best linear approximation for the following SBox F(x) : (F_2)^m -> (F_2)^n, where m = n = 4
// f(0x0) = 0xE, f(0x1) = 0x4, f(0x2) = 0xD, f(0x3) = 0x1,
// f(0x4) = 0x2, f(0x5) = 0xF, f(0x6) = 0xB, f(0x7) = 0x8,
// f(0x8) = 0x3, f(0x9) = 0xA, f(0xA) = 0x6, f(0xB) = 0xC,
// f(0xC) = 0x5, f(0xD) = 0x9, f(0xE) = 0x0, f(0xF) = 0x7
//The SBox f(x) can be represented by four component functions, f(x) = (F_1(x), F_2(x), F_3(x), F_4(x)), presented below in algebraic normal form:
// F_1(x) = x1x2x3 XOR x2x3x4 XOR x1x2 XOR x2x3 XOR x1 XOR x2 XOR x4 XOR 1
// F_2(x) = x1x2x3 XOR x1x3 XOR x2x4 XOR x3x4 XOR x1 XOR x2 XOR 1
// F_3(x) = x1x2x3 XOR x1x2x4 XOR x1x2 XOR x1x3 XOR x1x4 XOR x2x3 XOR x2x4 XOR x3x4 XOR x3 XOR x4 XOR 1
// F_4(x) = x1x3x4 XOR x1x4 XOR x2x4 XOR x1 XOR x3

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "QuEST.h"

int good_linapx(int _l_x, int _l_fx, int input_width, int *approximates)
{
	int fx[16] = { 0xE, 0x4, 0xD, 0x1, 0x2, 0xF, 0xB, 0x8, 0x3, 0xA, 0x6, 0xC, 0x5, 0x9, 0x0, 0x7 };
	int c1 = 0, c2 = 0;
	for(int i = 0; i < (int)pow(2, input_width); i++)
	{
		int r1 = 0, r2 = 0;
		int x = i;
		int l_x = _l_x, l_fx = _l_fx;
		while (x)
		{
			r1 ^= ((x & 1) & (l_x & 1));
			x >>= 1;
			l_x >>= 1;
		}
		while (fx[i])
		{
			r2 ^= ((fx[i] & 1) & (l_fx & 1));
			fx[i] >>= 1;
			l_fx >>= 1;
		}
		if(r1 == r2) c1++;
		else c2++;
	}
	*approximates = c1;
	return (c1 - c2);
}

int main (int narg, char *varg[])
{
	if(narg != 2)
	{
		printf("\nUsage: %s <times>\n", varg[0]);
		printf("Replace <times> by a positive integer to repeat <times> times\n");
		return 0;
	}

	int repeat = atoi(varg[1]);

	int m = 4; //Input "x" size (in bits) of SBox F(x)
	int n = 4; //Output "F(x)" size (in bits) of SBox F(x)

	QuESTEnv env = createQuESTEnv();

    Qureg qubits = createQureg(m + n, env);
    initZeroState(qubits);
    reportQuregParams(qubits);
    reportQuESTEnv(env);

    //Single qubit unitary for XOR gate
    ComplexMatrix2 ux = {
    	.real={{0,1},{1,0}},
		.imag={{0,0},{0,0}}
    };

	int count_trivial = 0, count_bad_linapx = 0;
	printf("-----------------------------------------\n");
	printf("Measured States for linear approximation:\n");
	printf("-----------------------------------------\n");
	printf("  l\tprob.\t  Approx. # of x|f(x)\n");
	printf("-----------------------------------------\n");
	while(repeat)
	{
		initZeroState(qubits);

		//Qubits indexed from 0 to (n - 1) represents f_n(x), ..., f_2(x), f_1(x)
		//Qubits indexed from n to (m + n - 1) represents xm, ..., x2, x1

		for(int i = n; i < m + n; i++)
			hadamard(qubits, i);

		//Evaluating f_1(x)
		int f1_ctrls[] = { 0xE, 0xC, 0x8, 0x7, 0x6, 0x4, 0x1, 0x0 };
		int _size = 8; //f1_ctrls size
		pauliX(qubits, m - 1); //f1_ctrl = 0x0 implies XORing with 1
		for(int i = 0; i < _size - 1; i++) //Loop for (_size - 1) because 0x0 is already applied
		{
			int *ctrls = (int *)malloc(sizeof(int) * m); //Maximum of "m" input bits can be control bits
			int ctrl_size = 0;
			int term = f1_ctrls[i];
			int qb = m;
			while(term)
			{
				if(term & 1) ctrls[ctrl_size++] = qb;
				term >>= 1;
				qb++;
			}
			multiControlledUnitary(qubits, ctrls, ctrl_size, m - 1, ux);
			free(ctrls);
		}

		//Evaluating f_2(x)
		int f2_ctrls[] = { 0xB, 0xA, 0x8, 0x5, 0x4, 0x3, 0x0 };
		_size = 7; //f2_ctrls size
		pauliX(qubits, 2); //f2_ctrl = 0x0 implies XORing with 1
		for(int i = 0; i < _size - 1; i++) //Loop for (_size - 1) because 0x0 is already applied
		{
			int *ctrls = (int *)malloc(sizeof(int) * m); //Maximum of "m" input bits can be control bits
			int ctrl_size = 0;
			int term = f2_ctrls[i];
			int qb = m;
			while(term)
			{
				if(term & 1) ctrls[ctrl_size++] = qb;
				term >>= 1;
				qb++;
			}
			multiControlledUnitary(qubits, ctrls, ctrl_size, 2, ux);
			free(ctrls);
		}

		//Evaluating f_3(x)
		int f3_ctrls[] = { 0xE, 0xD, 0xC, 0xA, 0x9, 0x6, 0x5, 0x3, 0x2, 0x1, 0x0 };
		_size = 11; //f3_ctrls size
		pauliX(qubits, 1); //f3_ctrl = 0x0 implies XORing with 1
		for(int i = 0; i < _size - 1; i++) //Loop for (_size - 1) because 0x0 is already applied
		{
			int *ctrls = (int *)malloc(sizeof(int) * m); //Maximum of "m" input bits can be control bits
			int ctrl_size = 0;
			int term = f3_ctrls[i];
			int qb = m;
			while(term)
			{
				if(term & 1) ctrls[ctrl_size++] = qb;
				term >>= 1;
				qb++;
			}
			multiControlledUnitary(qubits, ctrls, ctrl_size, 1, ux);
			free(ctrls);
		}

		//Evaluating f_4(x)
		int f4_ctrls[] = { 0xB, 0x9, 0x8, 0x5, 0x2 };
		_size = 5; //f4_ctrls size
		for(int i = 0; i < _size; i++)
		{
			int *ctrls = (int *)malloc(sizeof(int) * m); //Maximum of "m" input bits can be control bits
			int ctrl_size = 0;
			int term = f4_ctrls[i];
			int qb = m;
			while(term)
			{
				if(term & 1) ctrls[ctrl_size++] = qb;
				term >>= 1;
				qb++;
			}
			multiControlledUnitary(qubits, ctrls, ctrl_size, 0, ux);
			free(ctrls);
		}

		for(int i = 0; i < m + n; i++)
			hadamard(qubits, i);

		int l = 0;
		for(int i = 0; i < m + n; i++)
		{
			int outcome = measure(qubits, i);
			l ^= (outcome << i);
		}
		if(l == 0) count_trivial++;
		int x_mask = (l >> n);
		int fx_mask = (l & ((1 << n) - 1));
		int approximates;
		int temp_count = good_linapx(x_mask, fx_mask, m, &approximates);
		double prob = (temp_count / (pow(2, m) * pow(2, n / 2))) * (temp_count / (pow(2, m) * pow(2, n / 2)));
		if(temp_count == 0) count_bad_linapx++;
		printf("0x%02X\t%f\t%d\n", l, prob, approximates);

		repeat--;
	}

	printf("-----------------------------------------\n");
	printf("\nNumber of times trivial solution measured is %d", count_trivial);
	printf("\nNumber of times bad linear approximation measured is %d\n", count_bad_linapx);

	destroyQureg(qubits, env);
    destroyQuESTEnv(env);

    return 0;
}