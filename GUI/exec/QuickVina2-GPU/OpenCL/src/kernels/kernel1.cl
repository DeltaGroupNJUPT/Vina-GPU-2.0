#include "kernel2.h"

#define EL_TYPE_SIZE 11
#define AD_TYPE_SIZE 20
#define XS_TYPE_SIZE 17
#define SY_TYPE_SIZE 18

int fl_to_sz(float x, float max_sz) {
	if (x <= 0) return 0;
	if (x >= max_sz) return max_sz;
	return (int)x;
}
int num_atom_types(int atu) {
	switch (atu) {
	case 0: return EL_TYPE_SIZE;
	case 1: return AD_TYPE_SIZE;
	case 2: return XS_TYPE_SIZE;
	case 3: return SY_TYPE_SIZE;
	default: printf("\nkernel1:num_atom_types ERROR!"); return INFINITY;// replace assert()
	}
}

const float vec_distance_sqr(const __global float* a, const float* b) {
	return pown(a[0] - b[0], 2) + pown(a[1] - b[1], 2) + pown(a[2] - b[2], 2);
}

const int triangular_matrix_index(int n, int i, int j) {
	if (j >= n) printf("\nkernel1:triangular_matrix_index ERROR!");// replace assert()
	//assert(j < n);
	if (i > j) printf("\nkernel1:triangular_matrix_index ERROR!");// replace assert()
	return i + j * (j + 1) / 2;
}

const int triangular_matrix_index_permissive(int n, int i, int j) {
	return (i <= j) ? triangular_matrix_index(n, i, j) : triangular_matrix_index(n, j, i);
}

float eval_fast(int type_pair_index, float r2, float m_cutoff_sqr, __global p_m_data_cl* m_data, float factor) {
	if (r2 > m_cutoff_sqr) printf("\nkernel1:eval_fast ERROR!");// replace assert()
	if (r2 * factor >= FAST_SIZE)printf("\nkernel1:eval_fast ERROR!");// replace assert()
	int i = (int)(factor * r2);
	if (i >= FAST_SIZE)printf("\nkernel1:eval_fast ERROR!");// replace assert()
	float res = m_data[type_pair_index].fast[i];
	return res;
}

const __global int* __private possibilities(						float*		coords,
												const	__global	int*		ig_m_data,
																	float		epsilon_fl, 
												const	__global	float*		ig_m_init, 
												const	__global	float*		ig_m_range,
												const	__global	int*		m_data_dims
) {
	int index[3];
	float temp_array[3];
	for (int i = 0; i < 3; i++) {
		if (coords[i] + epsilon_fl < ig_m_init[i]) printf("\nkernel1:possibilities ERROR!1");//替换了assert()
		if (coords[i] > ig_m_init[i] + ig_m_range[i] + epsilon_fl) printf("\nkernel1:possibilities ERROR!2");//替换了assert()
		const float tmp = (coords[i] - ig_m_init[i]) * m_data_dims[i] / ig_m_range[i];
		temp_array[i] = tmp;
		index[i] = fl_to_sz(tmp, m_data_dims[i] - 1);//若0<tmp<m_data.dim(i) - 1则 强制整形并输出
	}
	int temp = index[0] + m_data_dims[0] * (index[1] + m_data_dims[1] * index[2]);
	int temp2 = (int)temp * MAX_NUM_OF_EVERY_M_DATA_ELEMENT;
	const __global int* __private address;
	address = &(ig_m_data[temp2]);
	return address;
}

__kernel
void \nkernel1(			__global	grid_cl*		g,
						__global	grid_atoms_cl*	ga,
						__global	p_cl*			p,
				const				int				atu,
				const				int				nat,
				const				float			epsilon_fl,
				const	__global	int*			ig_m_data,
				const	__global	float*			ig_m_init,
				const	__global	float*			ig_m_range,
				const	__global	int*			m_data_dims,
									int				needed_size,
				const	__global	int*			needed,
						__global	affinities_cl	affinities[]
) {	
	int x = get_global_id(0);
	int y = get_global_id(1);
	int z = get_global_id(2);
	if (x >= GRID_MI)return;
	if (y >= GRID_MJ)return;
	if (z >= GRID_MK)return;

	float probe_coords[3] = {  g->m_init[0] + g->m_factor_inv[0] * x,
								g->m_init[1] + g->m_factor_inv[1] * y,
								g->m_init[2] + g->m_factor_inv[2] * z };

	const __global int* __private possibilities_ptr;
	possibilities_ptr = possibilities(probe_coords, ig_m_data, epsilon_fl, ig_m_init, ig_m_range, m_data_dims);

	for (int possibilities_i = 0; possibilities_i < MAX_NUM_OF_EVERY_M_DATA_ELEMENT; possibilities_i++) {
		if (possibilities_ptr[possibilities_i] == 0) break;
		int i = possibilities_ptr[possibilities_i];
		const int t1 = ga->atoms[i].types[atu];
		if (t1 >= nat) continue;
		const float r2 = vec_distance_sqr(&ga->atoms[i].coords, probe_coords);
		if (r2 <= p->m_cutoff_sqr) {
			for (int j = 0; j < needed_size; j++) {
				const int t2 = needed[j];
				if (t2 > nat) printf("\nkernel1:t2 ERROR!");// replace assert()
				const int type_pair_index = triangular_matrix_index_permissive(num_atom_types(atu), t1, t2);
				affinities[z * (GRID_MJ * GRID_MI) + y * (GRID_MI) + x].data[j] += \
						eval_fast(type_pair_index, r2, p->m_cutoff_sqr, p->m_data, p->factor);
			}
		}
	}
}