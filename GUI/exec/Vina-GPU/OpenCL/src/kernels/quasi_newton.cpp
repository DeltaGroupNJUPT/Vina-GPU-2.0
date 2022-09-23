#define EL_TYPE_SIZE 11
#define AD_TYPE_SIZE 20
#define XS_TYPE_SIZE 17
#define SY_TYPE_SIZE 18

//inline double norm3(const double* a) {
//	return sqrt(pown(a[0], 2) + pown(a[1], 2) + pown(a[2], 2));
//}


inline void change_cl_init(change_cl* g, const __global float* ptr) {
	for (int i = 0; i < 3; i++)g->position[i] = ptr[i];
	for (int i = 0; i < 3; i++)g->orientation[i] = ptr[i + 3];
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++)g->lig_torsion[i] = ptr[i + 3 + 3];
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)g->flex_torsion[i] = ptr[i + 3 + 3 + MAX_NUM_OF_LIG_TORSION];
	g->lig_torsion_size = ptr[3 + 3 + MAX_NUM_OF_LIG_TORSION + MAX_NUM_OF_FLEX_TORSION];
}

inline void change_cl_init_with_change(change_cl* g_new, const change_cl* g_old) {
	for (int i = 0; i < 3; i++)g_new->position[i] = g_old->position[i];
	for (int i = 0; i < 3; i++)g_new->orientation[i] = g_old->orientation[i];
	for (int i = 0; i < MAX_NUM_OF_LIG_TORSION; i++)g_new->lig_torsion[i] = g_old->lig_torsion[i];
	for (int i = 0; i < MAX_NUM_OF_FLEX_TORSION; i++)g_new->flex_torsion[i] = g_old->flex_torsion[i];
	g_new->lig_torsion_size = g_old->lig_torsion_size;
}

inline void output_type_cl_init(output_type_cl* out, __constant float* ptr); // Function prototype in mutate_conf.cpp
inline void output_type_cl_init_with_output(output_type_cl* out_new, const output_type_cl* out_old);// Function prototype in mutate_conf.cpp

void print_ouput_type(output_type_cl* x, int torsion_size) {
	for (int i = 0; i < 3; i++)printf("\nx.position[%d] = %0.16f", i, x->position[i]);
	for (int i = 0; i < 4; i++)printf("\nx.orientation[%d] = %0.16f", i, x->orientation[i]);
	for (int i = 0; i < torsion_size; i++)printf("\n x.torsion[%d] = %0.16f", i, x->lig_torsion[i]);
	printf("\n x.torsion_size = %f", x->lig_torsion_size);
}

void print_change(change_cl* g, int torsion_size) {
	for (int i = 0; i < 3; i++)printf("\ng.position[%d] = %0.16f", i, g->position[i]);
	for (int i = 0; i < 3; i++)printf("\ng.orientation[%d] = %0.16f", i, g->orientation[i]);
	for (int i = 0; i < torsion_size; i++)printf("\ng.torsion[%d] = %0.16f", i, g->lig_torsion[i]);
	printf("\ng.torsion_size = %f", g->lig_torsion_size);
}

inline int num_atom_types(int atu) {
	switch (atu) {
	case 0: return EL_TYPE_SIZE;
	case 1: return AD_TYPE_SIZE;
	case 2: return XS_TYPE_SIZE;
	case 3: return SY_TYPE_SIZE;
	default: printf("Kernel1:num_atom_types() ERROR!"); return INFINITY;
	}
}

inline void elementwise_product(float* out, const float* a, const float* b) {
	out[0] = a[0] * b[0];
	out[1] = a[1] * b[1];
	out[2] = a[2] * b[2];
}

inline float elementwise_product_sum(const float* a, const float* b) {
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline float access_m_data(__constant float* m_data, int m_i, int m_j, int i, int j, int k) {
	return m_data[i + m_i * (j + m_j * k)];
}

inline bool not_max(float x) {
	return (x < 0.1 * INFINITY);// Problem: replace max_fl with INFINITY?
}

inline void curl_with_deriv(float* e, float* deriv, float v, const float epsilon_fl) {
	if (*e > 0 && not_max(v)) {
		float tmp = (v < epsilon_fl) ? 0 : (v / (v + *e));
		*e *= tmp;
		for (int i = 0; i < 3; i++)deriv[i] *= pown(tmp, 2);
	}
}

inline void curl_without_deriv(float* e, float v, const float epsilon_fl) {
	if (*e > 0 && not_max(v)) {
		float tmp = (v < epsilon_fl) ? 0 : (v / (v + *e));
		*e *= tmp;
	}
}

float g_evaluate(			__constant	grid_cl*	g,
					const				float*		m_coords,				// double[3]
					const				float		slope,				// double
					const				float		v,					// double
										float*		deriv,				// double[3]
					const				float		epsilon_fl
) {
	int m_i = g->m_i;
	int m_j = g->m_j;
	int m_k = g->m_k;
	if(m_i * m_j * m_k == 0)printf("\nkernel2: g_evaluate ERROR!#1");
	float tmp_vec[3] = { m_coords[0] - g->m_init[0],m_coords[1] - g->m_init[1] ,m_coords[2] - g->m_init[2] };
	float tmp_vec2[3] = { g->m_factor[0],g->m_factor[1] ,g->m_factor[2] };
	float s[3];
	elementwise_product(s, tmp_vec, tmp_vec2); // 

	float miss[3] = { 0,0,0 };
	int region[3];
	int a[3];
	int m_data_dims[3] = { m_i,m_j,m_k };
	for (int i = 0; i < 3; i++){
		if (s[i] < 0) {
			miss[i] = -s[i];
			region[i] = -1;
			a[i] = 0;
			s[i] = 0;
		}
		else if (s[i] >= g->m_dim_fl_minus_1[i]) {
			miss[i] = s[i] - g->m_dim_fl_minus_1[i];
			region[i] = 1;
			if (m_data_dims[i] < 2)printf("\nKernel2: g_evaluate ERROR!#2");
			a[i] = m_data_dims[i] - 2;
			s[i] = 1;
		}
		else {
			region[i] = 0;
			a[i] = (int)s[i];
			s[i] -= a[i];
		}
		if (s[i] < 0)printf("\nKernel2: g_evaluate ERROR!#3");
		if (s[i] > 1)printf("\nKernel2: g_evaluate ERROR!#4");
		if (a[i] < 0)printf("\nKernel2: g_evaluate ERROR!#5");
		if (a[i] + 1 >= m_data_dims[i])printf("\nKernel2: g_evaluate ERROR!#5");
	}

	float tmp_m_factor_inv[3] = { g->m_factor_inv[0],g->m_factor_inv[1],g->m_factor_inv[2] };
	const float penalty = slope * elementwise_product_sum(miss, tmp_m_factor_inv);
	if (penalty <= -epsilon_fl)printf("\nKernel2: g_evaluate ERROR!#6");

	const int x0 = a[0];
	const int y0 = a[1];
	const int z0 = a[2];
		 
	const int x1 = x0 + 1;
	const int y1 = y0 + 1;
	const int z1 = z0 + 1;

	//const float f000 = access_m_data(g->m_data, m_i, m_j, x0, y0, z0);
	//const float f100 = access_m_data(g->m_data, m_i, m_j, x1, y0, z0);
	//const float f010 = access_m_data(g->m_data, m_i, m_j, x0, y1, z0);
	//const float f110 = access_m_data(g->m_data, m_i, m_j, x1, y1, z0);
	//const float f001 = access_m_data(g->m_data, m_i, m_j, x0, y0, z1);
	//const float f101 = access_m_data(g->m_data, m_i, m_j, x1, y0, z1);
	//const float f011 = access_m_data(g->m_data, m_i, m_j, x0, y1, z1);
	//const float f111 = access_m_data(g->m_data, m_i, m_j, x1, y1, z1);

	int base = (x0 + m_i * (y0 + m_j * z0)) * 8;
	__constant float* base_ptr = &g->m_data[base];

	const float f000 = *base_ptr;
	const float f100 = *(base_ptr + 1);
	const float f010 = *(base_ptr + 2);
	const float f110 = *(base_ptr + 3);
	const float f001 = *(base_ptr + 4);
	const float f101 = *(base_ptr + 5);
	const float f011 = *(base_ptr + 6);
	const float f111 = *(base_ptr + 7);

	const float x = s[0];
	const float y = s[1];
	const float z = s[2];
		  
	const float mx = 1 - x;
	const float my = 1 - y;
	const float mz = 1 - z;

	float f =
		f000 * mx * my * mz +
		f100 * x  * my * mz +
		f010 * mx * y  * mz +
		f110 * x  * y  * mz +
		f001 * mx * my * z	+
		f101 * x  * my * z	+
		f011 * mx * y  * z	+
		f111 * x  * y  * z  ;

	if (deriv) { // valid pointer
		const float x_g =
			f000 * (-1) * my * mz +
			f100 *   1  * my * mz +
			f010 * (-1) * y  * mz +
			f110 *	 1  * y  * mz +
			f001 * (-1) * my * z  +
			f101 *   1  * my * z  +
			f011 * (-1) * y  * z  +
			f111 *   1  * y  * z  ;


		const float y_g =
			f000 * mx * (-1) * mz +
			f100 * x  * (-1) * mz +
			f010 * mx *   1  * mz +
			f110 * x  *   1  * mz +
			f001 * mx * (-1) * z  +
			f101 * x  * (-1) * z  +
			f011 * mx *   1  * z  +
			f111 * x  *   1  * z  ;


		const float z_g =
			f000 * mx * my * (-1) +
			f100 * x  * my * (-1) +
			f010 * mx * y  * (-1) +
			f110 * x  * y  * (-1) +
			f001 * mx * my *   1  +
			f101 * x  * my *   1  +
			f011 * mx * y  *   1  +
			f111 * x  * y  *   1  ;

		float gradient[3] = { x_g, y_g, z_g };

		curl_with_deriv(&f, gradient, v, epsilon_fl);

		float gradient_everywhere[3];

		for (int i = 0; i < 3; i++) {
			gradient_everywhere[i] = ((region[i] == 0) ? gradient[i] : 0);
			deriv[i] = g->m_factor[i] * gradient_everywhere[i] + slope * region[i];
		}
		return f + penalty;
	}	
	else {  // none valid pointer
		printf("\nKernel2: g_evaluate ERROR!#7");
		curl_without_deriv(&f, v, epsilon_fl);
		return f + penalty;
	}
}


float ig_eval_deriv(						output_type_cl*		x,
											change_cl*			g, 
						const				float				v,
								__constant	ig_cl*				ig_cl_gpu,
											m_cl*				m_cl_gpu,
						const				float				epsilon_fl
) {
	float e = 0;
	int nat = num_atom_types(ig_cl_gpu->atu);
	for (int i = 0; i < m_cl_gpu->m_num_movable_atoms; i++) {
		int t = m_cl_gpu->atoms[i].types[ig_cl_gpu->atu];
		if (t >= nat) {
			for (int j = 0; j < 3; j++)m_cl_gpu->minus_forces.coords[i][j] = 0;
			continue;
		}
		float deriv[3];

		e = e + g_evaluate(&ig_cl_gpu->grids[t], m_cl_gpu->m_coords.coords[i], ig_cl_gpu->slope, v, deriv, epsilon_fl);

		for (int j = 0; j < 3; j++) m_cl_gpu->minus_forces.coords[i][j] = deriv[j];
	}
	return e;
}

inline void quaternion_to_r3(const float* q, float* orientation_m) {
	// Omit assert(quaternion_is_normalized(q));
	const float a = q[0];
	const float b = q[1];
	const float c = q[2];
	const float d = q[3];

	const float aa = a * a;
	const float ab = a * b;
	const float ac = a * c;
	const float ad = a * d;
	const float bb = b * b;
	const float bc = b * c;
	const float bd = b * d;
	const float cc = c * c;
	const float cd = c * d;
	const float dd = d * d;

	//Omit assert(eq(aa + bb + cc + dd, 1));
	matrix tmp;
	mat_init(&tmp, 0); // matrix with fixed dimension 3(here we treate this as a regular matrix(not triangular matrix!))

	matrix_set_element(&tmp, 3, 0, 0,		(aa + bb - cc - dd)	);
	matrix_set_element(&tmp, 3, 0, 1, 2 *	(-ad + bc)			);
	matrix_set_element(&tmp, 3, 0, 2, 2 *	(ac + bd)			);
							 
	matrix_set_element(&tmp, 3, 1, 0, 2 *	(ad + bc)			);
	matrix_set_element(&tmp, 3, 1, 1,		(aa - bb + cc - dd)	);
	matrix_set_element(&tmp, 3, 1, 2, 2 *	(-ab + cd)			);
							 
	matrix_set_element(&tmp, 3, 2, 0, 2 *	(-ac + bd)			);
	matrix_set_element(&tmp, 3, 2, 1, 2 *	(ab + cd)			);
	matrix_set_element(&tmp, 3, 2, 2,		(aa - bb - cc + dd)	);

	for (int i = 0; i < 9; i++) orientation_m[i] = tmp.data[i];
}

inline void local_to_lab_direction(			float* out,
									const	float* local_direction,
									const	float* orientation_m
) {
	out[0] =	orientation_m[0] * local_direction[0] +
				orientation_m[3] * local_direction[1] +
				orientation_m[6] * local_direction[2];
	out[1] =	orientation_m[1] * local_direction[0] +
				orientation_m[4] * local_direction[1] +
				orientation_m[7] * local_direction[2];
	out[2] =	orientation_m[2] * local_direction[0] +
				orientation_m[5] * local_direction[1] +
				orientation_m[8] * local_direction[2];
}

inline void local_to_lab(						float*		out,
							const				float*		origin,
							const				float*		local_coords,
							const				float*		orientation_m
) {
	out[0] = origin[0] + (	orientation_m[0] * local_coords[0] +
							orientation_m[3] * local_coords[1] +
							orientation_m[6] * local_coords[2]
							);			 
	out[1] = origin[1] + (	orientation_m[1] * local_coords[0] +
							orientation_m[4] * local_coords[1] +
							orientation_m[7] * local_coords[2]
							);			 
	out[2] = origin[2] + (	orientation_m[2] * local_coords[0] +
							orientation_m[5] * local_coords[1] +
							orientation_m[8] * local_coords[2]
							);
}

inline void angle_to_quaternion2(				float*		out,
									const		float*		axis,
												float		angle
) {
	if (norm3(axis) - 1 >= 0.001)printf("\nkernel2: angle_to_quaternion() ERROR!"); // Replace assert(eq(axis.norm(), 1));
	normalize_angle(&angle);
	float c = cos(angle / 2);
	float s = sin(angle / 2);
	out[0] = c;
	out[1] = s * axis[0];
	out[2] = s * axis[1];
	out[3] = s * axis[2];
}

void set(	const				output_type_cl* x,
								rigid_cl*		lig_rigid_gpu,
								m_coords_cl*	m_coords_gpu,	
			const				atom_cl*		atoms,				
			const				int				m_num_movable_atoms,
			const				float			epsilon_fl
) {
	//************** (root) node.set_conf **************// (CHECKED)
	for (int i = 0; i < 3; i++) lig_rigid_gpu->origin[0][i] = x->position[i]; // set origin
	for (int i = 0; i < 4; i++) lig_rigid_gpu->orientation_q[0][i] = x->orientation[i]; // set orientation_q
	quaternion_to_r3(lig_rigid_gpu->orientation_q[0], lig_rigid_gpu->orientation_m[0]);// set orientation_m
	// set coords
	int begin = lig_rigid_gpu->atom_range[0][0];
	int end =	lig_rigid_gpu->atom_range[0][1];
	for (int i = begin; i < end; i++) {
		local_to_lab(m_coords_gpu->coords[i], lig_rigid_gpu->origin[0], &atoms[i].coords[0], lig_rigid_gpu->orientation_m[0]);
	}
	//************** end node.set_conf **************//

	//************** branches_set_conf **************//
	//update nodes in depth-first order
	for (int current = 1; current < lig_rigid_gpu->num_children + 1; current++) { // current starts from 1 (namely starts from first child node)
		int parent = lig_rigid_gpu->parent[current];
		float torsion = x->lig_torsion[current - 1]; // torsions are all related to child nodes
		local_to_lab(	lig_rigid_gpu->origin[current],
						lig_rigid_gpu->origin[parent],
						lig_rigid_gpu->relative_origin[current],
						lig_rigid_gpu->orientation_m[parent]
						); // set origin
		local_to_lab_direction(	lig_rigid_gpu->axis[current],
								lig_rigid_gpu->relative_axis[current],
								lig_rigid_gpu->orientation_m[parent]
								); // set axis
		float tmp[4];
		float parent_q[4] = {	lig_rigid_gpu->orientation_q[parent][0],
								lig_rigid_gpu->orientation_q[parent][1] ,
								lig_rigid_gpu->orientation_q[parent][2] ,
								lig_rigid_gpu->orientation_q[parent][3] };
		float current_axis[3] = {	lig_rigid_gpu->axis[current][0],
									lig_rigid_gpu->axis[current][1],
									lig_rigid_gpu->axis[current][2] };

		angle_to_quaternion2(tmp, current_axis, torsion);
		angle_to_quaternion_multi(tmp, parent_q);
		quaternion_normalize_approx(tmp, epsilon_fl);

		for (int i = 0; i < 4; i++) lig_rigid_gpu->orientation_q[current][i] = tmp[i]; // set orientation_q
		quaternion_to_r3(lig_rigid_gpu->orientation_q[current], lig_rigid_gpu->orientation_m[current]); // set orientation_m

		// set coords
		begin = lig_rigid_gpu->atom_range[current][0];
		end =	lig_rigid_gpu->atom_range[current][1];
		for (int i = begin; i < end; i++) {
			local_to_lab(m_coords_gpu->coords[i], lig_rigid_gpu->origin[current], &atoms[i].coords[0], lig_rigid_gpu->orientation_m[current]);
		}
	}
	//************** end branches_set_conf **************//
}

void p_eval_deriv(						float*		out,
										int			type_pair_index,
										float		r2,
							__constant	p_cl*		p_cl_gpu,
					const				float		epsilon_fl
) {
	const float cutoff_sqr = p_cl_gpu->m_cutoff_sqr;
	if(r2 > cutoff_sqr) printf("\nkernel2: p_eval_deriv() ERROR!");
	__constant p_m_data_cl* tmp = &p_cl_gpu->m_data[type_pair_index];
	float r2_factored = tmp->factor * r2;
	if (r2_factored + 1 >= SMOOTH_SIZE) printf("\nkernel2: p_eval_deriv() ERROR!");
	int i1 = (int)(r2_factored);
	int i2 = i1 + 1;
	if (i1 >= SMOOTH_SIZE || i1 < 0)printf("\n kernel2: p_eval_deriv() ERROR!");
	if (i2 >= SMOOTH_SIZE || i2 < 0)printf("\n : p_eval_deriv() ERROR!");
	float rem = r2_factored - i1;
	if (rem < -epsilon_fl)printf("\nkernel2: p_eval_deriv() ERROR!");
	if (rem >= 1 + epsilon_fl)printf("\nkernel2: p_eval_deriv() ERROR!");
	float p1[2] = { tmp->smooth[i1][0], tmp->smooth[i1][1] };
	float p2[2] = { tmp->smooth[i2][0], tmp->smooth[i2][1] };
	float e = p1[0] + rem * (p2[0] - p1[0]);
	float dor = p1[1] + rem * (p2[1] - p1[1]);
	out[0] = e;
	out[1] = dor;
}

//(CHECKED)
inline void curl(float* e, float* deriv, float v, const float epsilon_fl) {
	if (*e > 0 && not_max(v)) {
		float tmp = (v < epsilon_fl) ? 0 : (v / (v + *e));
		(*e) = tmp * (*e);
		for (int i = 0; i < 3; i++)deriv[i] = deriv[i] * (tmp * tmp);
	}
}

//(CHECKED)
float eval_interacting_pairs_deriv(		__constant	p_cl*			p_cl_gpu,
									const				float			v,
									const				lig_pairs_cl*   pairs,
									const			 	m_coords_cl*		m_coords,
														m_minus_forces* minus_forces,
									const				float			epsilon_fl
) {
	float e = 0;

	for (int i = 0; i < pairs->num_pairs; i++) {
		const int ip[3] = { pairs->type_pair_index[i], pairs->a[i] ,pairs->b[i] };
		float coords_b[3] = { m_coords->coords[ip[2]][0], m_coords->coords[ip[2]][1], m_coords->coords[ip[2]][2] };
		float coords_a[3] = { m_coords->coords[ip[1]][0], m_coords->coords[ip[1]][1], m_coords->coords[ip[1]][2] };
		float r[3] = { coords_b[0] - coords_a[0], coords_b[1] - coords_a[1] ,coords_b[2] - coords_a[2] };
		float r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
	
		if (r2 < p_cl_gpu->m_cutoff_sqr) {
			float tmp[2];
			p_eval_deriv(tmp, ip[0], r2, p_cl_gpu, epsilon_fl);
			float force[3] = { r[0] * tmp[1], r[1] * tmp[1] ,r[2] * tmp[1] };
			curl(&tmp[0], force, v, epsilon_fl);
			e += tmp[0];
			for (int j = 0; j < 3; j++)minus_forces->coords[ip[1]][j] -= force[j];
			for (int j = 0; j < 3; j++)minus_forces->coords[ip[2]][j] += force[j];
		}
	}
	return e;
}

inline void product(float* res, const float*a,const float*b) {
	res[0] = a[1] * b[2] - a[2] * b[1];
	res[1] = a[2] * b[0] - a[0] * b[2];
	res[2] = a[0] * b[1] - a[1] * b[0];
}

void POT_deriv(	const					m_minus_forces* minus_forces,
				const					rigid_cl*		lig_rigid_gpu,
				const					m_coords_cl*		m_coords,
										change_cl*		g
) {
	int num_torsion = lig_rigid_gpu->num_children;
	int num_rigid = num_torsion + 1;
	float position_derivative_tmp[MAX_NUM_OF_RIGID][3];
	float position_derivative[MAX_NUM_OF_RIGID][3];
	float orientation_derivative_tmp[MAX_NUM_OF_RIGID][3];
	float orientation_derivative[MAX_NUM_OF_RIGID][3];
	float torsion_derivative[MAX_NUM_OF_RIGID]; // torsion_derivative[0] has no meaning(root node has no torsion)

	for (int i = 0; i < num_rigid; i++) {
		int begin = lig_rigid_gpu->atom_range[i][0];
		int end = lig_rigid_gpu->atom_range[i][1];
		for (int k = 0; k < 3; k++)position_derivative_tmp[i][k] = 0; 
		for (int k = 0; k < 3; k++)orientation_derivative_tmp[i][k] = 0;
		for (int j = begin; j < end; j++) {
			for (int k = 0; k < 3; k++)position_derivative_tmp[i][k] += minus_forces->coords[j][k];

			float tmp1[3] = {	m_coords->coords[j][0] - lig_rigid_gpu->origin[i][0],
								m_coords->coords[j][1] - lig_rigid_gpu->origin[i][1],
								m_coords->coords[j][2] - lig_rigid_gpu->origin[i][2] };
			float tmp2[3] = {  minus_forces->coords[j][0],
								minus_forces->coords[j][1],
								minus_forces->coords[j][2] };
			float tmp3[3];
			product(tmp3, tmp1, tmp2);
			for (int k = 0; k < 3; k++)orientation_derivative_tmp[i][k] += tmp3[k];
		}
	}

	// position_derivative 
	for (int i = num_rigid - 1; i >= 0; i--) {// from bottom to top
		for (int k = 0; k < 3; k++)position_derivative[i][k] = position_derivative_tmp[i][k]; // self
		// looking for chidren node
		for (int j = 0; j < num_rigid; j++) {
			if (lig_rigid_gpu->children_map[i][j] == true) {
				for (int k = 0; k < 3; k++)position_derivative[i][k] += position_derivative[j][k]; // self+children node
			}
		}
	}

	// orientation_derivetive
	for (int i = num_rigid - 1; i >= 0; i--) { // from bottom to top
		for (int k = 0; k < 3; k++)orientation_derivative[i][k] = orientation_derivative_tmp[i][k]; // self
		// looking for chidren node
		for (int j = 0; j < num_rigid; j++) {
			if (lig_rigid_gpu->children_map[i][j] == true) { // self + children node + product
				for (int k = 0; k < 3; k++)orientation_derivative[i][k] += orientation_derivative[j][k];
				float product_out[3];
				float origin_temp[3] = {	lig_rigid_gpu->origin[j][0] - lig_rigid_gpu->origin[i][0],
											lig_rigid_gpu->origin[j][1] - lig_rigid_gpu->origin[i][1],
											lig_rigid_gpu->origin[j][2] - lig_rigid_gpu->origin[i][2] };
				product(product_out, origin_temp, position_derivative[j]);
				for (int k = 0; k < 3; k++)orientation_derivative[i][k] += product_out[k];
			}
		}
	}

	// torsion_derivative
	for (int i = num_rigid - 1; i >= 0; i--) { // from bottom to top
		float sum = 0;
		for (int j = 0; j < 3; j++) sum += orientation_derivative[i][j] * lig_rigid_gpu->axis[i][j];
		torsion_derivative[i] = sum;
	}

	for (int k = 0; k < 3; k++)	g->position[k] = position_derivative[0][k];
	for (int k = 0; k < 3; k++) g->orientation[k] = orientation_derivative[0][k];
	for (int k = 0; k < num_torsion; k++) g->lig_torsion[k] = torsion_derivative[k + 1];// no meaning for node 0
}


float m_eval_deriv(					output_type_cl*		c,
										change_cl*			g,
										m_cl*				m_cl_gpu,
							__constant	p_cl*				p_cl_gpu,
							__constant	ig_cl*				ig_cl_gpu,
					const	__global	float*				v,
					const				float				epsilon_fl
) {
	set(c, &m_cl_gpu->ligand.rigid, &m_cl_gpu->m_coords, m_cl_gpu->atoms, m_cl_gpu->m_num_movable_atoms, epsilon_fl);

	float e = ig_eval_deriv(	c,
								g, 
								v[1],				
								ig_cl_gpu,
								m_cl_gpu,
								epsilon_fl							
							);
	
	e += eval_interacting_pairs_deriv(	p_cl_gpu,
										v[0],
										&m_cl_gpu->ligand.pairs,
										&m_cl_gpu->m_coords,
										&m_cl_gpu->minus_forces,
										epsilon_fl
									);

	POT_deriv(&m_cl_gpu->minus_forces, &m_cl_gpu->ligand.rigid, &m_cl_gpu->m_coords, g);

	return e;
}


// Only support one ligand, no flex !
inline float find_change_index_read(const change_cl* g, int index) {
	if (index < 3)return g->position[index];
	index -= 3;
	if (index < 3)return g->orientation[index];
	index -= 3;
	if (index < g->lig_torsion_size)return g->lig_torsion[index]; //
	printf("\nKernel2:find_change_index_read() ERROR!"); // Shouldn't be here
	return -1;
}

inline void find_change_index_write(change_cl* g, int index, float data) {
	if (index < 3) { g->position[index] = data; return; }
	index -= 3;
	if (index < 3) { g->orientation[index] = data; return; }
	index -= 3;
	if (index < g->lig_torsion_size) { g->lig_torsion[index] = data; return; } //
	printf("\nKernel2:find_change_index_write() ERROR!"); // Shouldn't be here
}

void minus_mat_vec_product(	const		matrix*		h,
							const		change_cl*	in,
										change_cl*  out
) {
	int n = h->dim;
	for (int i = 0; i < n; i++) {
		float sum = 0;
		for (int j = 0; j < n; j++) {
			sum += h->data[index_permissive(h, i, j)] * find_change_index_read(in, j);
		}
		find_change_index_write(out, i, -sum);
	}
}


inline float scalar_product(	const	change_cl*			a,
								const	change_cl*			b,
								int							n
) {
	float tmp = 0;
	for (int i = 0; i < n; i++) {
		tmp += find_change_index_read(a, i) * find_change_index_read(b, i);
	}
	return tmp;
}


float line_search(					 	m_cl*				m_cl_gpu,
							__constant	p_cl*				p_cl_gpu,
							__constant	ig_cl*				ig_cl_gpu,
										int					n,
					const				output_type_cl*		x,
					const				change_cl*			g,
					const				float				f0,
					const				change_cl*			p,
										output_type_cl*		x_new,
										change_cl*			g_new,
										float*				f1,
					const				float				epsilon_fl,
					const	__global	float*				hunt_cap
) {
	const float c0 = 0.0001;
	const int max_trials = 10;
	const float multiplier = 0.5;
	float alpha = 1;

	const float pg = scalar_product(p, g, n);

	for (int trial = 0; trial < max_trials; trial++) {

		output_type_cl_init_with_output(x_new, x);

		output_type_cl_increment(x_new, p, alpha, epsilon_fl);

		*f1 =  m_eval_deriv(x_new,
							g_new,
							m_cl_gpu,
							p_cl_gpu,
							ig_cl_gpu,
							hunt_cap,
							epsilon_fl
							);

		if (*f1 - f0 < c0 * alpha * pg)
			break;
		alpha *= multiplier;
	}
	return alpha;
}


bool bfgs_update(			matrix*			h,
					const	change_cl*		p,
					const	change_cl*		y,
					const	float			alpha,
					const	float			epsilon_fl
) {

	const float yp = scalar_product(y, p, h->dim);
	
	if (alpha * yp < epsilon_fl) return false;
	change_cl minus_hy;
	change_cl_init_with_change(&minus_hy, y);
	minus_mat_vec_product(h, y, &minus_hy);
	const float yhy = -scalar_product(y, &minus_hy, h->dim);
	const float r = 1 / (alpha * yp);
	const int n = 6 + p->lig_torsion_size;

	for (int i = 0; i < n; i++) {
		for (int j = i; j < n; j++) {
			float tmp = alpha * r * (find_change_index_read(&minus_hy, i) * find_change_index_read(p, j)
									+ find_change_index_read(&minus_hy, j) * find_change_index_read(p, i)) +
									+alpha * alpha * (r * r * yhy + r) * find_change_index_read(p, i) * find_change_index_read(p, j);

			h->data[i + j * (j + 1) / 2] += tmp;
		}
	}

	return true;
}



void bfgs(					output_type_cl*			x,
								change_cl*			g,
								m_cl*				m_cl_gpu,
					__constant	p_cl*				p_cl_gpu,
					__constant	ig_cl*				ig_cl_gpu,
			const	__global	float*				hunt_cap,
			const				float				epsilon_fl,
			const				int					max_steps
) 
{
	int n = 3 + 3 + x->lig_torsion_size; // the dimensions of matirx

	matrix h;
	matrix_init(&h, n, 0);
	matrix_set_diagonal(&h, 1);

	change_cl g_new;
	change_cl_init_with_change(&g_new, g);

	output_type_cl x_new;
	output_type_cl_init_with_output(&x_new, x);
	 
	float f0 = m_eval_deriv(	x,
								g,
								m_cl_gpu,
								p_cl_gpu,
								ig_cl_gpu,
								hunt_cap,
								epsilon_fl
							);

	float f_orig = f0;
	// Init g_orig, x_orig
	change_cl g_orig;
	change_cl_init_with_change(&g_orig, g);
	output_type_cl x_orig;
	output_type_cl_init_with_output(&x_orig, x);
	// Init p
	change_cl p;
	change_cl_init_with_change(&p, g);

	float f_values[MAX_NUM_OF_BFGS_STEPS + 1];
	f_values[0] = f0;

	for (int step = 0; step < max_steps; step++) {

		minus_mat_vec_product(&h, g, &p);
		float f1 = 0;

		const float alpha = line_search(	m_cl_gpu,
											p_cl_gpu,
											ig_cl_gpu,
											n,
											x,
											g,
											f0,
											&p,
											&x_new,
											&g_new,
											&f1,
											epsilon_fl,
											hunt_cap
										);

		change_cl y;
		change_cl_init_with_change(&y, &g_new);
		// subtract_change
		for (int i = 0; i < n; i++) {
			float tmp = find_change_index_read(&y, i) - find_change_index_read(g, i);
			find_change_index_write(&y, i, tmp);
		}
		f_values[step + 1] = f1;
		f0 = f1;
		output_type_cl_init_with_output(x, &x_new);
		if (!(sqrt(scalar_product(g, g, n)) >= 1e-5))break;
		change_cl_init_with_change(g, &g_new);

		if (step == 0) {
			float yy = scalar_product(&y, &y, n);
			if (fabs(yy) > epsilon_fl) {
				matrix_set_diagonal(&h, alpha * scalar_product(&y, &p, n) / yy);
			}
		}

		bool h_updated = bfgs_update(&h, &p, &y, alpha, epsilon_fl);
	}

	if (!(f0 <= f_orig)) {
		f0 = f_orig;
		output_type_cl_init_with_output(x, &x_orig);
		change_cl_init_with_change(g, &g_orig);
	}

	// write output_type_cl energy
	x->e = f0;
}