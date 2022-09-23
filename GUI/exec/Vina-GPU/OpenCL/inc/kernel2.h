#pragma once
// Macros below are shared in both device and host
#define TOLERANCE 1e-16
// kernel1 macros
#define MAX_NUM_OF_EVERY_M_DATA_ELEMENT 512
#define MAX_M_DATA_MI 16
#define MAX_M_DATA_MJ 16
#define MAX_M_DATA_MK 16
#define MAX_NUM_OF_TOTAL_M_DATA MAX_M_DATA_MI*MAX_M_DATA_MJ*MAX_M_DATA_MK*MAX_NUM_OF_EVERY_M_DATA_ELEMENT

//kernel2 macros
#define MAX_NUM_OF_LIG_TORSION 48
#define MAX_NUM_OF_FLEX_TORSION 1
#define MAX_NUM_OF_RIGID 48
#define MAX_NUM_OF_ATOMS 130 
#define SIZE_OF_MOLEC_STRUC ((3+4+MAX_NUM_OF_LIG_TORSION+MAX_NUM_OF_FLEX_TORSION+ 1)*sizeof(float) )
#define SIZE_OF_CHANGE_STRUC ((3+3+MAX_NUM_OF_LIG_TORSION+MAX_NUM_OF_FLEX_TORSION + 1)*sizeof(float))
#define MAX_HESSIAN_MATRIX_SIZE ((6 +  MAX_NUM_OF_LIG_TORSION + MAX_NUM_OF_FLEX_TORSION)*(6 +  MAX_NUM_OF_LIG_TORSION + MAX_NUM_OF_FLEX_TORSION + 1) / 2)
#define MAX_NUM_OF_LIG_PAIRS 4096
#define MAX_NUM_OF_BFGS_STEPS 64
#define MAX_NUM_OF_RANDOM_MAP 1000 // not too large (stack overflow!)
#define GRIDS_SIZE 17

#define MAX_NUM_OF_GRID_MI 128//55
#define MAX_NUM_OF_GRID_MJ 128//55
#define MAX_NUM_OF_GRID_MK 128//81

//#define GRID_MI 65//55
//#define GRID_MJ 71//55
//#define GRID_MK 61//81
#define MAX_P_DATA_M_DATA_SIZE 256
//#define MAX_NUM_OF_GRID_ATOMS 130
#define FAST_SIZE 2051
#define SMOOTH_SIZE 2051
#define MAX_CONTAINER_SIZE_EVERY_WI 5
#define EL_TYPE_H_CL 0


typedef struct {
	float data[GRIDS_SIZE];
} affinities_cl;

typedef struct {
	int types[4];
	float coords[3];
} atom_cl;

typedef struct {
	atom_cl atoms[MAX_NUM_OF_ATOMS];
} grid_atoms_cl;

typedef struct {
	float coords[MAX_NUM_OF_ATOMS][3];
} m_coords_cl;

typedef struct {
	float coords[MAX_NUM_OF_ATOMS][3];
} m_minus_forces;

typedef struct  { // namely molec_struc
	float position		[3];
	float orientation	[4];
	float lig_torsion	[MAX_NUM_OF_LIG_TORSION];
	float flex_torsion	[MAX_NUM_OF_FLEX_TORSION];
	float coords		[MAX_NUM_OF_ATOMS][3];
	float lig_torsion_size;
	float e;
} output_type_cl;

typedef struct  { // namely change_struc
	float position		[3];
	float orientation	[3];
	float lig_torsion	[MAX_NUM_OF_LIG_TORSION];
	float flex_torsion	[MAX_NUM_OF_FLEX_TORSION];
	float lig_torsion_size;
} change_cl;


typedef struct { // depth-first order
	int		atom_range		[MAX_NUM_OF_RIGID][2];
	float	origin			[MAX_NUM_OF_RIGID][3];
	float	orientation_m	[MAX_NUM_OF_RIGID][9]; // This matrix is fixed to 3*3
	float	orientation_q	[MAX_NUM_OF_RIGID][4];
	
	float	axis			[MAX_NUM_OF_RIGID][3]; // 1st column is root node, all 0s
	float	relative_axis	[MAX_NUM_OF_RIGID][3]; // 1st column is root node, all 0s
	float	relative_origin	[MAX_NUM_OF_RIGID][3]; // 1st column is root node, all 0s
	
	int		parent			[MAX_NUM_OF_RIGID]; // every node has only 1 parent node
	bool	children_map	[MAX_NUM_OF_RIGID][MAX_NUM_OF_RIGID]; // chidren_map[i][j] = true if node i's child is node j
	int		num_children;
	
} rigid_cl;

typedef struct {
	int type_pair_index	[MAX_NUM_OF_LIG_PAIRS];
	int a				[MAX_NUM_OF_LIG_PAIRS];
	int b				[MAX_NUM_OF_LIG_PAIRS];
	int num_pairs;
} lig_pairs_cl;

typedef struct {
	lig_pairs_cl pairs;
	rigid_cl rigid;
	int begin;
	int end;
} ligand_cl;

typedef struct {
	int		int_map		[MAX_NUM_OF_RANDOM_MAP];
	float	pi_map		[MAX_NUM_OF_RANDOM_MAP];
	float	sphere_map	[MAX_NUM_OF_RANDOM_MAP][3];
} random_maps;

typedef struct {
	atom_cl atoms[MAX_NUM_OF_ATOMS];
	m_coords_cl m_coords;
	m_minus_forces minus_forces;
	ligand_cl ligand;
	int m_num_movable_atoms;
} m_cl;

typedef struct {
	float m_init[3];
	float m_range[3];
	float m_factor[3];
	float m_dim_fl_minus_1[3];
	float m_factor_inv[3];
	int m_i;
	int m_j;
	int m_k;
	float m_data [(MAX_NUM_OF_GRID_MI) * (MAX_NUM_OF_GRID_MJ) * (MAX_NUM_OF_GRID_MK) * 8];
} grid_cl;

typedef struct {
	int atu;
	float slope;
	grid_cl grids[GRIDS_SIZE];
} ig_cl;

typedef struct {
	float fast[FAST_SIZE];
	float smooth[SMOOTH_SIZE][2];
	float factor;
} p_m_data_cl;

typedef struct {
	float m_cutoff_sqr;
	int n;
	float factor;
	p_m_data_cl m_data[MAX_P_DATA_M_DATA_SIZE];
} p_cl;

typedef struct  {
	int max_steps;
	float average_required_improvement;
	int over;
	int ig_grids_m_data_step;
	int	p_data_m_data_step;
	int	atu;
	int m_num_movable_atoms;
	float slope;
	float epsilon_fl;
	float epsilon_fl2;
	float epsilon_fl3;
	float epsilon_fl4;
	float epsilon_fl5;
}variables_bfgs;

typedef struct {
	output_type_cl container[MAX_CONTAINER_SIZE_EVERY_WI];
	int current_size;
}out_container;

