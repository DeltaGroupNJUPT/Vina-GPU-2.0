/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

#ifndef VINA_CACHE_H
#define VINA_CACHE_H

#include <string>
#include "igrid.h"
#include "grid.h"
#include "model.h"

//OpenCL Related
#include "commonMacros.h"
//#define MAX_NUM_OF_M_DATA_DIM0 99
//#define MAX_NUM_OF_M_DATA_DIM1 99
//#define MAX_NUM_OF_M_DATA_DIM2 99
//#define MAX_NUM_OF_M_DATA MAX_NUM_OF_M_DATA_DIM0*MAX_NUM_OF_M_DATA_DIM1*MAX_NUM_OF_M_DATA_DIM2
//#define MAX_M_DATA_MI 16
//#define MAX_M_DATA_MJ 16
//#define MAX_M_DATA_MK 16
//#define MAX_NUM_OF_EVERY_M_DATA_ELEMENT 512
//#define MAX_NUM_OF_TOTAL_M_DATA MAX_M_DATA_MI*MAX_M_DATA_MJ*MAX_M_DATA_MK*MAX_NUM_OF_EVERY_M_DATA_ELEMENT

//struct AtomConstant_gpu {
//	cl_mem atom_types;
//	cl_mem atom_coords;
//};
//#define MAX_NUM_OF_M_DATA 99
//end OpenCL Related

struct cache_mismatch {};                   
struct rigid_mismatch : public cache_mismatch {};     //�̳�cache_mismatch�ṹ��
struct grid_dims_mismatch : public cache_mismatch {};   //�̳�cache_mismatch�ṹ��
struct energy_mismatch : public cache_mismatch {};    //�̳�cache_mismatch�ṹ��


/*�̳�igrid�ṹ��
	scoring_function_version_��string��
	gd_��grid_dims�ṹ�壻
	slope_��double��
	atom_typing_used_��atom_type�ṹ���е�ö���ͱ���t ��
*/
struct cache : public igrid {
	cache(const std::string& scoring_function_version_, const grid_dims& gd_, fl slope_, atom_type::t atom_typing_used_);
	/*
	��������
	*/
	fl eval      (const model& m, fl v) const; // needs m.coords // clean up
	fl eval_deriv(      model& m, fl v) const; // needs m.coords, sets m.minus_forces // clean up
	std::vector<grid> get_grids()const;
	int get_atu()const;
	double get_slope() const;

#if 0 // no longer doing I/O of the cache
	void read(const path& name);         // can throw cache_mismatch���������
	void write(const path& name) const;  //д�������
#endif
	/*
	��������
	*/
	void populate(const model& m, const precalculate& p, const szv& atom_types_needed, bool display_progress = true);
	void populate_cl(const model& m, const precalculate& p, const szv& atom_types_needed, bool display_progress = true);
	bool m_data_check(const szv_grid ig)const;
public:
	std::string scoring_function_version;        //�ַ���
	atomv atoms; // for verification             //atom�ṹ������
	grid_dims gd;                                //grid_dims�ṹ��
	fl slope; // does not get (de-)serialized    //double
	atom_type::t atu;                            //atom_type�ṹ���е�ö���ͱ���t 
	std::vector<grid> grids;                     //grid�ṹ������
	/*
		���л�����save����;
		�����л�����load������
		ע����Ҫ����һ����Ҫ�ĺ꣬BOOST_SERIALIZATION_SPLIT_MEMBER 
	*/
	friend class boost::serialization::access;
	template<class Archive>
	void save(Archive& ar, const unsigned version) const;
	template<class Archive>
	void load(Archive& ar, const unsigned version);
	BOOST_SERIALIZATION_SPLIT_MEMBER()
};

#endif
