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

#include "naive_non_cache.h"
#include "curl.h"
/*  �����������pָ�븳��naive_non_cache�ṹ���е�pָ��
	input����p��precalculate�ṹ��ָ��
*/
naive_non_cache::naive_non_cache(const precalculate* p_) : p(p_) {}

/*
	input����m��model�ṹ�壬v��double
	inter����cutoff_sqr��this_e��r2��double��n��m_num_movable_atoms��t1��t2��unint
			a_coords��vec�ṹ�壬a��b��atom�ṹ��
	output����e��double

*/
fl naive_non_cache::eval(const model& m, fl v) const { // needs m.coords
	fl e = 0;
	const fl cutoff_sqr = p->cutoff_sqr();				//weighted_terms�ṹ�����ع���cutoff���������ķ���ֵcutoff_��ƽ��

	sz n = num_atom_types(p->atom_typing_used());		//precalculate�ṹ����ԭ�������ͣ�
														//EL�ͷ���11��AD�ͷ���20��XS�ͷ���17��SY�ͷ���18

	VINA_FOR(i, m.num_movable_atoms()) {				//for(i=0;m.m_num_movable_atoms<n;i++)
		fl this_e = 0;
		const atom& a = m.atoms[i];						 //atoms��atom�ṹ������������a
		sz t1 = a.get(p->atom_typing_used());			 // precalculate�ṹ����ԭ�������ͣ�
														 //el=11; ad=20; xs=17; sy=18
		if(t1 >= n) continue;                            //t1 >= n������һ��forѭ��
		
		
		const vec& a_coords = m.coords[i];               //coords��vec�ṹ������������a_coords

		VINA_FOR_IN(j, m.grid_atoms) {                   //for(i=0;i<grid_atoms.size;i++)
			const atom& b = m.grid_atoms[j];             //grid_atoms��atom�ṹ������������b
			sz t2 = b.get(p->atom_typing_used());        // precalculate�ṹ����ԭ�������ͣ�
														 //el=11; ad=20; xs=17; sy=18
			if(t2 >= n) continue;                        //t2 >= n������һ��forѭ��
			
			
			vec r_ba; r_ba = a_coords - b.coords;        //a_coords = { 2,2,3 };
														 //b.coords = { 2,2,2 };
														 //r_ba = a - b={0,0,1}
			fl r2 = sqr(r_ba);							 //r2=sqr(r_ba)=0^2+0^2+1^2=1
			if(r2 < cutoff_sqr) {
				sz type_pair_index = get_type_pair_index(p->atom_typing_used(), a, b);//��ȡ���Ͷ�����
				this_e +=  p->eval_fast(type_pair_index, r2);  //this_e=this_e+data(type_pair_index).fast[sz(factor * r2)}
			}
		}
			curl(this_e, v);           //���ж�this_e > 0����v < 0.1*max_fl(��ǰ��������double�ͱ������ֵ)
								       //���v<2.22045e-16 ,this_e=0������this_e=this_e*v / (v + this_e)
		e += this_e;                   //e=e+this_e
	}
	return e;
}

/*Added by Glinttsd
*/
std::vector<grid> naive_non_cache::get_grids()const {
	assert(false); // This function should not be called!
	std::vector<grid> g;
	return g;
};

int naive_non_cache::get_atu()const {
	assert(false); // This function should not be called!
	return 0;
}

double naive_non_cache::get_slope()const {
	assert(false); // This function should not be called!
	return 0;
}