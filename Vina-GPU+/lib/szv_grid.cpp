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

#include "szv_grid.h"
#include "brick.h"
/*
���幹�캯��szv_grid  ��ʼ����ֵ  Ϊ�ṹ����Ĳ��� m_init    m_range   m_data��ֵ
*/
szv_grid::szv_grid(const model& m, const grid_dims& gd, fl cutoff_sqr) : m_data(gd[0].n, gd[1].n, gd[2].n) {
	vec end;
	VINA_FOR_IN(i, gd) {
		m_init[i] = gd[i].begin;
		end   [i] = gd[i].end;
	}//��grid_dim��double����begin  ��end   ��ֵ
	m_range = end - m_init;//Ϊm_range��ֵ

	const sz nat = num_atom_types(m.atom_typing_used());	//atom_typing_usedö�����Ͳ��� //num_atom_types �õ�ö�ٶ�Ӧ����ֵ
	

 
//��get�����õ���t��Ӧ����ֵ uint  <   num_atom_types�����õ���t��Ӧ����ֵ uint  ����
//�������龭������brick_distance_sqr�õ����� double < cutoff_sqr
//��ô�� relevant_indexes��߼�Ԫ��i 
	szv relevant_indexes;  //uint�͵�����
	VINA_FOR_IN(i, m.grid_atoms) {//grid_atoms  �ṹ��atom���͵�����   
		const atom& a = m.grid_atoms[i];
		if(a.get(m.atom_typing_used()) < nat && brick_distance_sqr(m_init, end, a.coords) < cutoff_sqr)//get   �õ�t��Ӧ����ֵuint      atom_typing_used����ö����t��������
			relevant_indexes.push_back(i);
	}
	/*ѭ��ִ�� 
	x y zΪ    000   ��  a=grid_atoms��relevant_indexes0��  �ж�  �������龭������brick_distance_sqr���Ƿ�<cutoff_sqr
	                                                    m_data[i + m_i*(j + m_j*k)]�����Ԫ��relevant_indexes0
			   000   ��  a=grid_atoms��relevant_indexes1��     �ж�  �������龭������brick_distance_sqr���Ƿ�<cutoff_sqr
														m_data[i + m_i*(j + m_j*k)]�����Ԫ��relevant_indexes1        ��
			   000   ��  a=grid_atoms��relevant_indexes2��     �ж�  �������龭������brick_distance_sqr���Ƿ�<cutoff_sqr
														m_data[i + m_i*(j + m_j*k)]�����Ԫ��relevant_indexes2        ��
	x y zΪ	   001  .....
			   002  .....
			   00m_k  ....
			   m_im_jm_k....

   */
	VINA_FOR(x, m_data.dim0())
	VINA_FOR(y, m_data.dim1())
	VINA_FOR(z, m_data.dim2()) {
		VINA_FOR_IN(ri, relevant_indexes) {
			const sz i = relevant_indexes[ri];
			const atom& a = m.grid_atoms[i];
			if(brick_distance_sqr(index_to_coord(x, y, z), index_to_coord(x+1, y+1, z+1), a.coords) < cutoff_sqr)
				m_data(x, y, z).push_back(i);
		}
	}
}
/*����double ����average_num_possibilities   
* 
ѭ��ִ��  xyz��000  ��m_im_jm_k     counter += m_data(x, y, z).size()
����double����    counter /��m_i*m_j*m_k��

*/
fl szv_grid::average_num_possibilities() const {
	sz counter = 0;
	VINA_FOR(x, m_data.dim0())
	VINA_FOR(y, m_data.dim1())
	VINA_FOR(z, m_data.dim2()) {
		counter += m_data(x, y, z).size();
	}
	return fl(counter) / (m_data.dim0() * m_data.dim1() * m_data.dim2());
}
/*����ֻ������possibilities  ����uint������    ���볤��3������  coords
* inter     ����3uint��������  index
* function  ����index[i]��ֵ    ����   m_data[index[0] + m_i*(index[1] + m_j*index[2])]      m_data��array3d<szv>����
*/
const szv& szv_grid::possibilities(const vec& coords) const {
	boost::array<sz, 3> index;
	VINA_FOR_IN(i, index) {//ѭ��3��    
		assert(coords[i] + epsilon_fl >= m_init[i]);//�ж�����  epsilon_fl  ��2.2204460492503131e-016
		assert(coords[i] <= m_init[i] + m_range[i] + epsilon_fl);
		const fl tmp = (coords[i] - m_init[i]) * m_data.dim(i) / m_range[i];
		index[i] = fl_to_sz(tmp, m_data.dim(i) - 1);//��0<tmp<m_data.dim(i) - 1�� ǿ�����β����
	}
	return m_data(index[0], index[1], index[2]);//m_data��array3d<szv>����  ����m_data[index[0] + m_i*(index[1] + m_j*index[2])]
}

/*input     uint  i j  k
inter   index  tmp   ����3������
output  ����3������tmp
function  tmp[i] = m_init[i] + m_range[i] * index[i] / m_data.dim(i);    i= 0��1��2
����   index(i)= i  j   k   dim(i)=m_i;m_j;m_k;
*/
vec szv_grid::index_to_coord(sz i, sz j, sz k) const {
	vec index(i, j, k);
	vec tmp;
	VINA_FOR_IN(n, tmp) 
		tmp[n] = m_init[n] + m_range[n] * index[n] / m_data.dim(n);
	return tmp;
}
/*����grid_dims���͵ĺ���  szv_grid_dims    ����  ֻ��grid_dims���͵Ĳ���  gd   ���  grid_dims���͵�tmp
inter�� ����grid_dims���Ͳ���tmp   ��߰�������grid_dim���ͽṹ��
function����tmp[i]�Ĳ���begin   end    n��ֵ  Ȼ�󷵻�tmp
*/
grid_dims szv_grid_dims(const grid_dims& gd) {
	grid_dims tmp;    
	VINA_FOR_IN(i, tmp) {
		tmp[i].begin = gd[i].begin;//grid_dims�е�double���Ͳ���begin��end����
		tmp[i].end   = gd[i].end;
		fl n_fl = (gd[i].end - gd[i].begin) / 3; // 3A preferred size
		int n_int = int(n_fl);
		tmp[i].n     = (n_int < 1) ?  1 : sz(n_int);//��tmp��uint���Ͳ���n��ֵ
	}
	return tmp;
}

