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

#include "ssd.h"

 /*  input��
			m����model�ṹ�壻p����precalculate�ṹ�壻ig����igrid�ṹ�壻out����output_type�ṹ�壻
			g����change�ṹ�壻v����vec�ṹ��

			��eval_deriv������������:
			1.������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
			2.����һ��double����e,����鿴����ig�е�eval_deriv����,igӦ����igrid���������
			3.��p,v[2], other_pairs, coords,minus_forces����eval_interacting_pairs_deriv����������ֵ�ۼӵ�e��
			4.����ligands��Ԫ��
					��p,v[0],ligands��i��Ԫ�صĳ�Աpairs,coord,minus_forces����eval_interacting_pairs_deriv����������ֵ�ۼӵ�e��
			5.��coords, minus_forces, g.ligands����ligands.derivative�У��ֱ��ligands������Ԫ��������ligands[i].derivative(coords, minus_forces, g.ligands[i])
			6.��coords, minus_forces, g.flex����flex.derivative�У��ֱ��flex������Ԫ��������flex[i].derivative(coords, minus_forces, g.flex[i])
			7.����e��double��

 */

// clean up
void ssd::operator()(model& m, const precalculate& p, const igrid& ig, output_type& out, change& g, const vec& v) const { // g must have correct size
	out.e = m.eval_deriv(p, ig, v, out.c, g);			 //out�ṹ���е�double����e
	fl factor = initial_factor;							 //ssd�ṹ���е�initial_factor����
	VINA_U_FOR(i, evals) {								 //unsigned=unsigned int
		if(factor < min_factor) break;					 //���initial_factor<min_factor
		output_type candidate(out);						 //out�ṹ�帴�Ƹ�candidate��candidate���ı䲻��ʹout�仯
		candidate.c.increment(g, -factor);				 //������Ա仯��λ�á�����Ť���Լ�������Ť�ؽ��б�׼��
		change candidate_g(g);							 //g�ṹ�帴�Ƹ�candidate_g��candidate_g���ı䲻��ʹg�仯
		candidate.e = m.eval_deriv(p, ig, v, candidate.c, candidate_g);
		if(candidate.e <= out.e) {                
			out = candidate;                              //candidate����out
			g = candidate_g;				              //candidate_g����g	
			factor *= up;                                 //factor=factor* up
		}
		else {
			factor *= down;                              //factor=factor* down
		}
	}
	out.coords = m.get_heavy_atom_movable_coords();      //��i=0������m_num_movable_atoms(model�ṹ����),�����ж�atoms��atom�ṹ����������Ա����el�Ƿ����0��
													     //������out.coordsβ�����model��Աcoords������i��Ԫ��
}
