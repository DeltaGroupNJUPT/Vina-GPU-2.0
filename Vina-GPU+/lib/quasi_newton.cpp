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

#include "quasi_newton.h"
#include "bfgs.h"



/*
      m����model�ṹ�壻p����precalculate�ṹ��ָ�룻ig����igrid�ṹ��ָ�룻
	  v����vec�ṹ�壻c����conf�ṹ�壻g����change�ṹ�壻
	 ��eval_deriv������������:
			1.������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
			2.����һ��double����e,����鿴����ig�е�eval_deriv����,igӦ����igrid���������
			3.��p,v[2], other_pairs, coords,minus_forces����eval_interacting_pairs_deriv����������ֵ�ۼӵ�e��
			4.����ligands��Ԫ��
					��p,v[0],ligands��i��Ԫ�صĳ�Աpairs,coord,minus_forces����eval_interacting_pairs_deriv����������ֵ�ۼӵ�e��
			5.��coords, minus_forces, g.ligands����ligands.derivative�У��ֱ��ligands������Ԫ��������ligands[i].derivative(coords, minus_forces, g.ligands[i])
			6.��coords, minus_forces, g.flex����flex.derivative�У��ֱ��flex������Ԫ��������flex[i].derivative(coords, minus_forces, g.flex[i])
			7.����e��double��
	  ��operator��������
			input����c��conf�ṹ�壬g��change�ṹ��
			���ڷ���eval_deriv������������ֵ


*/
struct quasi_newton_aux {
	model* m;
	const precalculate* p;
	const igrid* ig;           
	const vec v;
	quasi_newton_aux(model* m_, const precalculate* p_, const igrid* ig_, const vec& v_) : m(m_), p(p_), ig(ig_), v(v_) {}
	fl operator()(const conf& c, change& g) {
		const fl tmp = m->eval_deriv(*p, *ig, v, c, g);
		return tmp;                                    //e
	}
};



/*
		m����model�ṹ�壻p����precalculate�ṹ�壻ig����igrid�ṹ�壻out����output_type��
		g����change�ṹ�壻v����vec�ṹ��
		out.e=res=bfgs��������ֵf0
*/
void quasi_newton::operator()(model& m, const precalculate& p, const igrid& ig, output_type& out, change& g, const vec& v) const { // g must have correct size
	quasi_newton_aux aux(&m, &p, &ig, v);                                      //auxΪ�ع�quasi_newton_aux�ṹ��
	fl res = bfgs(aux, out.c, g, max_steps, average_required_improvement, 10); //f0
	out.e = res;                                                               
}

