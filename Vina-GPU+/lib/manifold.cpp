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

#include "manifold.h"
#include "recent_history.h"
#include "coords.h"
#include "quasi_newton.h"

const vec authentic_v(1000, 1000, 1000); 
/*�������������    
input��m  model    p��precalculate�ṹ��  ig��ig_widened��igrid�ṹ��   p_widened��p_widened�ṹ��  corner1��corner2������   generator:�����������
out��output_type
function������tmp   output_container    Ȼ��ǰ����ִ��parallel_mc::operator()  ��tmp��ֵ    �ж�tmp��Ϊ��    ���ҷ���tmp���׸�Ԫ��
*/
output_type manifold::operator()(model& m, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator) const {
	output_container tmp;
	this->operator()(m, tmp, p, ig, p_widened, ig_widened, corner1, corner2, generator); // call the version that produces the whole container
	VINA_CHECK(!tmp.empty());
	return tmp.front();
}
/*����ṹ��ָ�����ͺ���  select_rs     ֻ��
* input   output_container  �ṹ������  mf  ֻ��     rstart��uint      ����������� generator
*�жϽṹ���������ղ���Ԫ�ظ������ڵ���rstart     ��  ��������ṹ�������еĵ�i���ṹ���ַ     iΪ��0-��Ԫ�ظ���-1��֮�������� i 
* ���� ���ؿ�ָ��
*/
const conf* select_rs(const output_container& mf, sz rstart, rng& generator) { // clean up
	if(!mf.empty() && mf.size() >= rstart) { //�������ղ���Ԫ�ظ������ڵ���rstart
		sz i = random_sz(0, mf.size()-1, generator);//����һ����0-��Ԫ�ظ���-1��֮�������� i uint����
		return &(mf[i].c);//��������ṹ�������еĵ�i���ṹ���ַ
	}
	return NULL;//����ָ��Ϊ��
}
/*
* �жϽṹ�������е��ڲ�����仯�Ƿ�����   ���ҷ������һ�ν�����෴ֵ
input��conf���ͽṹ��C   output_container���ͽṹ������mf    ������������internal_too_close    double  ��exclusion�������ж����ƺ���internal_too_close����Ϊ�жϱ�׼��
out�� tmp  bool
�жϽṹ������mf[i]�еĽṹ��conf�еĽṹ������torsions�Ƿ��ligands�ṹ�������е�torsions���ƣ������жϽ����internal_too_close��ֵ
���һ���жϵĽ�������ƵĻ�   ����fause   ���򷵻�ture

*/
bool determine_iu(const conf& c, const output_container& mf, std::vector<bool>& internal_too_close, fl exclusion) { // clean up
	bool tmp = true;
	VINA_FOR_IN(i, mf) {
		bool t = c.internal_too_close(mf[i].c, exclusion);/* �ж�ligands�ṹ�������е�torsions�Ƿ��mf[i].conf.�ṹ�������е�torsions�����ǵĻ�����true�����򷵻�false */
		internal_too_close[i] = t;//�жϽ����������ֵ
		if(t) 
			tmp = false;//���һ���жϵĽ�������ƵĻ�   ����fause   ���򷵻�ture
	}
	return tmp;
}
/*

�жϽṹ�������е��ڲ�����仯���ⲿ����仯�Ƿ�����  �ǵĻ�����false  ���򷵻�ture
input  c  conf�ṹ��   mf  output_container  ���ṹ������   internal_too_close  bool����  ���ڲ��������Ƶ��жϽ���� exclusion��scale�ṹ��
out    bool 

external_too_close���ж��ⲿƥ��仯�Ƿ�����   �����ж�ligands[i].rigid�е�position��orientation�Ƿ��c.ligands[i].rigid�е�position��orientation���ƣ����ǵĻ�����false��
�����ж�flex[i]�е�torsions������c.flex[i]�е�torsions�����Ƿ����ƣ����ǵĻ�����false��
Ҫ����󷵻�true
*/

bool conf_is_legal(const conf& c, const output_container& mf, const std::vector<bool>& internal_too_close, const scale& exclusion) {
	assert(mf.size() == internal_too_close.size());// �ж�mf�ṹ��������Ԫ�ظ����Ƿ��internal_too_close������Ԫ�ظ������
	VINA_FOR_IN(i, mf)
		if(internal_too_close[i] && c.external_too_close(mf[i].c, exclusion))//�ж��ڲ��仯���ⲿ�仯�Ƿ�����
			return false;
	return true;
}
/*input   internal_conf��conf�ṹ��ֻ��     mf�� �ṹ������  ֻ��    internal_too_close�������������� ֻ��   uniq�� bool    spread��exclusion �� scale�ṹ�� ֻ��     rp��double  rs�� conf�ṹ��ָ��ֻ��  num_attempts��uint   failed ��bool    generator�������������
out�� tmp  conf�ṹ�� 
function��  ������internal_conf���ɸ�tmp    
			ѭ��ִ��             num_attempts��  
						1����tmp�����ⲿ�仯�����ԣ�λ�á�����Ť�أ������ԣ�Ť�أ�
						2�����uniqΪ0��tmp�ڲ������ⲿ���嶼����  ����tmp   ��
						3�����ִ�е�(attempt + 1 >= num_attempts)  �򷵻�tmp  ��failed���ture
			���ӵ������򷵻�internal_conf
*/
conf generate_external(const conf& internal_conf, const output_container& mf, const std::vector<bool>& internal_too_close, bool uniq, const scale& spread, const scale& exclusion, fl rp, const conf* rs, unsigned num_attempts, bool& failed, rng& generator) {
	// FIXME if one ligand with no side chains, don't use the same (pos, orient) more than once 
	failed = false;
	VINA_U_FOR(attempt, num_attempts) {//for  ��attempt=0��attempt<num_attempts;attempt++��
		conf tmp(internal_conf);
		tmp.generate_external(spread, rp, rs, generator); //�����ⲿ�仯�����ԣ�λ�á�����Ť�أ������ԣ�Ť�أ���
		if(uniq || conf_is_legal(tmp, mf, internal_too_close, exclusion))//���uniqΪ0��tmp�ڲ������ⲿ���嶼����  ����tmp
			return tmp;
		if(attempt + 1 >= num_attempts) { 
			failed = true;
			return tmp;
		}
	}
	VINA_CHECK(false); 
	return internal_conf; // placating the compiler
}
/*input   double   from   to    progress
* out    double
* function���ж� from �Ƿ�< 2.2204460492503131e-016  
��   return  from   
����   return    from * e��(progress * log(to/from))�η�
*/
fl extrapolate_cap(fl from, fl to, fl progress) {
	if(from < epsilon_fl) return from;
	return from * std::exp(progress * std::log(to/from));
}


/*input   from to�� ֻ�� double����   progress   ��double
out   double����
function��  ��tmp[i]��ֵ������  
�ж� from[i] �Ƿ�< 2.2204460492503131e-016 
��  tmp[i]=from[i] 
��  tmp[i]=from[i] * e��(progress * log(to[i]/from[i]))�η�
*/
vec extrapolate_cap(const vec& from, const vec& to, fl progress) {
	vec tmp;
	VINA_FOR_IN(i, tmp)
		tmp[i] = extrapolate_cap(from[i], to[i], progress);
	return tmp;
}
/*
* ����:   phase  :  uint , out:����output_type�ṹ�� , m��������model, mf������output_container�ṹ������ ֻ��,  p�����ýṹ��precalculate ֻ��,  ig�����ýṹ��igrid ֻ��, corner2corner��init_manifold_factor��double , e_internal_stats�����ýṹ��manifold, par�����ýṹ��manifoldֻ��, rng& generator�������������)
  �����   void

function: ѭ��ִ��  step��0-num_steps
  ��
		1����������ĺ���out����candidate��Ȼ���candidate������confִ��generate_internal�����������������ڲ��仯��λ�á�����Ť�ء�
		2���жϽṹ�������е��ڲ�����仯�Ƿ����ƽ����������internal_too_close��   ���ҷ������һ�ν�����෴ֵ
		3��������ִ��seti����������c.ligands��Ԫ�أ���model��Աinternal_coords, coords, ligands[i].begin��ligand��Ա��, ligands[i].end��ligand��Ա�����뺯��c.ligands[i].rigid.apply���������coords   ������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords

	                       ���У�generate_internal�������������Ϊ��spread.torsion, rp, rs, generator��

								spread��scale���ͽṹ�壩������Ϊ��orner2corner�����룩 * manifold_factor,pi * manifold_factor,pi * manifold_factor��
								                         manifold_factor=init_manifold_factor(����) * std::exp(- par.manifold_lambda * step / par.num_steps)
								rp                     �� 1 - std::exp(std::log(1-par.max_prob) * phase/par.num_phases)

								rs                     �� mf�ṹ�������е� ��������ṹ���ַ

         4���ж�rigid_trials>0
		 5��ѭ��ִ��   i��0-rigid_trials
										 ��

											1������bool����  failed =false
											2������candidate2��candidateһ�����͵Ľṹ�壩���Ҹ����е�����conf��ֵ���ı�failed
																	��ֵ������


											3���жϣ� ���failed=ture  �� ��candidate_is_set=ture   ���� i+1 < rigid_trials ��  ��������ѭ��
												  ����ִ�к���sete������c.ligands��Ԫ�أ���model��Աinternal_coords, coords, ligands[i].begin��ligand��Ա��, ligands[i].end��ligand��Ա�����뺯��c.ligands[i].rigid.apply���������coords   ������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords

											4����candidate2��e��ֵ����ֵ������1.����һ��double����e,����鿴����ig�е�eval����,igӦ����igrid���������
																		 2.����ligands��Ԫ�� ��p,v[2],other_pairs, coords����eval_interacting_pairs����������ֵ�ۼӵ�e��

											5���жϣ� ���candidate_is_setΪ1����  candidate2.e < candidate.e  ��ô��ִ��
												   ���� ִ��  ��candidate_is_set=ture  candidate = candidate2 ��
										   ��

		6�����candidate_is_set    if��ture��
					                    ��    
					                      1������ 	fl e_internal = (3 * step < 2 * num_steps) ? 0 : m.evali(p, hunt_cap);
																	����  evali���ã�1.����һ��double����e=0
*																	                2.����ligands��Ԫ��   ��p,v[0],ligands��i��Ԫ�صĳ�Աpairs,internal_coords����eval_interacting_pairs����������ֵ�ۼӵ�e��
*																	                3.����e
					
					                       2����add�������������e_internal    ���ı�e_internal_stats�����룩���error_estimate_sqr��x_estimate 

																	this_error_sqr=��e_internal - x_estimate��^2;
																	
																	error_estimate_sqr=weight * this_error_sqr + (1 - weight) * error_estimate_sqr
																	x_estimate         = weight * x              + (1 - weight) * x_estimate;

											3���ı�candidate�еĲ���e    candidate.e += e_internal;
					
											4���жϣ�   if(step == 0 || candidate.e < out.e)   out = candidate;

						                     ��

	��
		
		
		*/



void manifold_phase(sz phase, output_type& out, model& m, const output_container& mf, const precalculate& p, const igrid& ig, fl corner2corner, fl init_manifold_factor, recent_history& e_internal_stats, const manifold& par, rng& generator) { // out.first is starting conf on input
	out.e = max_fl; // �ṹ�����e��ֵ

	const fl rp = 1 - std::exp(std::log(1-par.max_prob) * phase/par.num_phases); // 
	//rp=1-e��(log(1-max_prob)*phase/num_phases)�η�
	fl e_internal_cost = par.relative_pair_cost * m.num_internal_pairs();//e_internal_cost=relative_pair_cost*pairs������Ԫ�ظ���
	fl e_external_cost = par.relative_pair_cost * m.num_other_pairs() + m.num_movable_atoms();
	const sz rigid_trials = 1 +  fl_to_sz(1 + par.cost_factor * e_internal_cost / e_external_cost, 100);
	//ǰ��С��100����0 ��ǿ������   �����±ߵ��ж�ִ��
	unsigned num_steps = (std::max)(unsigned(1), par.num_steps);

	VINA_U_FOR(step, num_steps) {//ѭ��num_steps
		const fl manifold_factor = init_manifold_factor * std::exp(- par.manifold_lambda * step / par.num_steps);
		//�����±ߵ�scale�ṹ��ĳ�ʼ��
		vec hunt_cap; hunt_cap = extrapolate_cap(par.hunt_cap, authentic_v, fl(step) / num_steps);

		const scale spread(corner2corner * manifold_factor,
			               pi            * manifold_factor,
						   pi            * manifold_factor);//����scale�ṹ��  ����ʼ��   �õ���forѭ��������manifold_factor

		sz rstart = fl_to_sz(par.rstart_fraction * par.num_phases, par.num_phases);
		//���ڹ��ɽṹ��ָ��rs    ���ں����ж�ִ���õ�
		const conf* rs = select_rs(mf, rstart, generator);
		//����mf�ṹ�������е� ��������ṹ���ַ
		output_type candidate = out;
		candidate.c.generate_internal(spread.torsion, rp, rs, generator); // ���������ڲ��仯��λ�á�����Ť�ء�

		std::vector<bool> internal_too_close(mf.size());//��mfһ�����ȵ�����
		bool uniq = determine_iu(candidate.c, mf, internal_too_close, par.exclusion.torsion); // internal_too_close is output
		//�жϽṹ�������е��ڲ�����仯�Ƿ�����   ���ҷ������һ�ν�����෴ֵ
		m.seti(candidate.c);
/*SETI
* ����model�ĳ�Ա����seti
* input:conf�ೣ������c
* ouput:void
* ���ã���model��Աatoms, internal_coords, ��c��Աligands����ligands�ĳ�Ա����set_conf    (*this)[i].set_conf(atoms, coords, c[i]);
*       ������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���internal_coords
*/



		VINA_CHECK(rigid_trials > 0);
		bool candidate_is_set = false;//
		VINA_FOR(i, rigid_trials) {
			bool failed = false;
			output_type candidate2(generate_external(candidate.c, mf, internal_too_close, uniq, spread, par.exclusion, rp, rs, par.num_attempts, failed, generator), 0);
			//���Ըı�  generate_external�� failed   ���ҷ���conf�ṹ��        ��������  ��ʼ��   candidate2
			if(failed && 
				(candidate_is_set || i+1 < rigid_trials)) // candidate is set or a chance to set it remains
				continue;//���failed=ture  �� ��candidate_is_set=ture   ���� i+1 < rigid_trials ��  ��������ѭ��
			m.sete(candidate2.c);

			/*SETE
* ����model�ĳ�Ա����sete
* input:conf�ೣ������c
* ouput:void
* ���ã�����c.ligands��Ԫ�أ���model��Աinternal_coords, coords, ligands[i].begin��ligand��Ա��, ligands[i].end��ligand��Ա�����뺯��c.ligands[i].rigid.apply���������coords
* ������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords				
*/
			
			candidate2.e = m.evale(p, ig, hunt_cap);  

			/*
* model�ĳ�Ա����evale
* input:precalculate�ೣ������p��igrid�ೣ������ig��vec�ೣ������v
* output:double������
* ���ã�1.����һ��double����e,����鿴����ig�е�eval����,igӦ����igrid���������
*		2.����ligands��Ԫ��
*				��p,v[2],other_pairs, coords����eval_interacting_pairs����������ֵ�ۼӵ�e��
*		3.����e 
*/
			if(!candidate_is_set || candidate2.e < candidate.e) {
				candidate_is_set = true;
				candidate = candidate2; // this works because internal conf is the same in both // FIXME inefficient, all this
			}//���candidate_is_setΪ1����  candidate2.e < candidate.e  ��ô��ִ��   ����ִ��  ��candidate_is_set=ture  candidate = candidate2 ��
		}

		VINA_CHECK(candidate_is_set);
		if(true) {
			fl e_internal = (3 * step < 2 * num_steps) ? 0 : m.evali(p, hunt_cap);
			/*EVALI
* model�ĳ�Ա����evali
* input:precalculate�ೣ������p��vec�ೣ������v
* output:double������
* ���ã�1.����һ��double����e=0
*		2.����ligands��Ԫ��
*				��p,v[0],ligands��i��Ԫ�صĳ�Աpairs,internal_coords����eval_interacting_pairs����������ֵ�ۼӵ�e��
*		3.����e
*/
			e_internal_stats.add(e_internal);//�ı�e_internal_stats���error_estimate_sqr��x_estimate
			/*add��
	input����x:double��
	inter����this_error_sqr��double��
		this_error_sqr=��x - x_estimate��^2;
		error_estimate_sqr=weight * this_error_sqr + (1 - weight) * error_estimate_sqr
		x_estimate         = weight * x              + (1 - weight) * x_estimate;
*/
			candidate.e += e_internal;//�ı�e
			if(step == 0 || candidate.e < out.e)
				out = candidate;
		}
	}
}
/*input��	 m,					����model���ͽṹ��
			out,				���ýṹ������output_container
			p,                  ���ýṹ��precalculate   ֻ��
			ig_widened�� ig,    ���ýṹ��igrid   ֻ��
			p_widened,          ���ýṹ��precalculate   ֻ��
			corner2 corner1,    ��������   ֻ��
			rng& generator      �����������
	out��  void

	function��
				1������manifold���ͽṹ��tuning_par   ����authentic_v��ֵ���������hunt_cap
				2������corner2corner  ���ڸ�manifold_phase������
				3������conf_size����s��change���ͽṹ��gS   s��ֵ��g
				4������recent_history����e_internal_stats  ����ֵ
	            5��ִ��ѭ��  phase��0-np
							��      1������output_type����tmp����ʼ��
									2��tmp�������confִ�к���randomize������(corner1, corner2, generator)�����ɸ��Ե����λ�á�����Ť�غ������е����Ť��
									3��ִ��manifold_phase����
									4���жϣ�   ���use_ssd���ṹ��manifold�Ĳ�����Ϊ1����ô�Խṹ��ssd_parִ�����غ���   �ı�tmp��e����  tmp.e=res=bfgs��������ֵf0
											   ����ѽṹ��ssd_par��evals������ֵ�� quasi_newton_par��quasi_newton���͵Ľṹ�壩�����max_steps����
											   Ȼ���quasi_newton_parִ�����غ���   �ı�tmp��e����  tmp.e=res=bfgs��������ֵf0
									5��tmp������confִ�к���set          ������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
									6����tmp��coords��ֵ    ִ��get_heavy_atom_movable_coords����������
									7�����ƽṹ��find_closest��closest_rmsd
									8���жϣ�//�� closest_rmsd��һ������<out��Ԫ�ظ����� �ڶ�������<min_rmsd   ��tmp.e < out[closest_rmsd.first].e
											//��out[closest_rmsd.first] = tmp
											//������out�ĺ������tmp
							��
						6��	����output_container���� null_arrayֻ��
						7���ж�out��Ϊ��    �������������   �ж�������ȷ

	*/


void manifold::operator()(model& m, output_container& out, const precalculate& p, const igrid& ig, const precalculate& p_widened, const igrid& ig_widened, const vec& corner1, const vec& corner2, rng& generator) const {
	manifold tuning_par;
	tuning_par.hunt_cap = authentic_v;//��tuning_par��������ֵ

	fl corner2corner = 0.5 * (corner2 - corner1).norm();//  norm��sqrt��sqr(data[0]) + sqr(data[1]) + sqr(data[2])��
	//���ڸ�manifold_phase��ֵ
	conf_size s = m.get_size();//����S   ���ڸ�  g��ֵ
/*get_size
* ���ã� 1.���ligands��Ťת������������ʱtmp.ligands������
*       2.���flex��Ťת������������ʱtmp.flex������
*		3.����tmp��ֵ��  s
*/
	change g(s);//����change���ͽṹ�岢��ʼ��

	sz np = (std::max)(sz(1), num_phases);//����ѭ��  ȡ1��num_phases����
	recent_history e_internal_stats(0, 10, 10);//����ṹ��recent_history���ҳ�ʼ��    ���� ��manifold_phase��ֵ



	VINA_FOR(phase, np) {
		output_type tmp(s, 0);
		tmp.c.randomize(corner1, corner2, generator);
/*���ɸ��Ե����λ�á�����Ť�غ������е����Ť�أ�
ligands�ṹ�����������������position�ṹ���orientation��Ԫ���Լ�[-pi,pi]���������torsions��
flex�ṹ������������[-pi,pi]���������torsions
*/
		manifold_phase(phase, tmp, m, out, p_widened, ig_widened, corner2corner, 1, e_internal_stats, *this, generator); // 1 - initial manifold
		
		
		if(use_ssd)
			ssd_par(m, p, ig, tmp, g, authentic_v);//�������� �ı�tmp��e����  tmp.e=res=bfgs��������ֵf0
		else {
			quasi_newton quasi_newton_par; 
			quasi_newton_par.max_steps = ssd_par.evals;
			quasi_newton_par(m, p, ig, tmp, g, authentic_v);//�������� �ı�tmp��e����  tmp.e=res=bfgs��������ֵf0
		}


		m.set(tmp.c);
/*set
* ������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
*/
		tmp.coords = m.get_heavy_atom_movable_coords();//��tmp��coords��ֵ
/*
* ������get_heavy_atom_movable_coords����������
* function��1.����һ����ʱvecv����tmp
* 2.i��0������m_num_movable_atoms�������ж�atoms��i��Ԫ�صĳ�Աel�Ƿ����0��������tmpβ�����model��Աcoords������i��Ԫ��
* 3.����tmp
*/
		VINA_CHECK(tmp.coords.size() > 0);//�ж�coords�Ƿ�Ϊ��
		std::pair<sz, fl> closest_rmsd = find_closest(tmp.coords, out);
		//����find_closest
		

		//�� closest_rmsd��һ������<out��Ԫ�ظ����� �ڶ�������<min_rmsd   ��tmp.e < out[closest_rmsd.first].e
		//��out[closest_rmsd.first] = tmp
		//������out�ĺ������tmp
		if(closest_rmsd.first < out.size() && closest_rmsd.second < min_rmsd) { // have a very similar one
			if(tmp.e < out[closest_rmsd.first].e) { // the new one is better, apparently
				out[closest_rmsd.first] = tmp; // FIXME? slow
			}
		}
		else { // nothing similar
			out.push_back(new output_type(tmp));
		}
	}
	// final tunings
	const output_container null_array;//�ṹ������
	VINA_CHECK(!out.empty());
	out.sort();//��out�нṹ��Ĳ���������������
	VINA_CHECK(out.front().e <= out.back().e); // //��֤����   make sure the sorting worked in the correct order
}
