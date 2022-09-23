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

#include "terms.h"
#include "brick.h"

/*		outout����tmp:double ������
		������һ��double������Ȼ�󽫽ṹ��conf_independent_inputs�е�
		num_tors��num_rotors��num_heavy_atoms��num_hydrophobic_atoms��ligand_max_num_h_bonds��num_ligands��ligand_lengths_sum
		��˳����ӵ�tmp����������ĩβһ��Ԫ�غ��档
*/

conf_independent_inputs::operator flv() const {
	flv tmp;                                     //double����
	tmp.push_back(num_tors);                    
	tmp.push_back(num_rotors);
	tmp.push_back(num_heavy_atoms);
	tmp.push_back(num_hydrophobic_atoms);
	tmp.push_back(ligand_max_num_h_bonds);
	tmp.push_back(num_ligands);
	tmp.push_back(ligand_lengths_sum);
	return tmp;
}

/*	
	input����m��model�ṹ�壬i��atom_index��
	output����acc��unsigned
	bonds��bond�ṹ������������atom�ṹ���е� bonds��
	b��bond�ṹ�壻
	a��atom �ṹ�壻

	a����b���ԭ���������grid_atoms�����е�����
	��a��������acc��1
*/
unsigned conf_independent_inputs::num_bonded_heavy_atoms(const model& m, const atom_index& i) const { // FIXME? - could be static, but I don't feel like declaring function friends
	unsigned acc = 0;
	const std::vector<bond>& bonds = m.get_atom(i).bonds;
	VINA_FOR_IN(j, bonds) {                   //for(j=0;j<bonds.size;j++)
		const bond& b = bonds[j];
		const atom& a = m.get_atom(b.connected_atom_index);
		if(!a.is_hydrogen())    //�Ƿ�Ϊ��(ad==6||ad==12);
			++acc;
	}
	return acc;
}

/*  ������ԭ�ӵĿ���ת����
	input����m��model�ṹ�壻i��atom_index�ṹ��
	output����acc��unsigned
	bonds��bond�ṹ������������atom�ṹ���е� bonds��
	b��bond�ṹ�壻
	a��atom �ṹ�壻

	a����b���ԭ���������grid_atoms�����е�����
	��a�����⡢b.rotatableΪtrue����num_bonded_heavy_atoms�ṹ����acc��1
	��acc��1
*/
// FIXME? - could be static, but I don't feel like declaring function friends
unsigned conf_independent_inputs::atom_rotors(const model& m, const atom_index& i) const { // the number of rotatable bonds to heavy ligand atoms
	unsigned acc = 0;
	const std::vector<bond>& bonds = m.get_atom(i).bonds;
	VINA_FOR_IN(j, bonds) {                          //for(i=0;i<bonds.size;i++)
		const bond& b = bonds[j];
		const atom& a = m.get_atom(b.connected_atom_index);
		if(b.rotatable && !a.is_hydrogen() && num_bonded_heavy_atoms(m, b.connected_atom_index) > 1) // not counting CH_3, etc
			++acc;
	}
	return acc;
}


/*  ������ԭ�ӵĿ���ת����
input����m��model�ṹ�壻
	num_tors��num_rotors��num_heavy_atoms��num_hydrophobic_atoms��
	ligand_max_num_h_bonds��ligand_lengths_sum��Щdouble�������ݶ���ʼ��Ϊ0
	num_ligandsΪligands��Ԫ�ظ�����vector_mutable�̳�vector��
*/
conf_independent_inputs::conf_independent_inputs(const model& m) {
	num_tors = 0;
	num_rotors = 0;
	num_heavy_atoms = 0;
	num_hydrophobic_atoms = 0;
	ligand_max_num_h_bonds = 0;
	num_ligands = m.ligands.size();
	ligand_lengths_sum = 0;

	VINA_FOR_IN(ligand_i, m.ligands) {                     //for(ligand_i=0;ligand_i< m.ligands.size;ligand_i++)
		const ligand& lig = m.ligands[ligand_i];
		ligand_lengths_sum += m.ligand_length(ligand_i);
		VINA_RANGE(i, lig.begin, lig.end) {                //for(i=lig.begin;i<lig.end;i++)
			const atom& a = m.atoms[i];
			if(a.el != EL_TYPE_H) {
				unsigned ar = atom_rotors(m, atom_index(i, false));

				num_tors += 0.5 * ar;

				if(ar > 2) num_rotors += 0.5;
				else num_rotors += 0.5 * ar;

				++num_heavy_atoms;
				if(xs_is_hydrophobic(a.xs))
					++num_hydrophobic_atoms;

				if(xs_is_acceptor(a.xs) || xs_is_donor(a.xs))
					++ligand_max_num_h_bonds;
			}
		}
	}
}


/*		
		outout����tmp: string ������
		������һ��string������Ȼ�󽫽ṹ��conf_independent_inputs�е�
		��num_tors������num_rotors������num_heavy_atoms������num_hydrophobic_atoms����
		��ligand_max_num_h_bonds������num_ligands������ligand_lengths_sum���ַ���
		��˳����ӵ�tmp����������ĩβһ��Ԫ�غ��档
		���ת����tmpʹ���Ϊdouble�������жϽ������tmp double����Ԫ�ظ����Ƿ���ת��ǰ��tmp string����Ԫ�ظ�����ȣ�
		�ǵĻ�����tmp��ת��ǰ�ģ�
*/

std::vector<std::string> conf_independent_inputs::get_names() const { // FIXME should probably be static
	std::vector<std::string> tmp;
	tmp.push_back("num_tors");
	tmp.push_back("num_rotors");
	tmp.push_back("num_heavy_atoms");
	tmp.push_back("num_hydrophobic_atoms");
	tmp.push_back("ligand_max_num_h_bonds");
	tmp.push_back("num_ligands");
	tmp.push_back("ligand_lengths_sum");
	VINA_CHECK(static_cast<flv>(*this).size() == tmp.size()); // FIXME?  ת��֮������Ĳ���ԭ�ж���
															   //����һ����������һ����ʱ������
	return tmp;
}

/*
num_tors��num_rotors��num_heavy_atoms��num_hydrophobic_atoms��
ligand_max_num_h_bonds��num_ligands��ligand_lengths_sum��Щdouble�������ݶ���ʼ��Ϊ0
*/
conf_independent_inputs::conf_independent_inputs() : 
	num_tors(0), num_rotors(0), num_heavy_atoms(0), 
	num_hydrophobic_atoms(0), ligand_max_num_h_bonds(0), num_ligands(0), 
	ligand_lengths_sum(0) {}

/*
	input����a��b��double������
	inter����n��unint
	output����acc��double��
	nΪa��b����Ԫ�ظ�������Сֵ

	demo��a{2��3��4��5}��a{1��1��2}����n=3
		  acc=2*1+3*1+2*4=13

*/
inline fl inner_product_shortest(const flv& a, const flv& b) {
	sz n = (std::min)(a.size(), b.size());
	fl acc = 0;
	VINA_FOR(i, n)                             //for(i=0;i<n;i++)
		acc += a[i] * b[i];
	return acc;
}


/*
       input����weights��double������include_internal��bool��
	   inter����eΪfactors�ṹ���е�e��double��������iΪfactors�ṹ���е�i��double��������
	   output���� tmp��double

	   demo��
	   e{2��3��4��5}��weights{1��1��2}��i{2��3}
	   if��include_internal=false����e*weights��
			 n=3
			 tmp=2*1+3*1+2*4=13��
	   else     ��i*weights��
	         n=2
	         tmp=2*1+3*1=5

*/
fl factors::eval(const flv& weights, bool include_internal) const {
	fl tmp = inner_product_shortest(e, weights);     //eΪfactors�ṹ���е�e��double������
	if(include_internal)
		tmp += inner_product_shortest(i, weights);   //iΪfactors�ṹ���е�i��double������
	return tmp;
}


// terms
/*
 input����enabled_only��bool��
 inter����tmp���ַ���������
	�ú���ͨ�������bool����enabled_only��term_set�ṹ���е�enabled��bool�������ж��Ƿ���tmp����ĩβ���term_set�ṹ����fun[i].name
		��enabled_onlyΪ0����enabled[i]���ڵ���1��
*/
std::vector<std::string> terms::get_names(bool enabled_only) const { // does not include conf-independent
	std::vector<std::string> tmp;
	/*��ȷ��enabled������Ԫ�ظ����Ƿ��fun����Ԫ����ͬ���ǵĻ�����������ֹ��
		������enabled_only��Ϊ0����enabled[i]���ڵ���1����tmp����ĩβ���fun[i].name��term�ṹ����name��*/
	distance_additive_terms.get_names(enabled_only, tmp);//һ���̳���distance_additive_terms�ṹ���term_set�ṹ��
	           usable_terms.get_names(enabled_only, tmp);//һ���̳���usable�ṹ���term_set�ṹ�壬�����洦�����tmp�ַ�������ĩβ���
	         additive_terms.get_names(enabled_only, tmp);//һ���̳���additive�ṹ���term_set�ṹ�壬�����洦�����tmp�ַ�������ĩβ���
	   intermolecular_terms.get_names(enabled_only, tmp);//һ���̳���intermolecular�ṹ���term_set�ṹ�壬�����洦�����tmp�ַ�������ĩβ���
	return tmp;
}


/*
	�ú�������distance_additive_terms�ṹ����fun����Ԫ�ظ�����usable_terms�ṹ����fun����Ԫ�ظ�����
	additive_terms��fun����Ԫ�ظ������ܺ�

*/
sz terms::size_internal() const {
	return distance_additive_terms.size() + usable_terms.size() + additive_terms.size();
	//size()���ú�������fun����Ԫ�ظ�����term_set�ṹ����fun�����������͵�ָ�������������vector���Ʋ���Ч�ʸ���
	//
}



/*
	input����enabled_only��bool����
	inter����enabled��bool��������
	��enabled_only=false����term_set�ṹ���е�enabled[i]>0����¶�fun[i]��Ԫ�ظ��������ۼӣ�i=0-fun����Ԫ�ظ�����
*/
sz terms::size_conf_independent(bool enabled_only) const { // number of parameters does not necessarily equal the number of operators
	sz acc = 0;
	VINA_FOR_IN(i, conf_independent_terms)    //for(i=0;i<conf_independent_terms.size(fun����Ԫ�ظ���);i++)
		if(!enabled_only || conf_independent_terms.enabled[i])//enabled_onlyΪ0����enabled[i]���ڵ���1
			acc += conf_independent_terms[i].size();         //conf_independent_terms[i]=fun[i]
	return acc;
}


/*

	�����������Ѱ�� distance_additive_terms.fun[i].cutoff��usable_terms.fun[i].cutoff��
	additive_terms.fun[i].cutoff�е����ֵ

*/

fl terms::max_r_cutoff() const {
	fl tmp = 0;
	tmp = (std::max)(tmp, distance_additive_terms.max_cutoff());//max_cutoff():Ѱ��fun[i].cutoff�е����ֵ
	tmp = (std::max)(tmp,            usable_terms.max_cutoff());
	tmp = (std::max)(tmp,          additive_terms.max_cutoff());
	return tmp;
}

/*
	input����m��model�ṹ�壬i��j��atom_index�ṹ�壬r��double��out��double����
	inter����a��b��atom�ṹ��
	a��b���ж�i��j�ĳ�Աin_grid�Ƿ�Ϊ�棬���򷵻�grid_atoms�����ĵڣ�����i��j�ĳ�Աi��j����ֵ��Ԫ�س������ã����򷵻�atoms�����ĵڣ�����i��j�ĳ�Աi��j����ֵ��Ԫ�س�������
*/
void terms::eval_additive_aux(const model& m, const atom_index& i, const atom_index& j, fl r, flv& out) const { // out is added to
	const atom& a = m.get_atom(i);
	const atom& b = m.get_atom(j);

	sz offset = 0;
	VINA_FOR_IN(k, distance_additive_terms)                      //for(k=0;k<term_set.fun.size().size;k++)
		if(r < distance_additive_terms[k].cutoff)                //(r<term_set<distance_additive_terms>.fun[i].cutoff)
			out[k] += distance_additive_terms[k].eval(a, b, r);  //get���������ǽṹ��atom_base��Ա��
															     //��Ϊatom_typing_use=XS�����Է���XS_TYPE_SIZE

	offset += distance_additive_terms.size();                    
	VINA_FOR_IN(k, usable_terms)								 //for(k=0;k<term_set<usable_terms>.fun.size().size;k++)
		if(r < usable_terms[k].cutoff)
			out[offset + k] += usable_terms[k].eval(a, b, r);    //�������out[k]�������

	offset += usable_terms.size();
	VINA_FOR_IN(k, additive_terms)							      //for(k=0;k<term_set<additive_terms>.fun.size().size;k++)                            
		if(r < additive_terms[k].cutoff)
			out[offset + k] += additive_terms[k].eval(m, i, j);   //�������out[k]�������
	
	VINA_CHECK(offset + additive_terms.size() == size_internal());//size_internal�����ú�������distance_additive_terms�ṹ����fun����Ԫ�ظ�����usable_terms�ṹ����fun����Ԫ�ظ�����
																//additive_terms��fun����Ԫ�ظ������ܺ�
}


/*�ⲿ���أ����ܰɣ�
	input����m��model�ṹ��
	output����tmp��double����

*/
flv terms::evale(const model& m) const {
	VINA_CHECK(m.ligands.size() == 1);      //ligandsΪģ�����Ϊligand�ṹ���vector_mutable�ṹ��
	VINA_CHECK(m.flex.size() == 0);
	VINA_CHECK(m.atoms.size() == m.num_movable_atoms()); // no inflex

	flv tmp(size(), 0);            //size()�ú�������distance_additive_terms�ṹ����fun����Ԫ�ظ�����
								   //usable_terms�ṹ����fun����Ԫ�ظ�����
								   //additive_terms��fun����Ԫ�ظ�����
								   //intermolecular_terms������Ԫ�ظ����ܺ�
	fl max_r_cutoff_sqr = sqr(max_r_cutoff()); //��Ѱ�� distance_additive_terms.fun[i].cutoff��usable_terms.fun[i].cutoff��
											   //additive_terms.fun[i].cutoff�е����ֵ��ƽ��

	grid_dims box = m.movable_atoms_box(0); // add nothing����3��Ԫ�ص�grid_dim�ṹ������
	vec box_begin = grid_dims_begin(box);   //������grid_dim�ṹ���е�begin��������
	vec box_end   = grid_dims_end  (box);   //������grid_dim�ṹ���е�end��������

	szv relevant_atoms;                     //unint����
	VINA_FOR_IN(j, m.grid_atoms)           //for(j=0;j< m.grid_atoms.size;j++)
		if(brick_distance_sqr(box_begin, box_end, m.grid_atoms[j].coords) < max_r_cutoff_sqr)
			/*
			demo��
			box_begin��1��2��3����	box_end��7��8��9����m.grid_atoms[j].coords��0��9��7��
			����
			closest(1��8��7)��closest��m.grid_atoms[j].coords��ӽ�box_begin��box_end�е�ֵ��
			����
			(1 - 0) ^ 2 + (8 - 9) ^ 2 + (7 - 7) ^ 2 = 2
			max_r_cutoff_sqΪdistance_additive_terms.fun[i].cutoff��usable_terms.fun[i].cutoff��
			additive_terms.fun[i].cutoff�е����ֵ��ƽ��
			*/
			relevant_atoms.push_back(j);    //���if��������j�ӵ�relevant_atoms�����������

	VINA_FOR(i, m.num_movable_atoms()) {                 //for(i=0;i<m.num_movable_atoms;i++)
		const vec& coords = m.coords[i];                 
		VINA_FOR_IN(relevant_j, relevant_atoms) {        //for(relevant_j=0;relevant_j<relevant_atoms.size;relevant_j++)
			const sz j = relevant_atoms[relevant_j];     //��relevant_atoms�е�Ԫ�ر���Ϊj
			const atom& b = m.grid_atoms[j];             //��grid_atoms�е�Ԫ�ر���Ϊb
			fl d2 = vec_distance_sqr(coords, b.coords); //��������vec�����data����Ԫ�ز��ƽ��

			if(d2 > max_r_cutoff_sqr) continue; // most likely scenario   if����ֱ����һ��VINA_FOR_IN
			fl d = std::sqrt(d2);               //��d2������
			eval_additive_aux(m, atom_index(i, false), atom_index(j, true), d, tmp);
		}
	}
	sz offset = size_internal();       //�ú�������distance_additive_terms�ṹ����fun����Ԫ�ظ�����usable_terms�ṹ����fun����Ԫ�ظ�����
										//additive_terms��fun����Ԫ�ظ������ܺ�
	VINA_FOR_IN(k, intermolecular_terms)
		tmp[offset + k] += intermolecular_terms[k].eval(m);

	VINA_CHECK(offset + intermolecular_terms.size() == tmp.size()); 
	return tmp;
}


/*�ڲ����أ����ܰɣ�
	input����m��model�ṹ��
	output����tmp��double����
*/

flv terms::evali(const model& m) const {
	VINA_CHECK(m.ligands.size() == 1);
	VINA_CHECK(m.flex.size() == 0);
	VINA_CHECK(m.atoms.size() == m.num_movable_atoms());

	flv tmp(size_internal(), 0);                         //size_internal�����ú�������distance_additive_terms�ṹ����fun����Ԫ�ظ�����usable_terms�ṹ����fun����Ԫ�ظ�����
	                                                     //	additive_terms��fun����Ԫ�ظ������ܺ�
	fl max_r_cutoff_sqr = sqr(max_r_cutoff());          //��Ѱ�� distance_additive_terms.fun[i].cutoff��usable_terms.fun[i].cutoff��
											            //additive_terms.fun[i].cutoff�е����ֵ��ƽ��
	const interacting_pairs& pairs = m.ligands.front().pairs;    //����interacting_pair�ṹ������
	VINA_FOR_IN(i, pairs) {                                      //for(i=0;i<pairs.size;i++)
		const interacting_pair& ip = pairs[i];                   //����interacting_pair�ṹ��
		fl d2 = vec_distance_sqr(m.coords[ip.a], m.coords[ip.b]);//��������vec�����data����Ԫ�ز��ƽ��
		if(d2 > max_r_cutoff_sqr) continue; // most likely scenario    if����ֱ����һ��ѭ��
		fl d = std::sqrt(d2);                                    //��d2������
		eval_additive_aux(m, atom_index(ip.a, false), atom_index(ip.b, false), d, tmp);
	}
	return tmp;
}


/*�ⲿ������ǿ��ֻ�е�������ϵͳ��Ҫ�õ��������
	input����m��model�ṹ��
	output����tmp��double����
*/
flv terms::evale_robust(const model& m) const {
	VINA_CHECK(m.ligands.size() == 1); // only single-ligand systems are supported by this procedure

	flv tmp(size(), 0);                             

	fl max_r_cutoff_sqr = sqr(max_r_cutoff());

	grid_dims box = m.movable_atoms_box(0);            //add nothing����3��Ԫ�ص�grid_dim�ṹ��
	vec box_begin = grid_dims_begin(box);	           //������grid_dim�ṹ���е�begin��������
	vec box_end = grid_dims_end(box);                  //������grid_dim�ṹ���е�end��������

	const sz n  = num_atom_types(m.atom_typing_used());/*n=
														ԭ��������
														EL�ͷ���11;
														AD�ͷ���20;
														XS�ͷ���17;
														SY�ͷ���18
														*/

	std::vector<atom_index> relevant_atoms;              //atom_index�ṹ������

	VINA_FOR_IN(j, m.grid_atoms) {
		const atom& a = m.grid_atoms[j];
		const sz t = a.get(m.atom_typing_used());        //ö����t��Ӧ����ֵ
		if(brick_distance_sqr(box_begin, box_end, a.coords) < max_r_cutoff_sqr && t < n) // exclude, say, Hydrogens
	          /*
	          demo��
	          box_begin��1��2��3����	box_end��7��8��9����m.grid_atoms[j].coords��0��9��7��
	          ����
	          closest(1��8��7)��closest��m.grid_atoms[j].coords��ӽ�box_begin��box_end�е�ֵ��
	          ����
	          (1 - 0) ^ 2 + (8 - 9) ^ 2 + (7 - 7) ^ 2 = 2
	          max_r_cutoff_sqΪdistance_additive_terms.fun[i].cutoff��usable_terms.fun[i].cutoff��
	          additive_terms.fun[i].cutoff�е����ֵ��ƽ��
	          */
			relevant_atoms.push_back(atom_index(j, true));//��(atom_index(j, true)����relevant_atoms�ṹ�������������
	}

	VINA_FOR_IN(j, m.atoms) {
		const atom& a = m.atoms[j];
		const vec& a_coords = m.coords[j];
		if(m.find_ligand(j) < m.ligands.size()) continue; // skip ligand atoms, add only flex/inflex  ����ֱ�ӽ�����һ��ѭ��
		const sz t = a.get(m.atom_typing_used());
		if(brick_distance_sqr(box_begin, box_end, a_coords) < max_r_cutoff_sqr && t < n) // exclude, say, Hydrogens
			relevant_atoms.push_back(atom_index(j, false));   //��(atom_index(j, false)���������relevant_atoms�ṹ�������������
	}

	VINA_FOR_IN(lig_i, m.ligands) {											
		const ligand& lig = m.ligands[lig_i];								
		VINA_RANGE(i, lig.begin, lig.end) {									
			const vec& coords = m.coords[i];								//��m�ṹ�������е�coords��������
			const atom& a = m.atoms[i];										//��m�ṹ�������е�atoms�ṹ�����
			const sz t = a.get(m.atom_typing_used());						//��m�ṹ�������е�m_atom_typing_used��������

			if(t < n) { // exclude, say, Hydrogens
				VINA_FOR_IN(relevant_j, relevant_atoms) {
					const atom_index& j = relevant_atoms[relevant_j];
					fl d2 = vec_distance_sqr(coords, m.atom_coords(j));    //��������vec�����data����Ԫ�ز��ƽ��
					if(d2 > max_r_cutoff_sqr) continue; // most likely scenario
					fl d = std::sqrt(d2);
					eval_additive_aux(m, atom_index(i, false), j, d, tmp);
				}
			}
		}
	}

	sz offset = size_internal();
	VINA_CHECK(intermolecular_terms.size() == 0);
	VINA_CHECK(offset + intermolecular_terms.size() == tmp.size());

	return tmp;
}


/*	��ȡ�ⲿ�ڲ�Ԫ��
	input����m��model�ṹ��
	output����tmp��factors�ṹ��
*/
factors terms::eval(const model& m) const {
	factors tmp;
	tmp.e = evale(m);     //�ⲿ���أ����ܰɣ�
	tmp.i = evali(m);     //�ڲ����أ����ܰɣ�
	return tmp;
}



/*
	input����in��conf_independent_inputs��it����������
	inter����conf_independent_terms��ģ���β�Ϊconf_independent��term_set�ṹ��
	output����x��double��
	if(conf_independent�ṹ����enabled[i]��ֵ����0)��i=0-conf_independent�ṹ����fun����Ԫ�ظ�����
		xΪ  x + read_iterator(i) * in.******   // read_iterator(i)Ϊ������ָ����һ��Ԫ�� ��******Ӧ����ȡ���ڵ������������ǰһ��everthing.cpp�еĺ���������
*/
fl terms::eval_conf_independent(const conf_independent_inputs& in, fl x, flv::const_iterator& it) const { // evaluates enabled only
	VINA_FOR_IN(i, conf_independent_terms)   //for(i=0;i<conf_independent_terms.size(fun����Ԫ�ظ���);i++ ��
		if(conf_independent_terms.enabled[i]) 
			x = conf_independent_terms[i].eval(in, x, it);     // x + read_iterator(i) * in.ligand_lengths_sum    
	return x;
}

/*
	input����v:double����
	��distance_additive_terms��usable_terms��additive_terms��intermolecular_terms�ṹ���е�enabled��ֵ���ڵ���1ʱ
	��ʱ����������Ԫ�����Σ���v�е�Ԫ�أ���ӵ�tmp������ĩβ��������ӣ�
	
*/
flv terms::filter_external(const flv& v) const {
	flv tmp;
	flv::const_iterator i = v.begin();      //��������i��ʼ��Ϊָ������v�����ĵ�һ��Ԫ��
	distance_additive_terms.filter(i, tmp); //��term_set�ṹ����enabled��ֵ���ڵ���1ʱ����ʱ��������Ԫ����ӵ�tmp������ĩβ
	           usable_terms.filter(i, tmp);
	         additive_terms.filter(i, tmp);
	   intermolecular_terms.filter(i, tmp);
	VINA_CHECK(i == v.end());              //�жϵ�����i�Ƿ��������v����
	return tmp;
}


/*
         input����v:double����
         ��distance_additive_terms��usable_terms��additive_terms��intermolecular_terms�ṹ���е�enabled��ֵ���ڵ���1ʱ
         ��ʱ����������Ԫ�����Σ���v�е�Ԫ�أ���ӵ�tmp������ĩβ���������
         ͬfilter_external������������һ��

*/
flv terms::filter_internal(const flv& v) const {
	flv tmp;
	flv::const_iterator i = v.begin();
	distance_additive_terms.filter(i, tmp);
	           usable_terms.filter(i, tmp);
	         additive_terms.filter(i, tmp);
	VINA_CHECK(i == v.end());
	return tmp;
}



/*		
		input����f��factors�ṹ�壻
		inter����tmp��factors�ṹ�壻
		tmp.e�õ�f.e���˺��������tmp.i�õ�f.i���˺����������󷵻�factors�ṹ��ĺ�����
*/
factors terms::filter(const factors& f) const {
	factors tmp;
	tmp.e = filter_external(f.e);
	tmp.i = filter_internal(f.i);
	return tmp;
}


/*
	����������ڴ�ӡterm_set�ṹ����fun[i].name�ͼ̳���conf_independent_terms�ṹ��
	��term_set�ṹ�����fun[i].name���Լ�fun����Ԫ�ظ����ܺͺ�cutoff�е����ֵ
*/

void terms::display_info() const {
	std::vector<std::string> enabled_names = get_names(true);//���ذ���term_set�ṹ����fun[i].name��string������
	std::cout << "Enabled terms: \n";  
	VINA_FOR_IN(i, enabled_names)                 //for(i=0;i<enabled_names.size;i++)
		std::cout << enabled_names[i] << '\n';   //��fun[i].name��ӡ����
	std::cout << '\n';

	std::vector<std::string> enabled_operators;
	conf_independent_terms.get_names(true, enabled_operators);   //out���ذ����̳���conf_independent_terms�ṹ��
																 //��term_set�ṹ�����fun[i].name��string������
	std::cout << "Enabled conf-independent operators: \n"; 
	VINA_FOR_IN(i, enabled_operators)						//for(i=0;i<enabled_operators.size;i++)
		std::cout << enabled_operators[i] << '\n';				//��fun[i].name��ӡ����
	std::cout << '\n';

	VINA_SHOW(size());                          //��ӡһ��size����=����Ԫ�ظ����ܺ�
												//�ú�������distance_additive_terms�ṹ����fun����Ԫ�ظ�����
								                //usable_terms�ṹ����fun����Ԫ�ظ�����
								                //additive_terms��fun����Ԫ�ظ�����
									            //intermolecular_terms��fun����Ԫ�ظ����ܺ�

	VINA_SHOW(size_internal());               //��ӡһ��size_internal()=fun����Ԫ�ظ������ܺ�
											  //�ú�������distance_additive_terms�ṹ����fun����Ԫ�ظ�����
										      //usable_terms�ṹ����fun����Ԫ�ظ�����
										      //additive_terms��fun����Ԫ�ظ������ܺ�

	VINA_SHOW(max_r_cutoff());              //��ӡһ��max_r_cutoff()=cutoff�����ֵ
											//distance_additive_terms.fun[i].cutoff��usable_terms.fun[i].cutoff��
										    //additive_terms.fun[i].cutoff�е����ֵ
}
