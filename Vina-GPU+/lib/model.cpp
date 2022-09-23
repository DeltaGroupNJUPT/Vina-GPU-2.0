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

#include "model.h"
#include "file.h"
#include "curl.h"


/*
* ����ģ��get_atom_range,��ȡԭ�Ӿ���
* input:T���Ͷ���t���磺ligand�ࣩ
* output��atom_range��
* ���ã���һ���ڵ��ϣ�����һ��t�ĳ�Աnode(atom_range��),����t�ĳ�Աchildren����������Ԫ�صĳ�Աbegin����Сֵ��ֵ���ÿ����ĳ�Աbegin����������Ԫ�صĳ�Աend�����ֵ��ֵ���ÿ����ĳ�Աend
*        ����˵����ͼһ
*/

template<typename T>
atom_range get_atom_range(const T& t) {
	atom_range tmp = t.node;//����һ����ʱatom_range����tmp�����t�ĳ�Աnode
	//����t�ĳ�Աchildren����,�����������Ԫ�صĳ�Աbegin����Сֵ����ֵ��tmp�ĳ�Աbegin�������������Ԫ�صĳ�Աend�����ֵ����ֵ��tmp�ĳ�Աend
	VINA_FOR_IN(i, t.children) {
		atom_range r = get_atom_range(t.children[i]);
		if(tmp.begin > r.begin) tmp.begin = r.begin;
		if(tmp.end   < r.end  ) tmp.end   = r.end;
	}
	return tmp;
}


//�ṹ��branch_metrics��Ĭ�ϳ�ʼ��unsigned int���͵����ݳ�Աlength��corner2corner��Ϊ0
struct branch_metrics {
	sz length;
	sz corner2corner;
	branch_metrics() : length(0), corner2corner(0) {}
};


/*
* ����ģ��get_branch_metrics
* input:T���Ͷ���t
* output��branch_metrics��
* ���ã���ͼһ����
*/
template<typename T>
branch_metrics get_branch_metrics(const T& t) {
	//����һ����ʱbranch_metrics�����tmp�����淵��ֵ
	branch_metrics tmp;
	//�ж�t�ĳ�Աchildren�Ƿ�Ϊ�գ�Ϊ��ֱ�ӷ���Ĭ�ϳ�ʼ����branch_metrics����󣬷������ִ�����´��룬
	if(!t.children.empty()) {
		//����unsigned int��������corner2corner_max����ʼ��Ϊ0
		sz corner2corner_max = 0;
		//����unsigned int������������lengths
		szv lengths;

		//����t�ĳ�Աchildren����
		VINA_FOR_IN(i, t.children) {
			branch_metrics res = get_branch_metrics(t.children[i]);//�ݹ����
			if(corner2corner_max < res.corner2corner)//�ж�corner2corner_max�Ƿ�С�ڴ��ӷ�֧��corner2corner������ֵcorner2corner��corner2corner_max
				corner2corner_max = res.corner2corner;
			lengths.push_back(res.length + 1); //���ӷ�֧�ĳ��ȼ�һ��������lengths��   FIXME? weird compiler warning (sz -> unsigned)
		}

		//����������lengths������Ԫ��
		std::sort(lengths.begin(), lengths.end());

		//����lengths������ĩβԪ�ص����ã���ֵ��tmp��length��Ա
		tmp.length = lengths.back();

		//tmp��length��Ա������tmp��corner2corner��Ա
		tmp.corner2corner = tmp.length;
		//�������lengths�Ĵ�С���ڵ���2��tmp��corner2corner��Ա��ֵΪlengths���Ԫ�ص�����
		if(lengths.size() >= 2)
			tmp.corner2corner += lengths[lengths.size() - 1];//�����壬�Ҿ��ô˴�Ӧ��lengths.size() - 2����֪���ǲ��ǳ���д����
		//���tmp��corner2corner��ԱС��corner2corner_max����corner2corner_max��ֵ��ֵ��corner2corner
		if(tmp.corner2corner < corner2corner_max)
			tmp.corner2corner = corner2corner_max;
	}
	return tmp;
}

/*
* model�ĳ�Ա����ligand_longest_branch
* input:unsigned int����ligand_number
* output: model�ĳ�Աligands������ligand_number��Ԫ�ض�Ӧ��branch_metrics���length����
*/
sz model::ligand_longest_branch(sz ligand_number) const {
	return get_branch_metrics(ligands[ligand_number]).length;
}


/*
* model�ĳ�Ա����length
* input:unsigned int����ligand_number
* output: model�ĳ�Աligands������ligand_number��Ԫ�ض�Ӧ��branch_metrics���corner2corner����
*/
sz model::ligand_length(sz ligand_number) const {
	return get_branch_metrics(ligands[ligand_number]).corner2corner;
}

/*
* ligand�ĳ�Ա����set_range
* input:��
* output: ��
* ���ã���ȡligand�ĳ�Աnode��begin���������к���ĳ�Աnode��begin���ҵ���Сֵ��ֵ��begin
*       ��ȡligand�ĳ�Աnode��end��  �������к���ĳ�Աnode��end��  �ҵ����ֵ��ֵ��end
*/
void ligand::set_range() {
	atom_range tmp = get_atom_range(*this);
	begin = tmp.begin;
	end   = tmp.end;
}

/////////////////// begin MODEL::APPEND /////////////////////////

// FIXME hairy code - needs to be extensively commented, asserted, reviewed and tested

/*
* �ṹ�壺appender_info
* ��Ա��1.unsigend int����grid_atoms_size��m_num_movable_atoms��atoms_size
*       2.���캯��������
*           ��ʼ��ʱ����һ��model����
*           model����ĳ�Աgrid_atoms������С��ֵ��grid_atoms_size��
*           model����ĳ�Աm_num_movable_atoms��ֵ��m_num_movable_atoms��
*           model����ĳ�Աatoms������С��ֵ��atoms_size
*/
struct appender_info {
	sz grid_atoms_size;
	sz m_num_movable_atoms;
	sz atoms_size;

	appender_info(const model& m) : grid_atoms_size(m.grid_atoms.size()), m_num_movable_atoms(m.m_num_movable_atoms), atoms_size(m.atoms.size()) {}
};

/*
* �ࣺappender
* ��Ա��private:1.appender_info��a_info
*               2.appender_info��b_info
*               3.��������new_grid_index��
*							input:unsigned int����x
*							output:unsigned int����
*							���ã�����appender��Ĳ���ֵis_a��Ϊ�淵��x��Ϊ�ٷ���a_info�ĳ�Աgrid_atoms_size+x
*       
*		public: 1.bool����is_a
*			    2.���캯������������model����󣬷ֱ�󶨵�a_info��b_info��is_a��ʼ��Ϊ��
*               3.������ò���������
*							input:unsigned int����x
*							output:unsigned int����
*							���ã�����appender��Ĳ���ֵis_a��Ϊ���ж�x�Ƿ�С��a_info�����m_num_movable_atoms��Ϊ�淵��x
*																												Ϊ�ٷ���b_info�����m_num_movable_atoms+x
*															  Ϊ���ж�x�Ƿ�С��b_info�����m_num_movable_atoms��Ϊ�淵��a_info�����m_num_movable_atoms+x															
*                											  												    Ϊ�ٷ���a_info�����atoms_size+x
*               4.���ص��ò���������
* 							input:atom_index�����x
* 							output:atom_index��
* 							���ã���x������tmp���ж�tmp�ĳ�Աin_grid�Ƿ�Ϊ�棬������tmp�ĳ�Աi���ڽ�����tmp�ĳ�Աinew_grid_index��ķ���ֵ
*																			  ������tmp�ĳ�Աi���ڽ�����tmp�ĳ�Աoperator��ķ���ֵ
*                                 ����tmp
*				5.��������update
*							input:interacting_pair������ip
* 							output:��
* 							���ã���ip�ĳ�Աa����operator()���������µ�ip�ĳ�Աa
*                 				  ��ip�ĳ�Աb����operator()���������µ�ip�ĳ�Աb
*				6.����update
*							input:ʸ��vec������
* 							output:��
* 							���ã���
*				7.����update
*							input:ligand������lig
*							output:��
*							���ã���1��	1.����appender�����(*this)����transform������
*									    2.��ȡlig�ĳ�Աbegin��end������diff = end + begin
*									    3.������������Ķ�����ò�������������*thisʹ�õ��ò��������㣨*this)(begin)������ֵ��begin
*									    4.����end = begin + diff  
* 								  ��2�� 1.����lig�ĳ�Աnode�ĳ�Աbegin��end���������к���ĳ�Աnode�ĳ�Աbegin��end
* 									    2.�����ȡ��ʽΪ���磺��һ��lig.children[i]�������lig.children[i][j],������lig.children[i][j][k]........�ȵȣ�ֱ�����еĺ��û�к���Ϊֹ
* 									    3.���и�����·�ʽ����ͬ����һ��˵������������ͬ
* 									    4.����diff= end - begin;��begin���뱾appender����ò�������������ȡ��begin;���¼���end = begin + diff
* 								  ��3������lig��Աpairs����pairs��Ԫ������update(interacting_pair& ip)����ȡ���º��pairs
* 								  ��4������lig��Աcont����cont��Ԫ������update(parsed_line& p)����ȡ���º��cont
*				8.����update
*							input:residue������r
* 							output:��
* 							���ã�1.����r�ĳ�Աnode�ĳ�Աbegin��end���������к���ĳ�Աnode�ĳ�Աbegin��end
* 							      2.�����ȡ��ʽΪ���磺��һ��r.children[i]�������r.children[i][j],������r.children[i][j][k]........�ȵȣ�ֱ�����еĺ��û�к���Ϊֹ
* 								  3.���и�����·�ʽ����ͬ����һ��˵������������ͬ
* 							      4.����diff= end - begin;��begin���뱾appender����ò�������������ȡ��begin;���¼���end = begin + diff
*               9.����update
*							input:parsed_line������p
*                           output:��
*							���ã��ж�p��optional������Ƿ�Ϊ�գ���Ϊ����optional������Ԫ�����뺯��operator()�У�����ֵ���¸�ֵ��optional�����
*               10.����update
*							input:atom������a
* 							output:��
* 							����:����a�ĳ�Աbonds��Ԫ�أ���Ԫ�صĳ�Աconnected_atom_index�������ص�operator����������connected_atom_index
*				11.����ģ�壺append
*                  input:T���͵���������a��T���͵ĳ�����������b         
*				   output:��				
* 				   ���ã�1.��ȡa�����Ĵ�Сa_sz				
* 					     2.������b��Ԫ����ӵ�����aβ��				
* 				         3.��is_a��ֵ��Ϊ�棬��������aԭ����Ԫ�أ�����ÿ��Ԫ�ظ���T���δ�����Ӧ���ذ汾update����
* 				         4.��is_a��ֵ��Ϊ�٣���������a����ӵ�Ԫ�أ�����ÿ��Ԫ�ظ���T���δ�����Ӧ���ذ汾update����
*				12.����ģ�壺coords_append
* 				   input:T���͵���������a��T���͵ĳ�����������b         
* 				   output:��				
* 				   ���ã�1.����b�Ŀ���b_copy
*						 2.��is_a��ֵ��Ϊ�棬��������a��Ԫ�أ�����ÿ��Ԫ�ظ���T���δ�����Ӧ���ذ汾update����
*						 3.��is_a��ֵ��Ϊ�٣���������b_copy��Ԫ�أ�����ÿ��Ԫ�ظ���T���δ�����Ӧ���ذ汾update����
*					     4.��b_copy����ͷ��һ����Ԫ�أ�b_info��Աm_num_movable_atoms��С��Ŀ������a��ĳһԪ�أ��±�Ϊa_info��Աm_num_movable_atomsb_info��С��ǰ
* 					     5.��b_copy����ʣ��Ԫ�ز���a��β��
*/

class appender {
	appender_info a_info;
	appender_info b_info;
	sz new_grid_index(sz x) const {
		return (is_a ? x : (a_info.grid_atoms_size + x)); // a-grid_atoms spliced before b-grid_atoms
	}
public:
	bool is_a;

	appender(const model& a, const model& b) : a_info(a), b_info(b), is_a(true) {}

	sz operator()(sz x) const { // transform coord index
		if(is_a) {
			if(x < a_info.m_num_movable_atoms)  return x; // a-movable unchanged
			else                                return x + b_info.m_num_movable_atoms; // b-movable spliced before a-inflex
		}
		else {
			if(x < b_info.m_num_movable_atoms)  return x + a_info.m_num_movable_atoms; // a-movable spliced before b-movable
			else                                return x + a_info.atoms_size; // all a's spliced before b-inflex
		}
	}
	atom_index operator()(const atom_index& x) const { // transform atom_index
		atom_index tmp(x);
		if(tmp.in_grid) tmp.i = new_grid_index(tmp.i);
		else            tmp.i = operator()(tmp.i);
		return tmp;
	}

	// type-directed old -> new transformations
	void update(interacting_pair& ip) const {
		ip.a = operator()(ip.a);
		ip.b = operator()(ip.b);
	}
	void update(vec& v) const { // coordinates & forces - do nothing
	}
	void update(ligand& lig) const {
		lig.transform(*this); // ligand as an atom_range subclass
		transform_ranges(lig, *this);
		VINA_FOR_IN(i, lig.pairs)
			this->update(lig.pairs[i]);
		VINA_FOR_IN(i, lig.cont)
			this->update(lig.cont[i]); // parsed_line update, below
	}
	void update(residue& r) const {
		transform_ranges(r, *this);
	}
	void update(parsed_line& p) const {
		if(p.second)
			p.second = operator()(p.second.get());
	}
	void update(atom& a) const {
		VINA_FOR_IN(i, a.bonds) {
			bond& b = a.bonds[i];
			b.connected_atom_index = operator()(b.connected_atom_index); // atom_index transformation, above
		}
	}

	// ligands, flex, flex_context, atoms; also used for other_pairs
	template<typename T>
	void append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbbbbbbb
		sz a_sz = a.size();
		vector_append(a, b);

		is_a = true;
		VINA_FOR(i, a_sz)
			update(a[i]);

		is_a = false;
		VINA_RANGE(i, a_sz, a.size())
			update(a[i]);
	}

	// internal_coords, coords, minus_forces, atoms
	template<typename T>
	void coords_append(std::vector<T>& a, const std::vector<T>& b) { // first arg becomes aaaaaaaabbbbbbbbbaab
		std::vector<T> b_copy(b); // more straightforward to make a copy of b and transform that than to do piecewise transformations of the result

		is_a = true;
		VINA_FOR_IN(i, a)
			update(a[i]);

		is_a = false;
		VINA_FOR_IN(i, b_copy)
			update(b_copy[i]);

		// interleave 
		typedef typename std::vector<T>::const_iterator const cci;
		cci b1 = b_copy.begin();
		cci b2 = b_copy.begin() + b_info.m_num_movable_atoms;
		cci b3 = b_copy.end();

		a.insert(a.begin() + a_info.m_num_movable_atoms , b1 , b2);
		a.insert(a.end()                                , b2 , b3);
	}
};


/*
* ����model�ĳ�Ա����append
* input:model�ೣ������m
* ouput:��
* ���ã�1.��鱾model������m_atom_typing_used�������model������m_atom_typing_used�Ƿ���ȣ�������򱨴���ֹ����ִ��
*       2.�ñ�model�����������model����󴴽�һ��appender��t
*		3.������model������Աother_pairs������Ԫ����ӵ���model������Աother_pairs����β����������model.append�����ĸ����㷨�Ա�model������Աother_pairs������Ԫ�ؽ��и���	
*		4.i��0ѭ��������model��Աatoms�Ĵ�С-1��
*				j��0ѭ����������model��Աatoms�Ĵ�С-1��
*						���i���ڵ��ڱ�model��Աm_num_movable_atoms����j���ڵ�������model��Աm_num_movable_atoms����ֱ��ִ����һ��ѭ��
*						����model��Աatoms�ĵ�i��Ԫ�ذ󶨵�a������model��Աatoms�ĵ�j��Ԫ�ذ󶨵�b
*						����model��Աm_atom_typing_used����a��get��������÷���ֵt1����m_atom_typing_used��ӦСд��ʽ��Ӧ��unsignd int��
*						������model��Աm_atom_typing_used����b��get��������÷���ֵt2����m_atom_typing_used��ӦСд��ʽ��Ӧ��unsignd int��
*						����model��Աm_atom_typing_used���뺯��num_atom_types��÷���ֵn
*						���ti��t2��С��n��
*									��t�ĳ�Աis_a����Ϊ�棬��i����t����ȡ����ֵnew_i;�ٽ�t�ĳ�Աis_a����Ϊ�٣���j����t����ȡ����ֵnew_j
*									���t1<=t2,�õ�type_pair_index = t1 + t2*(t2+1)/2������type_pair_index = t2 + t1*(t1+1)/2
*									����interacting_pair������type_pair_index��new_i, new_j��ʼ���������ö�����ӵ�other_pairsβ��
*
*		5.��鱾model�����Աminus_forces�Ĵ�С��coords�Ĵ�С�Ƿ���ȣ�������򱨴���ֹ�����ִ��
*		6.�������model�����Աminus_forces�Ĵ�С��coords�Ĵ�С�Ƿ���ȣ�������򱨴���ֹ�����ִ��
* 
*		7.������model������Աinternal_coords������Ԫ�ز��뵽��model������Աinternal_coordss�����ʵ�λ�ã�����coords_append�����ĸ����㷨��internal_coords���и��£�
*		8.������model������Աcoords         ������Ԫ�ز��뵽��model������Աcoords          �����ʵ�λ�ã�����coords_append�����ĸ����㷨��coords         ���и��£�
*		9.������model������Աminus_forces   ������Ԫ�ز��뵽��model������Աminus_forces    �����ʵ�λ�ã�����coords_append�����ĸ����㷨��minus_forces   ���и��£�
*		
*		10.������model������Աligands     ������Ԫ����ӵ���model������Աligands     ����β����������model.append�����ĸ����㷨�Ա�model������Աligands     ������Ԫ�ؽ��и���
*		11.������model������Աflex        ������Ԫ����ӵ���model������Աflex        ����β����������model.append�����ĸ����㷨�Ա�model������Աflex        ������Ԫ�ؽ��и���
*		12.������model������Աflex_context������Ԫ����ӵ���model������Աflex_context����β����������model.append�����ĸ����㷨�Ա�model������Աflex_context������Ԫ�ؽ��и���
*       
*		13.������model������Աgrid_atoms������Ԫ����ӵ���model������Աgrid_atoms����β����������model.append�����ĸ����㷨�Ա�model������Աgrid_atoms������Ԫ�ؽ��и���
*		14.������model������Աatoms������Ԫ�ز��뵽��model������Աatoms�����ʵ�λ�ã�����coords_append�����ĸ����㷨��atoms���и��£�
*		15.����������m_num_movable_atoms׷�ӵ��������m_num_movable_atoms
*/


void model::append(const model& m) {
	VINA_CHECK(atom_typing_used() == m.atom_typing_used());

	appender t(*this, m);

	t.append(other_pairs, m.other_pairs);

	VINA_FOR_IN(i, atoms)
		VINA_FOR_IN(j, m.atoms) {
			if(i >= m_num_movable_atoms && j >= m.m_num_movable_atoms) continue; // no need for inflex-inflex interactions

			const atom& a =   atoms[i];
			const atom& b = m.atoms[j];

			sz t1 = a.get(atom_typing_used());
			sz t2 = b.get(atom_typing_used());
			sz n = num_atom_types(atom_typing_used());

			if(t1 < n && t2 < n) {
				t.is_a =  true;
				sz new_i = t(i);
				t.is_a = false;
				sz new_j = t(j);
				sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
				other_pairs.push_back(interacting_pair(type_pair_index, new_i, new_j));
			}
		}

	VINA_CHECK(  minus_forces.size() ==   coords.size());
	VINA_CHECK(m.minus_forces.size() == m.coords.size());

	t.coords_append(internal_coords, m.internal_coords);
	t.coords_append(         coords, m         .coords);
	t.coords_append(   minus_forces, m   .minus_forces); // for now, minus_forces.size() == coords.size() (includes inflex)

	t.append(ligands,         m.ligands);
	t.append(flex,            m.flex);
	t.append(flex_context,    m.flex_context);

	t       .append(grid_atoms, m.grid_atoms);
	t.coords_append(     atoms, m     .atoms);

	m_num_movable_atoms += m.m_num_movable_atoms;
}

///////////////////  end  MODEL::APPEND /////////////////////////


/////////////////// begin MODEL::INITIALIZE /////////////////////////


/*
* ����model�ĳ�Ա����sz_to_atom_index
* input:unsigned int����i
* ouput:atom_index��
* ���ã��ж�i�Ƿ�С��model��Աgrid_atoms�Ĵ�С�����򷵻�һ��atom_index�࣬���ҳ�ʼ������ĳ�ԱiΪi��in_gridΪtrue
*                                               ���򷵻�һ��atom_index�࣬���ҳ�ʼ������ĳ�ԱiΪ(i-model��Աgrid_atoms�Ĵ�С)��in_gridΪfalse
*/
atom_index model::sz_to_atom_index(sz i) const {
	if(i < grid_atoms.size()) return atom_index(i                    ,  true);
	else                      return atom_index(i - grid_atoms.size(), false);
}

/*
* ����model�ĳ�Ա����distance_type_between
* input:distance_type_matrix�ೣ������mobility��atom_index�ೣ������i��atom_index�ೣ������j
* output:ö������distance_type 
*        ������˳��ѡ�����
*        1.i�ĳ�Աin_gridΪ�沢��j�ĳ�Աin_gridΪ�淵��DISTANCE_FIXED
*        2.ֻ��i�ĳ�Աin_gridΪ�棬�����j�ĳ�Աi�Ƿ�С��m_num_movable_atoms���Ƿ���DISTANCE_VARIABLE���񷵻�DISTANCE_FIXED
*		 3.ֻ��j�ĳ�Աin_gridΪ�棬�����i�ĳ�Աi�Ƿ�С��m_num_movable_atoms���Ƿ���DISTANCE_VARIABLE���񷵻�DISTANCE_FIXED
*        4.���i�ĳ�ԱiС��model�ĳ�Աatoms�����Ĵ�С����j�ĳ�ԱiС��model�ĳ�Աatoms�����Ĵ�С��һ��Ϊ�٣��򱨴���ֹ����ִ��
*        5.���i�ĳ�Աi����j�ĳ�Աi������DISTANCE_FIXED
*        6.���i�ĳ�ԱiС��j�ĳ�Աi����a,b�ֱ�����mobility�ĵ��������������mobility��Աm_data�ĵڣ�a + b*(b-1)/2��Ԫ��
*        7.���i�ĳ�Աi����j�ĳ�Աi����b,a�ֱ�����mobility�ĵ��������������mobility��Աm_data�ĵڣ�b + a*(a-1)/2��Ԫ��
*/
distance_type model::distance_type_between(const distance_type_matrix& mobility, const atom_index& i, const atom_index& j) const {
	if(i.in_grid && j.in_grid) return DISTANCE_FIXED;
	if(i.in_grid) return (j.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
	if(j.in_grid) return (i.i < m_num_movable_atoms) ? DISTANCE_VARIABLE : DISTANCE_FIXED;
	assert(!i.in_grid);
	assert(!j.in_grid);
	assert(i.i < atoms.size());
	assert(j.i < atoms.size());
	sz a = i.i;
	sz b = j.i;
	if(a == b) return DISTANCE_FIXED;
	return (a < b) ? mobility(a, b) : mobility(b, a);//////////////////////////////////////��֪��ʲô��˼
}

/*
* ����model�ĳ�Ա����atom_coords
* input:atom_index������i
* ouput:vec�ೣ������
* ���ã��ж�i�ĳ�Աin_grid�Ƿ�Ϊ�棬���򷵻�model�ĳ�Աgrid_atoms�����ĵڣ�i�ĳ�Աi����ֵ����Ԫ�صĳ�Աcoordsʸ��
*                                   ���򷵻�model�ĳ�Աcoords�����ĵڣ�i�ĳ�Աi����ֵ����Ԫ��
*/
const vec& model::atom_coords(const atom_index& i) const {
	return i.in_grid ? grid_atoms[i.i].coords : coords[i.i];
}

/*
* ����model�ĳ�Ա����distance_sqr_between
* input:atom_index�ೣ������a,b
* ouput:double������
* ���ã���a,b�ֱ���atom_coords����������������άʸ��������������vec�����data����Ԫ�ز��ƽ����                                   
*/
fl model::distance_sqr_between(const atom_index& a, const atom_index& b) const {
	return vec_distance_sqr(atom_coords(a), atom_coords(b));
}


/*
* �ṹ��bond_less
* ������ò���������
* input:bond�ೣ������a,b
* ouput:����ֵ
* ���ã��Ƚ�a�ĳ�Աconnected_atom_index�ĳ�Աi��b�ĳ�Աconnected_atom_index�ĳ�Աi,ǰ�ߴ󷵻�false�����ߴ󷵻�True
*/
struct bond_less { // FIXME rm!?
	bool operator()(const bond& a, const bond& b) const {
		return a.connected_atom_index.i < b.connected_atom_index.i;
	}
};

/*
* ����model�ĳ�Ա����atom_exists_between
* input:distance_type_matrix�ೣ������mobility��atom_index�ೣ������a��atom_index�ೣ������b,unsigend int����vector��������relevant_atoms
* output:boolֵ��Ѱ���������������ԭ���Ƿ����
* ���ã�1.����a��b��in_grid�Ƿ�Ϊ�棬�� r2 = (grid_atoms[a.i].coords[0]-grid_atoms[b.i].coords[0])^2+(grid_atoms[a.i].coords[1]-grid_atoms[b.i].coords[1])^2+(grid_atoms[a.i].coords[2]-grid_atoms[b.i].coords[2])^2------------------------a.in_grid = TRUE  , b.in_grid = TRUE
*									    
*									    r2 = (coords[a.i][0]-grid_atoms[b.i].coords[0])^2+(coords[a.i][1]-grid_atoms[b.i].coords[1])^2+(coords[a.i][2]-grid_atoms[b.i].coords[2])^2---------------------------------------------------------a.in_grid = FALSE , b.in_grid = TRUE 
* 
*										r2 = (grid_atoms[a.i].coords[0]-coords[b.i][0])^2+(grid_atoms[a.i].coords[1]-coords[b.i][1])^2+(grid_atoms[a.i].coords[2]-coords[b.i][2])^2---------------------------------------------------------a.in_grid = TRUE  , b.in_grid = FALSE
* 
*										r2 = (coords[a.i][0]-coords[b.i][0])^2+(coords[a.i][1]-coords[b.i][1])^2+(coords[a.i][2]-coords[b.i][2])^2------------------------------------------------------------------------------------------a.in_grid = FALSE , b.in_grid = FALSE
*		2.����relevant_atoms�ڵ�Ԫ��
*				��relevant_atoms�ĵ�relevant_atoms_iԪ�ظ�ֵ��i
*				�ж�i�Ƿ�С��model��Աgrid_atoms�Ĵ�С�����򷵻�һ��atom_index��c�����ҳ�ʼ������ĳ�ԱiΪi��in_gridΪtrue;���򷵻�һ��atom_index��c�����ҳ�ʼ������ĳ�ԱiΪ(i-model��Աgrid_atoms�Ĵ�С)��in_gridΪfalse
*		 		��� a����c����b����c ֱ�ӽ�����һ��ѭ��                                        
*				��mobility, a, c����distance_type_between�������صõ�distance_typeö������ac
*				��mobility, b, c����distance_type_between�������صõ�distance_typeö������bc
*				���ac������DISTANCE_VARIABLE������bc������DISTANCE_VARIABLE������distance_sqr_between(a, c)С��r2������distance_sqr_between(b, c)С��r2�����棬�˳�����
*				����distance_sqr_between(a, c)����a��c��in_grid�Ƿ�Ϊ�棬�� r2 = (grid_atoms[a.i].coords[0]-grid_atoms[c.i].coords[0])^2+(grid_atoms[a.i].coords[1]-grid_atoms[c.i].coords[1])^2+(grid_atoms[a.i].coords[2]-grid_atoms[c.i].coords[2])^2------------------------a.in_grid = TRUE  , c.in_grid = TRUE
* 											  						    																																																							
* 											  						        r2 = (coords[a.i][0]-grid_atoms[c.i].coords[0])^2+(coords[a.i][1]-grid_atoms[c.i].coords[1])^2+(coords[a.i][2]-grid_atoms[c.i].coords[2])^2---------------------------------------------------------a.in_grid = FALSE , c.in_grid = TRUE 
* 																		    																																																						
* 											  							    r2 = (grid_atoms[a.i].coords[0]-coords[c.i][0])^2+(grid_atoms[a.i].coords[1]-coords[c.i][1])^2+(grid_atoms[a.i].coords[2]-coords[c.i][2])^2---------------------------------------------------------a.in_grid = TRUE  , c.in_grid = FALSE
* 																		    																																																						
* 											  							    r2 = (coords[a.i][0]-coords[c.i][0])^2+(coords[a.i][1]-coords[c.i][1])^2+(coords[a.i][2]-coords[c.i][2])^2------------------------------------------------------------------------------------------a.in_grid = FALSE , c.in_grid = FALSE
* 
*				����distance_sqr_between(b, c)����a��c��in_grid�Ƿ�Ϊ�棬�� r2 = (grid_atoms[b.i].coords[0]-grid_atoms[c.i].coords[0])^2+(grid_atoms[b.i].coords[1]-grid_atoms[c.i].coords[1])^2+(grid_atoms[b.i].coords[2]-grid_atoms[c.i].coords[2])^2------------------------b.in_grid = TRUE  , c.in_grid = TRUE
* 											  						    																																																							
* 											  						        r2 = (coords[b.i][0]-grid_atoms[c.i].coords[0])^2+(coords[b.i][1]-grid_atoms[c.i].coords[1])^2+(coords[b.i][2]-grid_atoms[c.i].coords[2])^2---------------------------------------------------------b.in_grid = FALSE , c.in_grid = TRUE 
* 																		    																																																						
* 											  							    r2 = (grid_atoms[b.i].coords[0]-coords[c.i][0])^2+(grid_atoms[b.i].coords[1]-coords[c.i][1])^2+(grid_atoms[b.i].coords[2]-coords[c.i][2])^2---------------------------------------------------------b.in_grid = TRUE  , c.in_grid = FALSE
* 																		    																																																						
* 											  							    r2 = (coords[b.i][0]-coords[c.i][0])^2+(coords[a.i][1]-coords[b.i][1])^2+(coords[b.i][2]-coords[c.i][2])^2------------------------------------------------------------------------------------------b.in_grid = FALSE , c.in_grid = FALSE
*		3.���ϱ����Ҳ�����ص�ԭ�ӣ����ؼ�
* 
* 
* 
* 
* 
*        
*/
bool model::atom_exists_between(const distance_type_matrix& mobility, const atom_index& a, const atom_index& b, const szv& relevant_atoms) const { // there is an atom closer to both a and b then they are to each other and immobile relative to them
	fl r2 = distance_sqr_between(a, b);
	VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
		sz i = relevant_atoms[relevant_atoms_i];
		atom_index c = sz_to_atom_index(i);
		if(a == c || b == c) continue;
		distance_type ac = distance_type_between(mobility, a, c);
		distance_type bc = distance_type_between(mobility, b, c);
		if(ac != DISTANCE_VARIABLE &&
		   bc != DISTANCE_VARIABLE &&
		   distance_sqr_between(a, c) < r2 &&
		   distance_sqr_between(b, c) < r2)
			return true;
	}
	return false;
}

/*
* �ṹ��beads
* ��Ա���£�
*			1.double����radius_sqr
*			2.vector��data��Ԫ��Ϊpair�࣬pair���first���ݳ�ԱΪvec�࣬second���ݳ�ԱΪunsignd int����vector
*			3.���캯��beads��input:unsigned int����reserve_size,double��radius_sqr
*							 ���ã���ʼ��beads���radius_sqrΪ����radius_sqr��data��������������reserve_size��Ԫ�ص��ڴ�ռ�
*			4.����add��input:unsigned int����index,vec�ೣ������coords
* 					   ���ã�1.����vector��data
*									���(coords)��(data�ĵ�i��Ԫ�ص�first���ݳ�Ա)�����data����Ԫ�ز��ƽ��<radius_sqr
*									��data�ĵ�i��Ԫ�ص�second���ݳ�Աβ�����Ԫ��index
*									�˳�����
*							 2.��һ��û���ҶԺ���Ԫ�ص�����£�����һ��pair��tmp,��first���ݳ�ԱΪvec�࣬second���ݳ�ԱΪunsignd int����vector
*							   ��tmp��first���ݳ�Ա��ֵΪcoords��second���ݳ�Աβ�����index����tmp��ӵ�dataβ��
*/
struct beads {
	fl radius_sqr;
	std::vector<std::pair<vec, szv> > data;
	beads(sz reserve_size, fl radius_sqr_) : radius_sqr(radius_sqr_) { data.reserve(reserve_size); }
	void add(sz index, const vec& coords) {
		VINA_FOR_IN(i, data) {
			if(vec_distance_sqr(coords, data[i].first) < radius_sqr) {
				data[i].second.push_back(index);
				return;
			}
		}
		// not found
		std::pair<vec, szv> tmp;
		tmp.first = coords;
		tmp.second.push_back(index);
		data.push_back(tmp);
	}
};
/*
* ����model�ĳ�Ա����assign_bonds
* input:distance_type_matrix�ೣ������mobility
* ouput:��
* ���ã�1.����һ��double����ֵbond_length_allowance_factor = 1.1 ճ�᳤������ϵ��
*		  n = grid_atoms �� atoms ������С֮��
*		  bead_radius = 1.5 
*		2.����һ��beads�����beads_instance����ʼ��beads_instance��radius_sqr = 225 ,data��Ԥ���ռ�Ϊn
*		3.i��0ѭ����n
*				����һ��atom_index��i_atom_index����0 <= i < grid_atoms.size()ʱ����ʼ������ĳ�ԱiΪi��in_gridΪtrue����i >= grid_atoms.size()ʱ����ʼ������ĳ�ԱiΪ(i-model��Աgrid_atoms�Ĵ�С)��in_gridΪfalse
*				��0 <= i < grid_atoms.size()ʱ����i��grid_atoms[i_atom_index.i].coords����beads_instance��Ա����add����i >= grid_atoms.size()ʱ����i��coords[i_atom_index.i]����beads_instance��Ա����add
* 
*				��0 <= i < grid_atoms.size()ʱ����beads_instance.data[i](��i���ڱ���data����ǰ��i��0ѭ����n)��Ѱ�ҵ�һ������[data[i].first[0] - grid_atoms[i_atom_index.i].coords[0]]^2+[data[i].first[1] - grid_atoms[i_atom_index.i].coords[1]]^2+[data[i].first[2] - grid_atoms[i_atom_index.i].coords[2]]^2<225��beads_instance.data��Ԫ�ز���data[i].second��β������i
*				���data����������������Ԫ�أ�����һ��pair��tmp��first���ݳ�Ա����grid_atoms[i_atom_index.i].coords��second���ݳ�Աβ�����i����tmp����beads_instance.dataβ��
* 
*				��i >= grid_atoms.size()ʱ����beads_instance.data[i](��i���ڱ���data����ǰ��i��0ѭ����n)��Ѱ�ҵ�һ������[data[i].first[0] - coords[i_atom_index.i][0]]^2+[data[i].first[1] - coords[i_atom_index.i][1]]^2+[data[i].first[2] - coords[i_atom_index.i][2]]^2<225��beads_instance.data��Ԫ�ز���data[i].second��β������i
*				���data����������������Ԫ�أ�����һ��pair��tmp��first���ݳ�Ա����coords[i_atom_index.i]��second���ݳ�Աβ�����i����tmp����beads_instance.dataβ��
* 
*		4.i��0ѭ����n
*				����һ��atom_index��i_atom_index����0 <= i < grid_atoms.size()ʱ����ʼ������ĳ�ԱiΪi��in_gridΪtrue����i >= grid_atoms.size()ʱ����ʼ������ĳ�ԱiΪ(i-model��Աgrid_atoms�Ĵ�С)��in_gridΪfalse
*				��0 <= i < grid_atoms.size()ʱ����grid_atoms[i_atom_index.i].coords�󶨵�i_atom_coords����i >= grid_atoms.size()ʱ����coords[i_atom_index.i]�󶨵�i_atom_coords
*				��0 <= i < grid_atoms.size()ʱ����grid_atoms[i_atom_index.i]�󶨵�i_atom����i >= grid_atoms.size()ʱ����atoms[i_atom_index.i]�󶨵�i_atom
*				
*				����һ��double��������ֵmax_covalent_r = 1.74----------->atom_kind_data[i].covalent_radius���ֵ     ���۰뾶
*				����һ��double����ֵi_atom_covalent_radius = max_covalent_r = 1.74
*				���i_atom��adС��20��i_atom_covalent_radius��ֵΪatom_kind_data��i_atom.ad����covalent_radius
* 
*				����һ��unsigned int����vector -->relevant_atoms
*				����һ��double��������ֵbead_cutoff_sqr = (15 + 1.1*(i_atom_covalent_radius+1.74))^2
* 
*				b��0ѭ����beads_instance.data�Ĵ�С-1������Ԫ�أ�
*						����beads_instance.data[b].first[0] - i_atom_coords[0])^2 + ��beads_instance.data[b].first[1] - i_atom_coords[1])^2+��beads_instance.data[b].first[2] - i_atom_coords[2])^2> bead_cutoff_sqrʱ��ֱ�Ӽ�����һ��b����
*						��beads_instance.data[b].second�󶨵�bead_elements��unsigned int����vector��
*						b��bead_elements_iѭ����bead_elements�Ĵ�С-1������Ԫ�أ�
*								����j = bead_elements�ĵ�bead_elements_i��Ԫ��
*								����һ��atom_index��j_atom_index����0 <= j < grid_atoms.size()ʱ����ʼ������ĳ�ԱiΪj��in_gridΪtrue����i >= grid_atoms.size()ʱ����ʼ������ĳ�ԱiΪ(j-model��Աgrid_atoms�Ĵ�С)��in_gridΪfalse
*								��0 <= j < grid_atoms.size()ʱ����grid_atoms[j_atom_index.i]�󶨵�j_atom����j >= grid_atoms.size()ʱ����atoms[j_atom_index.i]�󶨵�j_atom
*								����i_atom��j_atom�Ĺ��۰뾶��bond_length
*								��mobility, i_atom_index, j_atom_index����distance_type_between��������ȡi_atom_index��j_atom_index�ľ�������dt
*								���dt != DISTANCE_VARIABLE����i != j��
*										��0 <= i < grid_atoms.size()����0 <= j < grid_atoms.size()ʱ����r2 = [grid_atoms[i_atom_index.i].coords[0]-grid_atoms[j_atom_index.i].coords[0]]^2+[grid_atoms[i_atom_index.i].coords[1]-grid_atoms[j_atom_index.i].coords[1]]^2+[grid_atoms[i_atom_index.i].coords[2]-grid_atoms[j_atom_index.i].coords[2]]^2
*										��0 <= i < grid_atoms.size()����j > grid_atoms.size()ʱ����r2 = [grid_atoms[i_atom_index.i].coords[0]-coords[j_atom_index.i][0]]^2+[grid_atoms[i_atom_index.i].coords[1]-coords[j_atom_index.i][1]]^2+[grid_atoms[i_atom_index.i].coords[2]-coords[j_atom_index.i][2]]^2
*										��i > grid_atoms.size()����0 <= j < grid_atoms.size()ʱ����r2 = [coords[i_atom_index.i][0]-grid_atoms[j_atom_index.i].coords[0]]^2+[coords[i_atom_index.i][1]-grid_atoms[j_atom_index.i].coords[1]]^2+[coords[i_atom_index.i][2]-grid_atoms[j_atom_index.i].coords[2]]^2
*										��i > grid_atoms.size()����j > grid_atoms.size()ʱ����r2 = [coords[i_atom_index.i][0]-coords[j_atom_index.i][0]]^2+[coords[i_atom_index.i][1]-coords[j_atom_index.i][1]]^2+[coords[i_atom_index.i][2]-coords[j_atom_index.i][2]]^2
*										���r2 < (1.1*(i_atom_covalent_radius+1.74))^2
*										relevant_atomsβ�����j
*				
*				relevant_atoms_i��0ѭ����relevant_atoms�Ĵ�С-1������Ԫ�أ�
*						����j����relevant_atoms��relevant_atoms_i��Ԫ��
*						���j <= i,ֱ�Ӽ�����һ��relevant_atoms_i����
*						����һ��atom_index��j_atom_index����0 <= j < grid_atoms.size()ʱ����ʼ������ĳ�ԱiΪj��in_gridΪtrue����i >= grid_atoms.size()ʱ����ʼ������ĳ�ԱiΪ(j-model��Աgrid_atoms�Ĵ�С)��in_gridΪfalse
* 						��0 <= j < grid_atoms.size()ʱ����grid_atoms[j_atom_index.i]�󶨵�j_atom����j >= grid_atoms.size()ʱ����atoms[j_atom_index.i]�󶨵�j_atom
*						����i_atom��j_atom�Ĺ��۰뾶��bond_length
*						��mobility, i_atom_index, j_atom_index����distance_type_between��������ȡi_atom_index��j_atom_index�ľ�������dt
*						��0 <= i < grid_atoms.size()����0 <= j < grid_atoms.size()ʱ����r2 = [grid_atoms[i_atom_index.i].coords[0]-grid_atoms[j_atom_index.i].coords[0]]^2+[grid_atoms[i_atom_index.i].coords[1]-grid_atoms[j_atom_index.i].coords[1]]^2+[grid_atoms[i_atom_index.i].coords[2]-grid_atoms[j_atom_index.i].coords[2]]^2								
* 						��0 <= i < grid_atoms.size()����j > grid_atoms.size()ʱ����r2 = [grid_atoms[i_atom_index.i].coords[0]-coords[j_atom_index.i][0]]^2+[grid_atoms[i_atom_index.i].coords[1]-coords[j_atom_index.i][1]]^2+[grid_atoms[i_atom_index.i].coords[2]-coords[j_atom_index.i][2]]^2
* 						��i > grid_atoms.size()����0 <= j < grid_atoms.size()ʱ����r2 = [coords[i_atom_index.i][0]-grid_atoms[j_atom_index.i].coords[0]]^2+[coords[i_atom_index.i][1]-grid_atoms[j_atom_index.i].coords[1]]^2+[coords[i_atom_index.i][2]-grid_atoms[j_atom_index.i].coords[2]]^2
* 						��i > grid_atoms.size()����j > grid_atoms.size()ʱ����r2 = [coords[i_atom_index.i][0]-coords[j_atom_index.i][0]]^2+[coords[i_atom_index.i][1]-coords[j_atom_index.i][1]]^2+[coords[i_atom_index.i][2]-coords[j_atom_index.i][2]]^2
*						���r2 < (1.1*bond_length)^2����mobility, i_atom_index, j_atom_index, relevant_atoms ����atom_exists_between��Ϊ�棨relevant_atoms�в�������i_atom_index��j_atom_index����һ�������ϵ��ԭ�ӣ�
*								���dt ���� DISTANCE_ROTOR����rotatable����1������Ϊ0
*								length ���� r2�Ŀ�����
*								��i_atom.bondsβ�����һ��bondԪ��,����j_atom_index, length, rotatable������г�ʼ��
*								��j_atom.bondsβ�����һ��bondԪ��,����j_atom_index, length, rotatable������г�ʼ��
*/
void model::assign_bonds(const distance_type_matrix& mobility) { // assign bonds based on relative mobility, distance and covalent length
	const fl bond_length_allowance_factor = 1.1;
	sz n = grid_atoms.size() + atoms.size();

	// construct beads
	const fl bead_radius = 15;
	beads beads_instance(n, sqr(bead_radius));
	VINA_FOR(i, n) {
		atom_index i_atom_index = sz_to_atom_index(i);
		beads_instance.add(i, atom_coords(i_atom_index));
	}
	// assign bonds
	VINA_FOR(i, n) {
		atom_index i_atom_index = sz_to_atom_index(i);
		const vec& i_atom_coords = atom_coords(i_atom_index);
		atom& i_atom = get_atom(i_atom_index);

		const fl max_covalent_r = max_covalent_radius(); // FIXME mv to atom_constants
		fl i_atom_covalent_radius = max_covalent_r;
		if(i_atom.ad < AD_TYPE_SIZE)
			i_atom_covalent_radius = ad_type_property(i_atom.ad).covalent_radius;

		//find relevant atoms
		szv relevant_atoms;
		const fl bead_cutoff_sqr = sqr(bead_radius + bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r));
		VINA_FOR_IN(b, beads_instance.data) {
			if(vec_distance_sqr(beads_instance.data[b].first, i_atom_coords) > bead_cutoff_sqr) continue;
			const szv& bead_elements = beads_instance.data[b].second;
			VINA_FOR_IN(bead_elements_i, bead_elements) {
				sz j = bead_elements[bead_elements_i];
				atom_index j_atom_index = sz_to_atom_index(j);
				atom& j_atom = get_atom(j_atom_index);
				const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
				distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
				if(dt != DISTANCE_VARIABLE && i != j) {
					fl r2 = distance_sqr_between(i_atom_index, j_atom_index);
					//if(r2 < sqr(bond_length_allowance_factor * bond_length))
					if(r2 < sqr(bond_length_allowance_factor * (i_atom_covalent_radius + max_covalent_r)))
						relevant_atoms.push_back(j);
				}
			}
		}
		// find bonded atoms
		VINA_FOR_IN(relevant_atoms_i, relevant_atoms) {
			sz j = relevant_atoms[relevant_atoms_i];
			if(j <= i) continue; // already considered
			atom_index j_atom_index = sz_to_atom_index(j);
			atom& j_atom = get_atom(j_atom_index);
			const fl bond_length = i_atom.optimal_covalent_bond_length(j_atom);
			distance_type dt = distance_type_between(mobility, i_atom_index, j_atom_index);
			fl r2 = distance_sqr_between(i_atom_index, j_atom_index);///////////////////////////////////////////////////////
			if(r2 < sqr(bond_length_allowance_factor * bond_length) && !atom_exists_between(mobility, i_atom_index, j_atom_index, relevant_atoms)) {
				bool rotatable = (dt == DISTANCE_ROTOR);
				fl length = std::sqrt(r2);
				i_atom.bonds.push_back(bond(j_atom_index, length, rotatable));
				j_atom.bonds.push_back(bond(i_atom_index, length, rotatable));
			}

		}
	}
}


/*
* ����model�ĳ�Ա����bonded_to_HD
* input:atom�ೣ������a
* ouput:boolֵ
* ���ã�1.����a�ĳ�Աbonds��Ԫ��
*       2.����һ����ʱbond�����������b����a�ĳ�Աbonds�ĵ�i��Ԫ�ذ󶨵�b��
*       3.��b�ĳ�Աconnected_atom_index����get_atom����������һ��atom������жϸö���ĳ�Աad�Ƿ����12����������˳���������true���������ѭ��
*       4.ѭ������������false
*/
bool model::bonded_to_HD(const atom& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).ad == AD_TYPE_HD) 
			return true;
	}
	return false;
}


/*
* ����model�ĳ�Ա����bonded_to_heteroatom
* input:atom�ೣ������a
* ouput:boolֵ
* ���ã�1.����a�ĳ�Աbonds��Ԫ��
*       2.����һ����ʱbond�����������b����a�ĳ�Աbonds�ĵ�i��Ԫ�ذ󶨵�a��
*       3.��b�ĳ�Աconnected_atom_index����get_atom����������һ��atom������жϸö���ĳ�Աis_heteroatom��������ֵ�Ƿ�Ϊ�棬��������˳���������true���������ѭ��
*       4.ѭ������������false
*/
bool model::bonded_to_heteroatom(const atom& a) const {
	VINA_FOR_IN(i, a.bonds) {
		const bond& b = a.bonds[i];
		if(get_atom(b.connected_atom_index).is_heteroatom())
			return true;
	}
	return false;
}

/*
* ����model�ĳ�Ա����assign_types
* input:��
* ouput:��
* ���ã�i��0ѭ������model��Աgrid_atoms�Ĵ�С��atoms�Ĵ�С֮��-1��
*				��i���뵽model��Ա����sz_to_atom_index�����ص�ֵ��ֵ��ai�ϣ���0<=i<=(grid_atoms�Ĵ�С-1)ʱ��aiΪatom_index(i,true)������Ϊatom_index(i - grid_atoms.size(), false)
*				��ai����model��Ա����get_atom��aiΪatom_index(i,true)ʱ�����model��Աgrid_atoms�ĵڣ�ai�ĳ�Աi����Ԫ�أ�aiΪatom_index(i - grid_atoms.size(), false)ʱ�����model��Աatoms�ĵ�i��ai�ĳ�Աi����Ԫ�أ����ص�ֵ�󶨵�a��
*				��a�ĳ�Աel��ֵ����a�ĳ�Աad==20��a�ĳ�Աxs==16 ��el=10������el=ad����ad_type_to_el_type������ķ���ֵ
*				��a�ĳ�Աxs�󶨵�x��
*				����һ��boolֵacceptor�����a�ĳ�Աad����9����10��acceptorΪtrue,����Ϊfalse
*				����һ��boolֵdonor_NorO�����a�ĳ�Աel����10���߽�a����model��Ա����bonded_to_HD����Ϊ�棬donor_NorOΪtrue,����Ϊfalse
*				����a�ĳ�Աel��СΪx��ֵ
*							��el =  EL_TYPE_H    =  0ʱ����ִ���κβ���
* 							��el =  EL_TYPE_C    =  1ʱ����a����model��Ա����bonded_to_heteroatom��������x��ֵΪ1������x��ֵΪ0
* 							��el =  EL_TYPE_N    =  2ʱ�����acceptor��donor_NorO��Ϊ�棬x��ֵΪ5��ֻ��acceptorΪ��x��ֵΪ4��ֻ��donor_NorOΪ��x��ֵΪ3����Ϊ��x��ֵΪ2
* 							��el =  EL_TYPE_O    =  3ʱ�����acceptor��donor_NorO��Ϊ�棬x��ֵΪ9��ֻ��acceptorΪ��x��ֵΪ8��ֻ��donor_NorOΪ��x��ֵΪ7����Ϊ��x��ֵΪ6
* 							��el =  EL_TYPE_S    =  4ʱ��x��ֵΪXS_TYPE_S_P   = 10
* 							��el =  EL_TYPE_P    =  5ʱ��x��ֵΪXS_TYPE_P_P   = 11
* 							��el =  EL_TYPE_F    =  6ʱ��x��ֵΪXS_TYPE_F_H   = 12
* 							��el =  EL_TYPE_Cl   =  7ʱ��x��ֵΪXS_TYPE_Cl_H  = 13
* 							��el =  EL_TYPE_Br   =  8ʱ��x��ֵΪXS_TYPE_Br_H  = 14
* 							��el =  EL_TYPE_I    =  9ʱ��x��ֵΪXS_TYPE_I_H   = 15
* 							��el =  EL_TYPE_Met  = 10ʱ��x��ֵΪXS_TYPE_Met_D = 16
* 							��el =  EL_TYPE_SIZE = 11ʱ����ִ���κβ���
*							�����ǣ�������ֹ����ִ��
*
*/
void model::assign_types() {
	VINA_FOR(i, grid_atoms.size() + atoms.size()) {
		const atom_index ai = sz_to_atom_index(i);
		atom& a = get_atom(ai);
		a.assign_el();
		sz& x = a.xs;

		bool acceptor   = (a.ad == AD_TYPE_OA || a.ad == AD_TYPE_NA); // X-Score forumaltion apparently ignores SA
		bool donor_NorO = (a.el == EL_TYPE_Met || bonded_to_HD(a));

		switch(a.el) {
			case EL_TYPE_H    : break;
			case EL_TYPE_C    : x = bonded_to_heteroatom(a) ? XS_TYPE_C_P : XS_TYPE_C_H; break;
			case EL_TYPE_N    : x = (acceptor && donor_NorO) ? XS_TYPE_N_DA : (acceptor ? XS_TYPE_N_A : (donor_NorO ? XS_TYPE_N_D : XS_TYPE_N_P)); break;
			case EL_TYPE_O    : x = (acceptor && donor_NorO) ? XS_TYPE_O_DA : (acceptor ? XS_TYPE_O_A : (donor_NorO ? XS_TYPE_O_D : XS_TYPE_O_P)); break;
			case EL_TYPE_S    : x = XS_TYPE_S_P; break;
			case EL_TYPE_P    : x = XS_TYPE_P_P; break;
			case EL_TYPE_F    : x = XS_TYPE_F_H; break;
			case EL_TYPE_Cl   : x = XS_TYPE_Cl_H; break;
			case EL_TYPE_Br   : x = XS_TYPE_Br_H; break;
			case EL_TYPE_I    : x = XS_TYPE_I_H; break;
			case EL_TYPE_Met  : x = XS_TYPE_Met_D; break;
			case EL_TYPE_SIZE : break;
			default: VINA_CHECK(false);
		}
	}
}

/*
* ����model�ĳ�Ա����find_ligand
* input:unsigned int����a
* ouput:unsigned int����ֵ
* ���ã�1.����model�ĳ�Աligands����
*       2.�ҵ���һ��a��ligands����Ԫ�صĳ�Աbegin��end��Χ֮�ڵ�Ԫ�أ��������±�
*       3.û�з����������򷵻�ligands�����Ĵ�С
*/
sz model::find_ligand(sz a) const {
	VINA_FOR_IN(i, ligands)
		if(a >= ligands[i].begin && a < ligands[i].end)
			return i;
	return ligands.size();
}

/*
* ����model�ĳ�Ա����bonded_to
* input:unsigned int����a��n��unsigned int��������out
* ouput:unsigned int����ֵ
* ���ã�1.�ж�a�Ƿ�������out�ڣ��������������������ִ���������
*       2.��a���뵽����outβ��
*       3.�ж�n�Ƿ����0���������ִ��������䣬����������
*		4.����atoms������a��Ԫ�صĳ�Աbonds����
*       5.��atoms������a��Ԫ�صĳ�Աbonds������i��Ԫ�ذ󶨵�b��
*       6.�ж�b�ĳ�Աconnected_atom_index��ĳ�Աin_grid�Ƿ�Ϊ�棬�������������������ִ���������
*       7.��b�ĳ�Աconnected_atom_index��ĳ�Աi��n-1,out����bonded_to
* Ŀ�ģ����������������atmos���±갴����˳�����δ���outβ���������ͼ��
*/

void model::bonded_to(sz a, sz n, szv& out) const {
	if(!has(out, a)) { // not found
		out.push_back(a);
		if(n > 0) 
			VINA_FOR_IN(i, atoms[a].bonds) {
				const bond& b = atoms[a].bonds[i];
				if(!b.connected_atom_index.in_grid)
					bonded_to(b.connected_atom_index.i, n-1, out);
			}
	}
}

/*
* ����model�ĳ�Ա����bonded_to
* input:unsigned int����a��n
* ouput:unsigned int����vector
* ���ã�1.����һ����ʱunsigned int����vector����tmp
*       2.����bonded_to����������汾������a,n,tmp
*       3.����tmp
*/
szv model::bonded_to(sz a, sz n) const {
	szv tmp;
	bonded_to(a, n, tmp);
	return tmp;
}

/*
* ����model�ĳ�Ա����initialize_pairs
* input:distance_type_matrix�ೣ������mobility
* ouput:��
* ���ã�����model��Աatmos�ڵ�Ԫ��
*			�ҵ���һ��i��ligands����Ԫ�صĳ�Աbegin��end��Χ֮�ڵ�Ԫ�أ��������±긳ֵ��i_lig
*			��i,3���뺯��bonded_to������atoms����i��������һ�������Ԫ�ص��±���ɵ�vector--->bonded_atoms�������ͼ��
*			j��i+1ѭ����atoms.size()-1----------->atoms.size()Ϊatmos�����Ĵ�С
*				���i��j������model��m_num_movable_atoms��model���ƶ�ԭ�ӵ���������ֱ�Ӽ�����һ��ѭ��
*				��i,j����model�У�����mobility�ĳ�Աm_data[i + j*(j-1)/2]���жϷ���ֵ�Ƿ�ΪDISTANCE_VARIABLE����j����bonded_atoms����������㣬����ִ��������䣬���������һ��ѭ��
*				��m_atom_typing_used����atmos��i��Ԫ�صĺ���get�еõ�m_atom_typing_used��Ӧ������ֵ����t1
*				��m_atom_typing_used����atmos��j��Ԫ�صĺ���get�еõ�m_atom_typing_used��Ӧ������ֵ����t2
*				��m_atom_typing_used����num_atom_types�еõ�m_atom_typing_used��Ӧ������ֵ����n
*				���t1��t2��С��n,����ִ��������䣬���������һ��ѭ��
*				��n, t1, t2����triangular_matrix_index_permissive�����Ǿ�������������t1<t2ʱ������t1 + t2*(t2+1)/2����ֵ��type_pair_index
*																				   ��t1>t2ʱ������t2 + t1*(t1+1)/2����ֵ��type_pair_index
*				����һ��interacting_pair��ip����ʼ����ʼ��type_pair_index��a��bΪtype_pair_index, i, j
*				���i_ligС��ligands��������������С���ҵ�һ��j��ligands����Ԫ�صĳ�Աbegin��end��Χ֮�ڵ�Ԫ�أ��������±����i_lig
*				��ip����ligands��i_lig��Ԫ�صĳ�Աpairsβ�����������model��other_pairsβ��
* 
*/
void model::initialize_pairs(const distance_type_matrix& mobility) {
	VINA_FOR_IN(i, atoms) {
		sz i_lig = find_ligand(i);
		szv bonded_atoms = bonded_to(i, 3);
		VINA_RANGE(j, i+1, atoms.size()) {
			if(i >= m_num_movable_atoms && j >= m_num_movable_atoms) continue; // exclude inflex-inflex
			if(mobility(i, j) == DISTANCE_VARIABLE && !has(bonded_atoms, j)) {
				sz t1 = atoms[i].get  (atom_typing_used());
				sz t2 = atoms[j].get  (atom_typing_used());
				sz n  = num_atom_types(atom_typing_used());
				if(t1 < n && t2 < n) { // exclude, say, Hydrogens
					sz type_pair_index = triangular_matrix_index_permissive(n, t1, t2);
					interacting_pair ip(type_pair_index, i, j);
					if(i_lig < ligands.size() && find_ligand(j) == i_lig)
						ligands[i_lig].pairs.push_back(ip);
					else
						other_pairs.push_back(ip);
				}
			}
		}
	}
}

/*
* ����model�ĳ�Ա����initialize
* input:distance_type_matrix�ೣ������mobility
* ouput:��
* ���ã�����model��Աligands�ڵ�Ԫ��
*			��ȡligand[i]�ĳ�Աnode��begin���������к���ĳ�Աnode��begin���ҵ���Сֵ��ֵ��ligand[i].begin
*			��ȡligand[i]�ĳ�Աnode��end��  �������к���ĳ�Աnode��end��  �ҵ����ֵ��ֵ��ligand[i].end
*		��mobility����assign_bonds��Ѱ�Ҽ���ԭ��
*		��grid_atoms��Ԫ�أ�atoms��Ԫ�ط�������
*		��mobilityinitialize_pairs����ʼ�������Ԫ�ص�pair
*/
void model::initialize(const distance_type_matrix& mobility) {
	VINA_FOR_IN(i, ligands)
		ligands[i].set_range();
	assign_bonds(mobility);
	assign_types();
	initialize_pairs(mobility);
}

///////////////////  end  MODEL::INITIALIZE /////////////////////////

/*
* ����model�ĳ�Ա����num_internal_pairs
* input:��
* ouput:unsigned int����ֵ
* ���ã�1.����һ����ʱunsigned int����ֵtmp����ʼ��Ϊ0
*       2.��model��Աligands����Ԫ�صĳ�Աpairs�Ĵ�С֮�ͣ���ֵ��tmp������
*/
sz model::num_internal_pairs() const {
	sz tmp = 0;
	VINA_FOR_IN(i, ligands)
		tmp += ligands[i].pairs.size();
	return tmp;
}

/*
* ����model�ĳ�Ա����get_movable_atom_types����ȡ���ƶ���ԭ������
* input:atom_type�ඨ���ö������t
* ouput:unsigned int��vector
* ���ã�1.����һ����ʱunsigned int��vector tmp
*       2.��atom_typing_used_ת���ɶ�Ӧ�ģ�atom_typing_used_��_TYPE_SIZE��ֵ��n
*		3.��i=1ѭ����i=m_num_movable_atoms-1
*               ��1����model��Աatoms�ĵ�i��Ԫ�ذ󶨵�a��
*               ��2����atom_typing_used_ת����a��Աatom_typing_used_��ӦСд��ʽ�������������ֵ��t
*               ��3�����tС��n����t����tmp�ڣ���t����tmpβ��
*       4.����tmp
*/
szv model::get_movable_atom_types(atom_type::t atom_typing_used_) const {
	szv tmp;
	sz n = num_atom_types(atom_typing_used_);
	VINA_FOR(i, m_num_movable_atoms) {
		const atom& a = atoms[i];
		sz t = a.get(atom_typing_used_);
		if(t < n && !has(tmp, t))
			tmp.push_back(t);
	}
	return tmp;
}

/*
* ����model�ĳ�Ա����get_size
* input:��
* ouput:conf_size��
* ���ã�1.���ligands��Ťת������������ʱtmp.ligands������
*       2.���flex��Ťת������������ʱtmp.flex������
*		3.����tmp
*/
conf_size model::get_size() const {
	conf_size tmp;
	tmp.ligands = ligands.count_torsions();
	tmp.flex    = flex   .count_torsions();
	return tmp;
}

/*
* ����model�ĳ�Ա����get_initial_conf
* input:��
* ouput:conf��
* ���ã�1.���ligands��Ťת������������ʱcs.ligands������
*       2.���flex��Ťת������������ʱcs.ligands������
*		3.����tmp
*		4.����һ��conf����tmp
*		tmp��		ligands������Ԫ�ظ���Ϊcs�е�ligands unint������Ԫ�ظ���;
*					flex������Ԫ�ظ���Ϊcs�е�flex uint������Ԫ�ظ���;
*
*					ligands[i].torsions double������Ԫ�ظ���Ϊs.ligands[ligands.size]������ֵ��Ϊ0��
*					flex[i].torsions double������Ԫ�ظ���Ϊs.flex[flex.size]������ֵ��Ϊ0��
*		5.	������Ϊnull
*		    tmp.ligands�ṹ�������У�
*		    data[0] = 0��data[1] = 0;data[2] = 0;
*		    orientation��һ��ʵ��Ϊ1���鲿��Ϊ0����Ԫ��
*		    torsions����Ԫ�ض���ֵΪ0��
*		    tmp.flex�ṹ�������У�
*		    torsions����Ԫ�ض���ֵΪ0
*		6.����ligands��Ԫ��
*			��tmp.ligands[i].rigid.position����Ϊligands[i].node.origin
*		7.����tmp
*/
conf model::get_initial_conf() const { // torsions = 0, orientations = identity, ligand positions = current
	conf_size cs = get_size();
	conf tmp(cs);
	tmp.set_to_null();
	VINA_FOR_IN(i, ligands)
		tmp.ligands[i].rigid.position = ligands[i].node.get_origin();
	return tmp;
}

/*
* ����model�ĳ�Ա����movable_atoms_box
* input:double����add_to_each_dimension��double����granularity�����ȣ�
* ouput:grid_dims����
* ���ã�1.��������vec��corner1��corner2�����Ҿ���ʼ�����ڵ�����Ԫ��Ϊ0
*       2.i��0ѭ����m_num_movable_atoms��model�ĳ�Ա��
*				��model��Աcoords�ĵ�i��Ԫ�ذ󶨵�v��
*				j��0ѭ����2
*						��i=0ʱ����v�ĵ�j��Ԫ��С��corner1�ĵ�j��Ԫ�أ���v�ĵ�j��Ԫ�ظ�ֵ��corner1�ĵ�j��Ԫ��
*						��i=0ʱ����v�ĵ�j��Ԫ�ش���corner1�ĵ�j��Ԫ�أ���v�ĵ�j��Ԫ�ظ�ֵ��corner1�ĵ�j��Ԫ��
*		 �����model��Աcoords����Ԫ�ص������������������Сֵ�����ֵ������corner1, corner2��		
*		3.corner1[0] = corner1[0]-add_to_each_dimension / 2, corner1[1] = corner1[1]-add_to_each_dimension / 2, corner1[2] = corner1[2]-add_to_each_dimension / 2
* 		  corner2[0] = corner2[0]+add_to_each_dimension / 2, corner2[1] = corner2[1]+add_to_each_dimension / 2, corner2[2] = corner2[2]+add_to_each_dimension / 2
*		4.����һ��grid_dims����gd
*		  ����corner1��corner2������ center[0] = 0.5 * (corner2[0] + corner1[0]) center[1] = 0.5 * (corner2[1] + corner1[1]) center[2] = 0.5 * (corner2[2] + corner1[2])
*		  i��0ѭ����2
*				gd[i].n = ��(corner2[i] - corner1[i]) / granularity)����ȡ����         gd[i].n---------gd[i]�ĳ�Աn
*				real_span = granularity*gd[i].n��       gd[i].begin = center[i] - real_span/2��          gd[i].end = center[i] + real_span/2             gd[i].begin---------gd[i]�ĳ�Աbegin            gd[i].end---------gd[i]�ĳ�Աend
*		5.����gd
* 
* 
*/

grid_dims model::movable_atoms_box(fl add_to_each_dimension, fl granularity) const {
	vec corner1(0, 0, 0), corner2(0, 0, 0);
	VINA_FOR(i, num_movable_atoms()) {
		const vec& v = movable_coords(i);
		VINA_FOR_IN(j, v) {
			if(i == 0 || v[j] < corner1[j]) corner1[j] = v[j];
			if(i == 0 || v[j] > corner2[j]) corner2[j] = v[j];
		}
	}
	corner1 -= add_to_each_dimension / 2;
	corner2 += add_to_each_dimension / 2;

	grid_dims gd;
	{ // always doing this now FIXME ?
		vec center; center = 0.5 * (corner2 + corner1);
		VINA_FOR_IN(i, gd) {
			gd[i].n = sz(std::ceil((corner2[i] - corner1[i]) / granularity));
			fl real_span = granularity * gd[i].n;
			gd[i].begin = center[i] - real_span/2;
			gd[i].end = gd[i].begin + real_span;
		}
	}
	return gd;
}

/*
* ����string_write_coord
* input:unsigned int����i��double����x��string������str
* ouput:��
* ���ã�1.�ȼ������i�Ƿ����0���������������룬���򱨴���ֹ����ִ��
*       2.i�Լ�1
*		3.����һ��string�������out
*       4.��x���и�ʽ���ñ�����out�У���ʽΪzzzz.zzz��������λ��ЧС��
*		5.��������out�е�string������str�У���ʼ��λ��Ϊi��������i��ȥ1��
*/
void string_write_coord(sz i, fl x, std::string& str) {
	VINA_CHECK(i > 0);
	--i;
	std::ostringstream out;
	//�趨�������ʽ,С�������6λ����,��ʾ��������С�����С��λ��
	out.setf(std::ios::fixed, std::ios::floatfield);
	out.setf(std::ios::showpoint);
	//setw(8)���������ȣ�8��ռλ�����Ҷ���
	//���ø���ֵ��С������Ϊ3
	out << std::setw(8) << std::setprecision(3) << x; 
	//���string�������string�Ĵ�С�Ƿ����8���������򱨴���ֹ����ִ��
	VINA_CHECK(out.str().size() == 8); 
	//���str�Ĵ�С�Ƿ����i+8��С�ڵ����򱨴���ֹ����ִ��
	VINA_CHECK(str.size() > i + 8);
	//��string�������out�е�Ԫ�ظ�ֵ��str�У���ʼ��λ��Ϊi
	VINA_FOR(j, 8)
		str[i+j] = out.str()[j];
}

/*
* ����coords_to_pdbqt_string
* input:vec�ೣ������coords��string�ೣ������str
* ouput:string��
* ���ã�����һ��strΪtmp������coords������Ԫ��ת���ɸ�ʽΪzzzz.zzz��������λ��ЧС�������ν�0��1��2λ�õ�Ԫ�ش���tmp��ʼλ��Ϊ30��38��46λ�ô�������tmp	
*/
std::string coords_to_pdbqt_string(const vec& coords, const std::string& str) {
	std::string tmp(str);
	string_write_coord(31, coords[0], tmp);
	string_write_coord(39, coords[1], tmp);
	string_write_coord(47, coords[2], tmp);
	return tmp;
}

/*
* ����model�ĳ�Ա����write_context
* input:parsed_line����vector��������c��ofile��out
* ouput:void
* ���ã�1.���ú���verify_bond_lengths��֤���ӳ���
*		2.����c�е�Ԫ��
*				��c�ĵ�i��Ԫ�ص�first���ݳ�Ա�󶨵�str��
*				���c�ĵ�i��Ԫ�ص�second���ݳ�Ա��Ϊ�գ�����c�ĵ�i��Ԫ�ص�second���ݳ�Ա�ڱ����ֵa����ȡcoords�ĵ�a��Ԫ�أ�������str���뵽coords_to_pdbqt_string��������ӡ������Ϣ�����У������ӡstr������
*/
void model::write_context(const context& c, ofile& out) const {
	verify_bond_lengths();
	VINA_FOR_IN(i, c) {
		const std::string& str = c[i].first;
		if(c[i].second) {
			out << coords_to_pdbqt_string(coords[c[i].second.get()], str) << '\n';
		}
		else
			out << str << '\n';
	}
}

/*
* ����model�ĳ�Ա����seti
* input:conf�ೣ������c
* ouput:void
* ���ã���model��Աatoms, internal_coords, ��c��Աligands����ligands�ĳ�Ա����set_conf    (*this)[i].set_conf(atoms, coords, c[i]);
*       ������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���internal_coords
*/
void model::seti(const conf& c) {
	ligands.set_conf(atoms, internal_coords, c.ligands);
}


/*
* ����model�ĳ�Ա����sete
* input:conf�ೣ������c
* ouput:void
* ���ã�����c.ligands��Ԫ�أ���model��Աinternal_coords, coords, ligands[i].begin��ligand��Ա��, ligands[i].end��ligand��Ա�����뺯��c.ligands[i].rigid.apply���������coords
* ������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords				
*/
void model::sete(const conf& c) {
	VINA_FOR_IN(i, ligands)
		c.ligands[i].rigid.apply(internal_coords, coords, ligands[i].begin, ligands[i].end);
	flex.set_conf(atoms, coords, c.flex);
}

/*
* ����model�ĳ�Ա����set
* input:conf�ೣ������c
* ouput:void
* ������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
*/
void model::set         (const conf& c) {
	ligands.set_conf(atoms, coords, c.ligands);
	flex   .set_conf(atoms, coords, c.flex);
}


/*
* model�ĳ�Ա����gyration_radius
* input:unsigned int����ligand_number
* output:double������
* ���ã���������ligand_number���õ�ligand��ligand_number��Ԫ�أ���ȡ���ĳ�Աbegin,end������begin,end-1��Χ����atoms��Ԫ�صĳ�Աel�Ƿ����EL_TYPE_H
* ��⹫ʽ1�����ع�ʽ1�����ֵ
*/
fl model::gyration_radius(sz ligand_number) const {
	//ligand_number
	VINA_CHECK(ligand_number < ligands.size()); //���ligand_number�Ƿ�С��ligands�����Ĵ�С����С���򱨴���ֹ����ִ��


	const ligand& lig = ligands[ligand_number]; //��ligands������ligand_number��Ԫ�ذ󶨵�lig��

	//�����ڲ�����double���ͱ���acc = 0��unsigned����counter = 0
	fl acc = 0;
	unsigned counter = 0;
	//��lig�ĳ�Աbegin��Сѭ����lig�ĳ�Աend��С-1
	VINA_RANGE(i, lig.begin, lig.end) {
		//���model��Աatoms�ĵ�i��Ԫ�صĳ�Աel������0
		//acc = acc+����coords��i��Ԫ�غ�lig�ĳ�Աnode�ĳ�Աorigin���data����Ԫ�ز��ƽ��
		//counter�Լ�1
		if(atoms[i].el != EL_TYPE_H) { // only heavy atoms are used
			acc += vec_distance_sqr(coords[i], lig.node.get_origin()); // FIXME? check!
			++counter;
		}
	}
	//���model��Աatomsû��һ��Ԫ�����ĳ�Աel����0������0�����򷵻ؼ���ʽ1
	return (counter > 0) ? std::sqrt(acc/counter) : 0;
}


/*
* model�ĳ�Ա����eval_interacting_pairs   �����໥���ö�
* input:precalculate�ೣ������p��interacting_pairs�ೣ������pairs��vecv��(vec����vector)��������coords
* output:double������
* ���ã�1.����p��m_cutoff_sqr��ֵ��cutoff_sqr����ֹ�����������һ��double����e=0
*		2.����pairs��Ԫ��
*				����i��Ԫ�ذ󶨵�ip��
*				����coords[ip.a]��ip.a-->ip�ĳ�Աa����coords[ip.b]��ip.b-->ip�ĳ�Աb����ƽ������(coords[ip.a][0]-coords[ip.b][0])^2+(coords[ip.a][1]-coords[ip.b][1])^2+(coords[ip.a][2]-coords[ip.b][2])^2����ֵ��r2
*				�ж�r2�Ƿ�С��cutoff_sqr���������ִ��������䣬���������һ��ѭ��
*				��ip.type_pair_index��r2����p��eval_fast������p.data.m_data[ip.type_pair_index].fast[sz(factor * r2)](����szΪǿ������ת����factorΪp.data.m_data[ip.type_pair_index]�����ݳ�Ա)��ֵ��tmp
*				��v < epsilon_fl�����б������ļ��������ʶ�����С���㸡������ʱ��tmp = 0������tmp = tmp*(v / (v + tmp))
*				e�ۼ�tmp
*		3.����e
* 
*/
fl eval_interacting_pairs(const precalculate& p, fl v, const interacting_pairs& pairs, const vecv& coords) { // clean up
	const fl cutoff_sqr = p.cutoff_sqr();
	fl e = 0;
	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		fl r2 = vec_distance_sqr(coords[ip.a], coords[ip.b]);
		if(r2 < cutoff_sqr) {
			fl tmp = p.eval_fast(ip.type_pair_index, r2);
			curl(tmp, v);
			e += tmp;
		}
	}
	return e;
}

/*
* model�ĳ�Ա����eval_interacting_pairs_deriv   �����໥���öԵ���
* input:precalculate�ೣ������p��double������v��interacting_pairs�ೣ������pairs��vecv��(vec����vector)��������coords��vecv��(vec����vector)����coords
* output:double������
* ���ã�1.����p��m_cutoff_sqr��ֵ��cutoff_sqr����ֹ�����������һ��double����e=0
*		2.����pairs��Ԫ��
*				����i��Ԫ�ذ󶨵�ip��
*				����coords[ip.a]��ip.a-->ip�ĳ�Աa����coords[ip.b]��ip.b-->ip�ĳ�Աb���Ĳ�r[0] = coords[ip.a][0]-coords[ip.b][0]��r[1] = coords[ip.a][1]-coords[ip.b][1]��r[2] = coords[ip.a][2]-coords[ip.b][2]
*				����r2[0] = r[0]^2��r2[1] = r[1]^2��r2[2] = r[2]^2
*				�ж�r2�Ƿ�С��cutoff_sqr���������ִ��������䣬���������һ��ѭ��
* 
*				��ip.type_pair_index��r2����p��eval_deriv�����ط���pair�ࣨfirs��second���ݳ�Ա��Ϊdouble��pr��
*				pr(p1.first  + ��factor * r2 -  sz(factor * r2)�� * (p2.first  - p1.first), p1.second + ��factor * r2 -  sz(factor * r2)��* (p2.second - p1.second))������p1 = smooth[sz(factor * r2)]��p2 = smooth[sz(factor * r2)+1]��Ϊpair��
*				(����szΪǿ������ת����factorΪp.data.m_data[ip.type_pair_index]�����ݳ�Ա��smoothΪp.data.m_data[ip.type_pair_index]�����ݳ�Ա��pair����vector��)��ֵ��pr
* 
*				force[0] = tmp.second * r[0]��force[1] = tmp.second * r[1]��force[2] = tmp.second * r[2]         tmp.second------------>tmp��second���ݳ�Ա
*				��tmp.first > 0 ���� v <  epsilon_fl�����б������ļ��������ʶ�����С���㸡������ʱ��tmp.first = 0��force[0] = 0��force[1] = 0��force[2] = 0��
*				��tmp.first > 0 ���� v >= epsilon_fl�����б������ļ��������ʶ�����С���㸡������ʱ��force[0] = force[0]*((v / (v + tmp.first))^2)��force[1] = force[1]*((v / (v + tmp.first))^2)��force[2] = force[2]*((v / (v + tmp.first))^2),
*																									   tmp.first = tmp.first*(v / (v + tmp.first))	
*				e�ۼ�tmp.first
*				forces[ip.a][0] = forces[ip.a][0] - force[0],forces[ip.a][1] = forces[ip.a][1] - force[1],forces[ip.a][2] = forces[ip.a][2] - force[2]
*				forces[ip.b][0] = forces[ip.b][0] + force[0],forces[ip.b][1] = forces[ip.b][1] + force[1],forces[ip.b][2] = forces[ip.b][2] + force[2]
*		3.����e
* 
*/
fl eval_interacting_pairs_deriv(const precalculate& p, fl v, const interacting_pairs& pairs, const vecv& coords, vecv& forces) { // adds to forces  // clean up
	const fl cutoff_sqr = p.cutoff_sqr();
	fl e = 0;
	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		vec r; r = coords[ip.b] - coords[ip.a]; // a -> b
		fl r2 = sqr(r);
		if(r2 < cutoff_sqr) {
			pr tmp = p.eval_deriv(ip.type_pair_index, r2);
			vec force; force = tmp.second * r;
			curl(tmp.first, force, v);
			e += tmp.first;
			// FIXME inefficient, if using hard curl
			forces[ip.a] -= force; // we could omit forces on inflex here
			forces[ip.b] += force;
		}
	}
	return e;
}

/*
* model�ĳ�Ա����evali
* input:precalculate�ೣ������p��vec�ೣ������v
* output:double������
* ���ã�1.����һ��double����e=0
*		2.����ligands��Ԫ��				
*				��p,v[0],ligands��i��Ԫ�صĳ�Աpairs,internal_coords����eval_interacting_pairs����������ֵ�ۼӵ�e��
*		3.����e
*/
fl model::evali(const precalculate& p,                                  const vec& v                          ) const { // clean up
	fl e = 0;
	VINA_FOR_IN(i, ligands) 
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, internal_coords); // probably might was well use coords here
	return e;
}

/*
* model�ĳ�Ա����evale
* input:precalculate�ೣ������p��igrid�ೣ������ig��vec�ೣ������v
* output:double������
* ���ã�1.����һ��double����e,����鿴����ig�е�eval����,igӦ����igrid���������
*		2.����ligands��Ԫ��
*				��p,v[2],other_pairs, coords����eval_interacting_pairs����������ֵ�ۼӵ�e��
*		3.����e
*/
fl model::evale(const precalculate& p, const igrid& ig, const vec& v                          ) const { // clean up
	fl e = ig.eval(*this, v[1]);
	e += eval_interacting_pairs(p, v[2], other_pairs, coords);
	return e;
}


/*
* model�ĳ�Ա����eval
* input:precalculate�ೣ������p��igrid�ೣ������ig��vec�ೣ������v��conf�ೣ������c
* output:double������
* ���ã�1.������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
*		2.��p,ig,v����evale����������ֵ��ֵ��e��
*		3.����ligands��Ԫ��
*				��p,v[0],ligands��i��Ԫ�صĳ�Աpairs,coords����eval_interacting_pairs����������ֵ�ۼӵ�e��
*		4.����e
*/
fl model::eval         (const precalculate& p, const igrid& ig, const vec& v, const conf& c           ) { // clean up
	set(c);
	fl e = evale(p, ig, v);
	VINA_FOR_IN(i, ligands) 
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords
	return e;
}

/*
* model�ĳ�Ա����eval_deriv
* input:precalculate�ೣ������p��igrid�ೣ������ig��vec�ೣ������v��conf�ೣ������c��change������g
* output:double������
* ���ã�1.������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
*		2.����һ��double����e,����鿴����ig�е�eval_deriv����,igӦ����igrid���������
*		3.��p,v[2], other_pairs, coords,minus_forces����eval_interacting_pairs_deriv����������ֵ�ۼӵ�e��
*		4.����ligands��Ԫ��
*				��p,v[0],ligands��i��Ԫ�صĳ�Աpairs,coord,minus_forces����eval_interacting_pairs_deriv����������ֵ�ۼӵ�e��
*		5.��coords, minus_forces, g.ligands����ligands.derivative�У��ֱ��ligands������Ԫ��������ligands[i].derivative(coords, minus_forces, g.ligands[i])
*		6.��coords, minus_forces, g.flex����flex.derivative�У��ֱ��flex������Ԫ��������flex[i].derivative(coords, minus_forces, g.flex[i])
*		7.����e
*/
fl model::eval_deriv  (const precalculate& p, const igrid& ig, const vec& v, const conf& c, change& g) { // clean up
	set(c);
	fl e = ig.eval_deriv(*this, v[1]); // sets minus_forces, except inflex
	e += eval_interacting_pairs_deriv(p, v[2], other_pairs, coords, minus_forces); // adds to minus_forces
	VINA_FOR_IN(i, ligands)
		e += eval_interacting_pairs_deriv(p, v[0], ligands[i].pairs, coords, minus_forces); // adds to minus_forces
	// calculate derivatives
	ligands.derivative(coords, minus_forces, g.ligands);
	flex   .derivative(coords, minus_forces, g.flex); // inflex forces are ignored
	return e;
}




/*
* model�ĳ�Ա����eval_intramolecular
* input:precalculate�ೣ������p��vec�ೣ������v��conf�ೣ������c
* output:double������
* ���ã�1.������������c��������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
*		2.����һ��double����e = 0
*		3.����ligands��Ԫ��
*				��p,v[0],ligands��i��Ԫ�صĳ�Աpairs,coords����eval_interacting_pairs_deriv����������ֵ�ۼӵ�e��
*		4.��m_atom_typing_used����num_atom_types�õ�m_atom_typing_used��Ӧ��TYPE_SIZE����ֵ��nat
*		  ����p.m_cutoff_sqr��ֵ��cutoff_sqr
*		5.i��0ѭ����m_num_movable_atoms-1
*				�ж�i�Ƿ���ligands��������һ��Ԫ�صĳ�Աbegin��end��Χ֮�ڣ�����ֱ����һ��iѭ�����������ִ���������		
*				��atmos�ĵ�i��Ԫ�ذ���a��
*				����m_atom_typing_used��Ӧ��a��unit������el, ad, xs, sy֮һ����ֵ��t1
*				��t1>=natʱֱ����һ��iѭ�����������ִ���������
*				j��0ѭ����grid_atoms.size()-1������grid_atoms��Ԫ�أ�
*						��atmos�ĵ�j��Ԫ�ذ���b��
*						����m_atom_typing_used��Ӧ��b��unit������el, ad, xs, sy֮һ����ֵ��t2
*						��t2>=natʱֱ����һ��jѭ�����������ִ���������
*						r2 = (coords[i][0]-b.coords[0])^2 +(coords[i][1]-b.coords[1])^2 + (coords[i][2]-b.coords[2])^2
*						���r2 < cutoff_sqr
*								���t1<=t2��type_pair_index = t1 + t2*(t2+1)/2�����t1>t2��type_pair_index = t2 + t1*(t1+1)/2
*								this_e = p.data.m_data[type_pair_index].fast[sz(factor * r2)]��factorΪp.data.m_data[type_pair_index]�ĳ�Ա��szΪǿ������ת��Ϊunsigned int		
*								��this_e>0����v[1]<0.1*max_fl(��ǰ��������double�ͱ������ֵ)ʱ�����v[1] < epsilon_fl=2.22045e-16��this_e=0������this_e=this_e*(v[1] / (v[1] + this_e))
*								��this_e�ۼӵ�e��
* 
*		6.i��0ѭ����other_pairs.size-1
*				��other_pairs�ĵ�interacting_pair��Ԫ�ذ���a��
*				���pair.a��ligands��������һ��Ԫ�صĳ�Աbegin��end��Χ֮�ڻ���pair.b��ligands��������һ��Ԫ�صĳ�Աbegin��end��Χ֮�ڣ���ֱ����һ��iѭ�����������ִ���������	
*				��r2 = ��coords[pair.a][0] - coords[pair.b][0]��^2 + ��coords[pair.a][1] - coords[pair.b][1]��^2 + ��coords[pair.a][2] - coords[pair.b][2]��^2
*				���r2 < cutoff_sqr
*						this_e = p.data.m_data[pair.type_pair_index].fast[sz(factor * r2)]��factorΪp.data.m_data[type_pair_index]�ĳ�Ա��szΪǿ������ת��Ϊunsigned int		
*						��this_e>0����v[2]<0.1*max_fl(��ǰ��������double�ͱ������ֵ)ʱ�����v[2] < epsilon_fl=2.22045e-16��this_e=0������this_e=this_e*(v[2] / (v[2] + this_e))
*						��this_e�ۼӵ�e��
*		7.����e
*/
fl model::eval_intramolecular(const precalculate& p, const vec& v, const conf& c) {
	set(c);
	fl e = 0;

	// internal for each ligand
	VINA_FOR_IN(i, ligands)
		e += eval_interacting_pairs(p, v[0], ligands[i].pairs, coords); // coords instead of internal coords

	sz nat = num_atom_types(atom_typing_used());
	const fl cutoff_sqr = p.cutoff_sqr();

	// flex-rigid
	VINA_FOR(i, num_movable_atoms()) {
		if(find_ligand(i) < ligands.size()) continue; // we only want flex-rigid interaction
		const atom& a = atoms[i];
		sz t1 = a.get(atom_typing_used());
		if(t1 >= nat) continue;
		VINA_FOR_IN(j, grid_atoms) {
			const atom& b = grid_atoms[j];
			sz t2 = b.get(atom_typing_used());
			if(t2 >= nat) continue;
			fl r2 = vec_distance_sqr(coords[i], b.coords);
			if(r2 < cutoff_sqr) {
				sz type_pair_index = triangular_matrix_index_permissive(nat, t1, t2);
				fl this_e = p.eval_fast(type_pair_index, r2);
				curl(this_e, v[1]);
				e += this_e;
			}
		}
	}

	// flex-flex
	VINA_FOR_IN(i, other_pairs) {
		const interacting_pair& pair = other_pairs[i];
		if(find_ligand(pair.a) < ligands.size() || find_ligand(pair.b) < ligands.size()) continue; // we only need flex-flex
		fl r2 = vec_distance_sqr(coords[pair.a], coords[pair.b]);
		if(r2 < cutoff_sqr) {
			fl this_e = p.eval_fast(pair.type_pair_index, r2);
			curl(this_e, v[2]);
			e += this_e;
		}
	}
	return e;
}

/*
* ����model�ĳ�Ա����eval_adjusted 
* input:scoring_function�ೣ������sf��precalculate�ೣ������p��igrid�ೣ������ig��vec�ೣ������v��conf�ೣ������c��double������intramolecular_energy��������������
* ouput:double������
* ���ã�1.��p, ig, v, c����eval�����л�ù���ֵe
*		2.����model��e-intramolecular_energy����sf��Ա����conf_independent�����ض�����������ֵ
*		
*/
fl model::eval_adjusted      (const scoring_function& sf, const precalculate& p, const igrid& ig, const vec& v, const conf& c, fl intramolecular_energy) {
	fl e = eval(p, ig, v, c); // sets c
	return sf.conf_independent(*this, e - intramolecular_energy);
}

/*
* ����model�ĳ�Ա����rmsd_lower_bound_asymmetric �½�ǶԳƾ��������
* input:model�ೣ������x,y
* ouput:double������
* ���ã�1.��ȡx�ĳ�Աm_num_movable_atoms��ֵ��n���ж�n�Ƿ����y�ĳ�Աm_num_movable_atoms��������򱨴���ֹ����ִ��
*		2.����һ��double������sum = 0,unsigned������counter = 0
*		3.i��0ѭ����n-1
*			��x�ĳ�Աatoms�ĵ�i��Ԫ�ذ󶨵�a��
*			���a�ĳ�Աel������0
*				��1������һ��double������r2���ڵ�ǰ��������double�ͱ������ֵ
*				��2��j��0ѭ����n-1
*						��y�ĳ�Աatoms�ĵ�j��Ԫ�ذ󶨵�b��
*						���a�ĳ�Աel��b�ĳ�Աel��Ȳ���b������ԭ�ӣ�!(b.ad==6||b.ad==12)��
*							[1]����һ��double������this_r2����x��Աcoords��i��Ԫ����y��Աcoords��j��Ԫ�����data����Ԫ�ز��ƽ��
*							[2]���this_r2 < r2����this_r2��ֵ��r2
*				�ж�r2<0.1*max_fl(��ǰ��������double�ͱ������ֵ),�����С�ڱ�����ֹ����ִ��
*				��r2�ۼӵ�sum�ϣ�counter����һ
*		4.���counter����0������0�����򷵻�sum / counter�Ŀ�����	
*/

fl model::rmsd_lower_bound_asymmetric(const model& x, const model& y) const { // actually static
	sz n = x.m_num_movable_atoms; 
	VINA_CHECK(n == y.m_num_movable_atoms);
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR(i, n) {
		const atom& a =   x.atoms[i];
		if(a.el != EL_TYPE_H) {
			fl r2 = max_fl;
			VINA_FOR(j, n) {
				const atom& b = y.atoms[j];
				if(a.same_element(b) && !b.is_hydrogen()) {
					fl this_r2 = vec_distance_sqr(x.coords[i], 
					                              y.coords[j]);
					if(this_r2 < r2)
						r2 = this_r2;
				}
			}
			assert(not_max(r2));
			sum += r2;
			++counter;
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

/*
* ����model�ĳ�Ա����rmsd_lower_bound  ���޾��������
* input:model�ೣ������m
* ouput:double������
* ���ã�1.��*this����model����m��˳�����뺯��rmsd_lower_bound_asymmetric����m��*this����model����˳�����뺯��rmsd_lower_bound_asymmetric����������������ֵ�ϴ��
*/
fl model::rmsd_lower_bound(const model& m) const {
	return (std::max)(rmsd_lower_bound_asymmetric(*this,     m),
		            rmsd_lower_bound_asymmetric(    m, *this));
}

/*
* ����model�ĳ�Ա����rmsd_upper_bound  ���޾��������
* input:model�ೣ������m
* ouput:double������
* ���ã�1.��鱾model��Աm_num_movable_atoms��m�ĳ�Աm_num_movable_atoms�Ƿ���ȣ�������򱨴���ֹ����ִ��
*		2.����һ��double������sum = 0,unsigned������counter = 0
*		3.i��0ѭ��������model��Աm_num_movable_atoms-1��
*			����model�ĳ�Աatoms�ĵ�i��Ԫ�ذ󶨵�a��
*			��m�ĳ�Աatoms�ĵ�i��Ԫ�ذ󶨵�b��
*			���a��Աad��b�ĳ�Աad�Ƿ���ȣ�������򱨴���ֹ����ִ��
*			���a��Աxs��b�ĳ�Աxs�Ƿ���ȣ�������򱨴���ֹ����ִ��
*			���a��Աel������0
*				��1������model��Աcoords��i��Ԫ����m��Աcoords��i��Ԫ�����data����Ԫ�ز��ƽ���ۼӵ�sum��
*				��2��counter����һ
*		4.���counter����0������0�����򷵻�sum / counter�Ŀ�����
*/
fl model::rmsd_upper_bound(const model& m) const {
	VINA_CHECK(m_num_movable_atoms == m.m_num_movable_atoms);
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR(i, m_num_movable_atoms) {
		const atom& a =   atoms[i];
		const atom& b = m.atoms[i];
		assert(a.ad == b.ad);
		assert(a.xs == b.xs);
		if(a.el != EL_TYPE_H) {
			sum += vec_distance_sqr(coords[i], m.coords[i]);
			++counter;
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

/*
* ����model�ĳ�Ա����rmsd_ligands_upper_bound �������޾��������
* input:model�ೣ������m
* ouput:double������
* ���ã�1.��鱾model��Աligands�����Ĵ�С��m�ĳ�Աligands�����Ĵ�С�Ƿ���ȣ�������򱨴���ֹ����ִ��
*		2.����һ��double������sum = 0,unsigned������counter = 0
*		3.ligand_i��0ѭ��������model��Աligands�Ĵ�С-1��
*			����model�ĳ�Աligands�ĵ�ligand_i��Ԫ�ذ󶨵�lig  ��
*			��  m    �ĳ�Աligands�ĵ�ligand_i��Ԫ�ذ󶨵�m_lig��
*			���lig��Աbegin��m_lig�ĳ�Աbegin �Ƿ���ȣ�������򱨴���ֹ����ִ��
*			���lig��Աend  ��m_lig�ĳ�Աend   �Ƿ���ȣ�������򱨴���ֹ����ִ��
*			i��lig��Աbeginѭ������lig��Աend-1��
*				 ����model�ĳ�Աatmos�ĵ�i��Ԫ�ذ󶨵�a��
* 				 ��  m    �ĳ�Աatmos�ĵ�i��Ԫ�ذ󶨵�b��
*				 ���a��Աad��b�ĳ�Աad�Ƿ���ȣ�������򱨴���ֹ����ִ��
* 				 ���a��Աxs��b�ĳ�Աxs�Ƿ���ȣ�������򱨴���ֹ����ִ��
*				 ���a��Աel������0
*				 	��1������model��Աcoords��i��Ԫ����m��Աcoords��i��Ԫ�����data����Ԫ�ز��ƽ���ۼӵ�sum��
* 				 	��2��counter����һ
*		4.���counter����0������0�����򷵻�sum / counter�Ŀ�����
*/
fl model::rmsd_ligands_upper_bound(const model& m) const {
	VINA_CHECK(ligands.size() == m.ligands.size());
	fl sum = 0;
	unsigned counter = 0;
	VINA_FOR_IN(ligand_i, ligands) {
		const ligand&   lig =   ligands[ligand_i];
		const ligand& m_lig = m.ligands[ligand_i];
		VINA_CHECK(lig.begin == m_lig.begin);
		VINA_CHECK(lig.end   == m_lig.end);
		VINA_RANGE(i, lig.begin, lig.end) {
			const atom& a =   atoms[i];
			const atom& b = m.atoms[i];
			assert(a.ad == b.ad);
			assert(a.xs == b.xs);
			if(a.el != EL_TYPE_H) {
				sum += vec_distance_sqr(coords[i], m.coords[i]);
				++counter;
			}
		}
	}
	return (counter == 0) ? 0 : std::sqrt(sum / counter);
}

/*
* ����model�ĳ�Ա����verify_bond_lengths ��֤���ӳ���
* input:��
* ouput:��
* ���ã�i��0ѭ������model��Աgrid_atoms�Ĵ�С��atoms�Ĵ�С֮��-1��
*				��i���뵽model��Ա����sz_to_atom_index�����ص�ֵ��ֵ��ai�ϣ���0<=i<=(grid_atoms�Ĵ�С-1)ʱ��aiΪatom_index(i,true)������Ϊatom_index(i - grid_atoms.size(), false)
* 				��ai����model��Ա����get_atom��aiΪatom_index(i,true)ʱ�����model��Աgrid_atoms�ĵ�i��ai�ĳ�Աi����Ԫ�أ�aiΪatom_index(i - grid_atoms.size(), false)ʱ�����model��Աatoms�ĵ�i��ai�ĳ�Աi����Ԫ�أ����ص�ֵ�󶨵�a��
* 				j��0ѭ����a��Աbonds�Ĵ�С-1
*					a��Աbonds�ĵ�j��Ԫ�ذ󶨵�b��
*					��ai��b�ĳ�Աconnected_atom_index�������ƽ���������õ�����ֵ�����㿪�������õ�d������ʽ3��
* 					�Ƚ�d��b�ĳ�Աlength���������okΪ�棬����Ϊ��
*					���okΪ��
*						�ڿ���̨��ʾd=d����ֵ
*						�ڿ���̨��ʾb.length=b.length����ֵ
*					���ok�Ƿ�Ϊ�棬��Ϊ���򱨴���ֹ����ִ��
*						
*		
*/
void model::verify_bond_lengths() const {
	VINA_FOR(i, grid_atoms.size() + atoms.size()) {
		const atom_index ai = sz_to_atom_index(i);
		const atom& a = get_atom(ai);
		VINA_FOR_IN(j, a.bonds) {
			const bond& b = a.bonds[j];
			fl d = std::sqrt(distance_sqr_between(ai, b.connected_atom_index));
			bool ok = eq(d, b.length);
			if(!ok) {
				VINA_SHOW(d);
				VINA_SHOW(b.length);
			}
			VINA_CHECK(ok);
		}
	}
}

/*
* ����model�ĳ�Ա����check_internal_pairs ��֤���ӳ���
* input:��
* ouput:��
* ���ã�����model��Ա����ligands��Ԫ��
*			��ligands�ĵ�i��Ԫ�ذ󶨵�lig��
* 			��lig�ĳ�Աpairs�󶨵�pairs��
*			��������pairs��Ԫ��
*				��pairs�ĵ�j��Ԫ�ذ󶨵�ip��
*				���ip�ĳ�Աa�Ƿ���ڵ���lig�ĳ�Աbegin��С���򱨴���ֹ����ִ��
* 				���ip�ĳ�Աb�Ƿ�С��lig�ĳ�Աend����С���򱨴���ֹ����ִ��
*				
*
* ��������begin�Ƿ�С�ڵ����໥���öԵ�a����������end�Ƿ�����໥���öԵ�b
* 
*/
void model::check_internal_pairs() const {
	VINA_FOR_IN(i, ligands) {
		const ligand& lig = ligands[i];
		const interacting_pairs& pairs = lig.pairs;
		VINA_FOR_IN(j, pairs) {
			const interacting_pair& ip = pairs[j];
			VINA_CHECK(ip.a >= lig.begin);
			VINA_CHECK(ip.b  < lig.end);
		}
	}
}

/*
* model�ĳ�Ա����:about
* input:��
* output:��
* ���ã��ڿ���̨��ʾatom_typing_used()   = model��Աm_atom_typing_used��ֵ
* 		�ڿ���̨��ʾnum_movable_atoms()  = model��Աm_num_movable_atoms��ֵ
*		�ڿ���̨��ʾnum_internal_pairs() = model��Աligands����Ԫ�صĳ�Աpairs�Ĵ�С֮��
* 		�ڿ���̨��ʾnum_other_pairs()	 = model��Աother_pairs�Ĵ�С
*		�ڿ���̨��ʾnum_ligands()		 = model��Աligands�Ĵ�С
* 		�ڿ���̨��ʾnum_flex()		     = model��Աflex�Ĵ�С
*
*/
void model::about() const {
	VINA_SHOW(atom_typing_used());
	VINA_SHOW(num_movable_atoms());
	VINA_SHOW(num_internal_pairs());
	VINA_SHOW(num_other_pairs());
	VINA_SHOW(num_ligands());
	VINA_SHOW(num_flex());
}

/*
* model�ĳ�Ա����:about
* input:��
* output:��
* ���ã�
*		����coords��Ԫ�أ��ڿ���̨��ӡa�ĳ�Ա����
* 		��ӡ
* 			coords:
* 			[coords[i][0],coords[i][1],coords[i][2]]
*			........................................
*		����internal_coords��Ԫ�أ��ڿ���̨��ӡa�ĳ�Ա����
* 		��ӡ
* 			internal_coords:
* 			[internal_coords[i][0],internal_coords[i][1],internal_coords[i][2]]
* 			...................................................................
*		����atoms��Ԫ�أ��ֱ����a�ϣ��ڿ���̨��ӡa�ĳ�Ա����
* 		��ӡ
* 			atoms:
* 			a.el a.ad a.xs a.sy a.charge
* 			a.bonds�Ĵ�С [a.coords[0],a.coords[1],a.coords[2]]       
*			...................................................
*		����grid_atoms��Ԫ�أ��ֱ����a�ϣ��ڿ���̨��ӡa�ĳ�Ա����
*		��ӡ
*			grid_atoms:
*			a.el a.ad a.xs a.sy a.charge
*			a.bonds�Ĵ�С [a.coords[0],a.coords[1],a.coords[2]]        
* 		    ...................................................
*/
void model::print_stuff() const {
	std::cout << "coords:\n";
	VINA_FOR_IN(i, coords)
		printnl(coords[i]);

	std::cout << "internal_coords:\n";
	VINA_FOR_IN(i, internal_coords)
		printnl(internal_coords[i]);

	std::cout << "atoms:\n";
	VINA_FOR_IN(i, atoms) {
		const atom& a = atoms[i];
		std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
		std::cout << a.bonds.size() << "  "; printnl(a.coords);
	}

	std::cout << "grid_atoms:\n";
	VINA_FOR_IN(i, grid_atoms) {
		const atom& a = grid_atoms[i];
		std::cout << a.el << " " << a.ad << " " << a.xs << " " << a.sy << "    " << a.charge << '\n';
		std::cout << a.bonds.size() << "  "; printnl(a.coords);
	}
	about();
}


/*
*input:double��r��covalent_r
*output���� r < 0 ���� covalent < ���б������ļ��������ʶ�����С���㸡���� ʱ����
*����r / covalent_r >2 ʱ����0
*    r / covalent_r <= 2 ʱ����1-��r / covalent_r��*��r / covalent_r��/4
*/
fl pairwise_clash_penalty(fl r, fl covalent_r) {
	// r = 0          -> max_penalty 
	// r = covalent_r -> 1
	// elsewhere      -> hyperbolic function
	assert(r >= 0);//r<0����
	assert(covalent_r > epsilon_fl);//covalent<���б������ļ��������ʶ�����С���㸡����ʱ����
	const fl x = r / covalent_r;
	if(x > 2) return 0;
	return 1-x*x/4;
}

/*
* model��Ա����clash_penalty_aux����ͻ�ͷ�
*input:interacting_pair��vector
*output:doule����
*���ã�
*	1.����һ��double����e = 0
*	2.����vector
*		(1)��vector�ĵ�i��Ԫ�ذ󶨵�ip��
*       (2)���㹫ʽ2���õ�r
*       (3)����atoms��ip.a��ip�ĳ�Աa����ip.b��ip�ĳ�Աb���Ĺ��۰뾶�ͣ��õ�covalent_r
*       (4)����ɶԳ�ͻ�ͷ�����r��covalent_r���뺯��pairwise_clash_penalty���������ۼӵ�e��
*   3.����e
*/
fl model::clash_penalty_aux(const interacting_pairs& pairs) const {
	fl e = 0;
	VINA_FOR_IN(i, pairs) {
		const interacting_pair& ip = pairs[i];
		const fl r = std::sqrt(vec_distance_sqr(coords[ip.a], coords[ip.b]));
		const fl covalent_r = atoms[ip.a].covalent_radius() + atoms[ip.b].covalent_radius();
		e += pairwise_clash_penalty(r, covalent_r);
	}
	return e;
}

/*
* model��Ա����clash_penalty����ͻ�ͷ�
*input:��
*output:doule����
*���ã�
*	1.����һ��double����e = 0
*	2.����model��Աligands��
*         ��ligands��i��Ԫ�صĳ�Աpairs����clash_penalty_aux���õ�����ۼӵ�e��
*   3.��model��Աother_pairs����clash_penalty_aux���õ�����ۼӵ�e��
*   4.����e
*/
fl model::clash_penalty() const {
	fl e = 0;
	VINA_FOR_IN(i, ligands) 
		e += clash_penalty_aux(ligands[i].pairs);
	e += clash_penalty_aux(other_pairs);
	return e;
}
