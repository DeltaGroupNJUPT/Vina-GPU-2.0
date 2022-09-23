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

#ifndef VINA_TERMS_H
#define VINA_TERMS_H

#include <boost/ptr_container/ptr_vector.hpp> 
#include "model.h"

/*
	����һ���ַ���name��
	virtual ~term()�麯����
	*/

struct term {
	std::string name;
	virtual ~term() {}
};


/*	distance_additive�̳нṹ��term�����Է��ʸ�ֵterm�е�public��Ա��
	����һ��double cutoff��
	������һ��double�麯��eval,����Ϊ����atom_base�ṹ��a��b��һ��double r;
	virtual ~distance_additive()�麯��
*/

struct distance_additive : public term {                                        //�̳�term�ṹ���public
	fl cutoff;
	distance_additive(fl cutoff_) : cutoff(cutoff_) {}
	virtual fl eval(const atom_base& a, const atom_base& b, fl r) const = 0;
	virtual ~distance_additive() {}
};


/*	usable�̳нṹ��distance_additive�����Է��ʺ͸�ֵdistance_additive��public��Ա��
	����ýṹ����ó�ʼ������double cuoff��ֵ��distance_additive�ṹ���е�double cuoff��
	���Զ���Ϊatom_type����ö������t��atom_typing_used��ֵΪtö��Ԫ���е�XS=2��
	������eval����������Ϊһ��atom_base�ṹ��a��atom_base�ṹ���Լ�double r���ṹ��atom_base�̳��˽ṹ��atom_type��
	����a.get����������atom_type�ĳ�Ա�����Ե�֪a.get(XS)=b.get(XS)=XS_TYPE_SIZE��double r=r��
	������һ��double�麯��eval������д,����Ϊ����unint t1��t2��һ��double r��
	virtual ~usable()�麯��

*/
struct usable : public distance_additive {
	atom_type::t atom_typing_used;             //���Ǳ�˵��Ϊatom_type����ö������t�����ͱ�����ȡֵֻ����ö��Ԫ���е�ĳһ��ֵ��
	usable(fl cutoff_) : distance_additive(cutoff_), atom_typing_used(atom_type::XS) {}     //XS=2
	fl eval(const atom_base& a, const atom_base& b, fl r) const { // should not be overriden
		return eval(a.get(atom_typing_used), b.get(atom_typing_used), r);       //get���������ǽṹ��atom_base��Ա��
																				//��Ϊatom_typing_use=XS�����Է���XS_TYPE_SIZE
	}
	virtual fl eval(sz t1, sz t2, fl r) const { return 0; } 
	virtual ~usable() {}
};



/* additive�̳нṹ��term
����һ��double cutoff��Ȼ�󱻸�ֵΪ�������ʶ���double���ɱ�ʾ�����ֵ00DB17EE��
������һ��double�麯��eval,����Ϊmodel�ṹ��m������atom_base�ṹ��i��j;
virtual ~additive() �麯��
*/

struct additive : public term {
	fl cutoff;
	additive() : cutoff(max_fl) {}                                    //�������ʶ���double���ɱ�ʾ�����ֵ���ҵļ������ʶ������double����00DB17EE
	virtual fl eval(const model& m, const atom_index& i, const atom_index& j) const = 0;
	virtual ~additive() {}
};



/* intermolecular�̳нṹ��term

������һ��double�麯��eval,����Ϊmodel�ṹ��m��

*/


struct intermolecular : public term {
	virtual fl eval(const model& m) const = 0;
};




/*
		������num_tors��num_rotors��num_heavy_atoms��num_hydrophobic_atoms��
		ligand_max_num_h_bonds��num_ligands��ligand_lengths_sum��double�������ݣ�
		�����˺���operator flv()�ͺ��� get_names()
		input����m��model�ṹ�壻
		������num_bonded_heavy_atoms������atom_rotors()����������û�ж��壻
		Ȼ���������ʽ���л������Խ�num_tors��num_rotors��num_heavy_atoms��num_hydrophobic_atoms��
		ligand_max_num_h_bonds��num_ligands��ligand_lengths_sum��Щ�ṹ���Ա�洢���ļ��л���ͨ�����緢��ȥ��
		����ʱ���������������͡�
		https://blog.csdn.net/zj510/article/details/8105408
*/
struct conf_independent_inputs {
	fl num_tors;
	fl num_rotors;
	fl num_heavy_atoms;
	fl num_hydrophobic_atoms;
	fl ligand_max_num_h_bonds;
	fl num_ligands;
	fl ligand_lengths_sum;
	operator flv() const;                        //��������operator flv()
	conf_independent_inputs(const model& m);
	std::vector<std::string> get_names() const;  //�������� get_names()
	conf_independent_inputs();
private:
	unsigned num_bonded_heavy_atoms(const model& m, const atom_index& i) const; // FIXME? - could be static, but I don't feel like declaring function friends
	unsigned atom_rotors(const model& m, const atom_index& i) const; // the number of rotatable bonds to heavy ligand atoms

	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & num_tors;
		ar & num_rotors;
		ar & num_heavy_atoms;
		ar & num_hydrophobic_atoms;
		ar & ligand_max_num_h_bonds;
		ar & num_ligands;
		ar & ligand_lengths_sum;
	}
};


/*
���������麯��eval()��size()�����麯���� 
����in:conf_independent_inputs�ṹ�壻x��double��it����������
flv::const_iterator ���������壬��������һ�ּ��vector��Ԫ�ز�����Ԫ�ص��������ͣ�
const_iteratorֻ�ܶ�ȡvector�е�Ԫ�أ��������ڸı���ֵ��
*/
struct conf_independent : public term {
	virtual fl eval(const conf_independent_inputs& in, fl x, flv::const_iterator& it) const = 0;
	virtual sz size() const = 0; // how many parameters does it take
}; 



/*
inter����enabled��������������fun���������͵�ָ�������������vector���Ʋ���Ч�ʸ��ߣ�
			 term_set�̳�distance_additive��
			 term_set�̳�usable��
			 term_set�̳�additive��
			 term_set�̳�intermolecular��
			 term_set�̳�conf_independent��
*/

template<typename T>
struct term_set {
	std::vector<bool> enabled;
	boost::ptr_vector<T> fun; // FIXME? const T?    //�������͵�ָ������fun,�����vector���Ʋ���Ч�ʸ���

	/*
		input����e:unsigned��f������������ָ�룬enabled��������������fun���������͵�ָ�������������vector���Ʋ���Ч�ʸ��ߣ�
		�ú�����������enabled������ĩβ����0��eС�ڵ���0����1��e����0����
		����fun��ĩβ�������f��
	*/
	void add(unsigned e, T* f) { // FIXME? const T* ?
		enabled.push_back(e > 0);      //������ĩβ����bool��������
		fun.push_back(f);
	}

	/*
	inter����enabled��������������
	output����tmp:unsigned int ��ʼֵΪ0��
	�ú�������enabled�д��ڵ���1��ֵ�ĸ�����
	*/
	sz num_enabled() const {
		sz tmp = 0;
		VINA_FOR_IN(i, enabled)  // for (uint VINA_MACROS_TMP = enabled.size, i = 0; i<VINA_MACROS_TMP; i++)
			if(enabled[i])      //����enabled������Ϊ1��ֵ
				++tmp;
		return tmp;
	}

/*
input����enabled_only��bool��out�����ַ���������
enabled��������������fun���������͵�ָ�������������vector���Ʋ���Ч�ʸ��ߣ�
��ȷ��enabled������Ԫ�ظ����Ƿ��fun����Ԫ����ͬ���ǵĻ�����������ֹ��
������enabled_onlyΪ0����enabled[i]���ڵ���1����out�ַ�������ĩβ���fun[i].name��term�ṹ����name��
*/
	void get_names(bool enabled_only, std::vector<std::string>& out) const { // appends to "out"
		VINA_CHECK(enabled.size() == fun.size());       //assert
		VINA_FOR_IN(i, fun)								// for (uint VINA_MACROS_TMP = fun.size, i = 0; i<VINA_MACROS_TMP; i++)
			if(!enabled_only || enabled[i])             
				out.push_back(fun[i].name);              //term�ṹ����name
	}



	/*
	input����in:const_iterator��������out��double����
	inter����enabled��������������fun���������͵�ָ�������������vector���Ʋ���Ч�ʸ��ߣ�
	��ȷ��enabled������Ԫ�ظ����Ƿ��fun����Ԫ����ͬ���ǵĻ�����������ֹ��
	���enabled�е�ֵ���ڵ���1ʱ�����������iԪ����ӵ�out������ĩβ��
	in++ѭ��enabled.size�Σ�����if�����
	*/
	void filter(flv::const_iterator& in, flv& out) const {
		VINA_CHECK(enabled.size() == fun.size());   // for (uint VINA_MACROS_TMP = enabled.size, i = 0; i<VINA_MACROS_TMP; i++)
		VINA_FOR_IN(i, enabled) {
			if(enabled[i])
				out.push_back(*in);                 //�����������iԪ����ӵ�out������ĩβ
			++in;
		}
	}


	/*
	inter����tmp��double
	enabled��������������fun���������͵�ָ�������������vector���Ʋ���Ч�ʸ��ߣ�
	ͨ��һ��ѭ��Ѱ��fun[i].cutoff�е����ֵ�����أ�additive�ṹ����cutoff��
	*/
	fl max_cutoff() const {
		fl tmp = 0; 
		VINA_FOR_IN(i, fun)                          // for (uint VINA_MACROS_TMP = fun.size, i = 0; i<VINA_MACROS_TMP; i++)
			tmp = (std::max)(tmp, fun[i].cutoff);   //additive�ṹ����cutoff
		return tmp;
	}

	/*
	    fun���������͵�ָ�������������vector���Ʋ���Ч�ʸ���
		�ú�������fun����Ԫ�ظ���
	*/
	sz size() const { return fun.size(); }

	/*
		input����unint i��
		fun���������͵�ָ�������������vector���Ʋ���Ч�ʸ���
		����fun[i]��ֵ��
	
	*/
	const T& operator[](sz i) const { return fun[i]; }
};




/*
	������double���͵�����e��i��
	����size()��������e������i������Ԫ���ܸ�����
	����num_weights��������������e��i����������Ԫ�ظ�����
	������eval������������������Ϊdouble��
	Ȼ���������ʽ���л������Խ�e��i��Щ�ṹ���Ա�洢���ļ��л���ͨ�����緢��ȥ��
	����ʱ���������������͡�
	https://blog.csdn.net/zj510/article/details/8105408

*/
struct factors {
	flv e; // external
	flv i; // internal
	sz size() const { return e.size() + i.size(); }       //����e������i������Ԫ���ܸ���
	//sz num_weights() const { return (std::max)(e.size(), i.size()); } // FIXME? compiler bug? getting warnings here
	sz num_weights() const { return (e.size() > i.size()) ? e.size() : i.size(); }
	fl eval(const flv& weights, bool include_internal) const;      //����
private:
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive & ar, const unsigned version) {
		ar & e;
		ar & i;
	}
};


/*
			 term_set�̳�distance_additive��Ȼ����һ��distance_additive_terms�ṹ��
			 term_set�̳�usable��Ȼ����һ��usable_terms�ṹ��
			 term_set�̳�additive��Ȼ����һ��additive_terms�ṹ��
			 term_set�̳�intermolecular��Ȼ����һ��intermolecular_terms�ṹ��
			 term_set�̳�conf_independent��Ȼ����һ��conf_independent_terms�ṹ��

*/

struct terms {
	term_set<distance_additive> distance_additive_terms;    //term_set�̳�distance_additive��Ȼ����һ��distance_additive_terms�ṹ��
	term_set<usable>            usable_terms;				//term_set�̳�usable��Ȼ����һ��usable_terms�ṹ��
	term_set<additive>          additive_terms;				//term_set�̳�additive��Ȼ����һ��additive_terms�ṹ��
	term_set<intermolecular>    intermolecular_terms;       //term_set�̳�intermolecular��Ȼ����һ��intermolecular_terms�ṹ��
	term_set<conf_independent>  conf_independent_terms;     //term_set�̳�conf_independent��Ȼ����һ��conf_independent_terms�ṹ��    

	// the class takes ownership of the pointer with 'add'
	void add(unsigned e, distance_additive* p) { distance_additive_terms.add(e, p); }//	��fun��ĩβ���p����enable��ĩβ���0��1
	void add(unsigned e, usable* p)            {            usable_terms.add(e, p); }//	��fun��ĩβ���p����enable��ĩβ���0��1
	void add(unsigned e, additive* p)          {          additive_terms.add(e, p); }//	��fun��ĩβ���p����enable��ĩβ���0��1
	void add(unsigned e, intermolecular* p)    {    intermolecular_terms.add(e, p); }//	��fun��ĩβ���p����enable��ĩβ���0��1
	void add(unsigned e, conf_independent* p)  {  conf_independent_terms.add(e, p); }//	��fun��ĩβ���p����enable��ĩβ���0��1

	std::vector<std::string> get_names(bool enabled_only) const; // does not include conf-independent ��������
	sz size_internal() const;      //�������ú�������distance_additive_terms�ṹ����fun����Ԫ�ظ�����
									//usable_terms�ṹ����fun����Ԫ�ظ�����
	                                //additive_terms��fun����Ԫ�ظ������ܺ�
	sz size() const { return size_internal() + intermolecular_terms.size(); }//// does not include conf-independent ��������
								   //�ú�������distance_additive_terms�ṹ����fun����Ԫ�ظ�����
								   //usable_terms�ṹ����fun����Ԫ�ظ�����
								   //additive_terms��fun����Ԫ�ظ�����
									//intermolecular_terms������Ԫ�ظ����ܺ�
	sz size_conf_independent(bool enabled_only) const; // number of parameters does not necessarily equal the number of operators
	fl max_r_cutoff() const;             //����
	flv evale(const model& m) const;    //����
	flv evali(const model& m) const;    //����
	flv evale_robust(const model& m) const;//����
	factors eval(const model& m) const;    //���ؽṹ��ĺ���shming
	fl eval_conf_independent(const conf_independent_inputs& in, fl x, flv::const_iterator& it) const;  //����
	flv filter_external(const flv& v) const;  //����
	flv filter_internal(const flv& v) const;  //����
	factors filter(const factors& f) const;   //����
	void display_info() const;      //���ڴ�ӡһЩ��Ҫ����Ϣ
private:
	void eval_additive_aux(const model& m, const atom_index& i, const atom_index& j, fl r, flv& out) const; // out is added to
	//����
};

#endif
