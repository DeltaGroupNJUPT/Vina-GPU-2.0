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

#include "weighted_terms.h"


/*
	input����t��terms�ṹ�塢weights��double������double����cutoff_��ʼ��Ϊ0�����봫��weighted_terms�ṹ���У�
	inter����enabled_usable_terms��unint���� ��atom_typing_used_��atom_type::tö�����ͣ�
		�����ж�ģ�����Ϊdistance_additive��term_set�ṹ����enabled�д���0��ֵ�ĸ�����ģ�����Ϊadditive��term_set�ṹ����enabled�д���0��ֵ�ĸ���
	��ģ�����Ϊintermolecular��term_set�ṹ����enabled�д���0��ֵ�ĸ��� �Ƿ񶼵���0���ǵĻ��������򣬷����˳���
		if��ģ�����Ϊusable��term_set�ṹ����enabled[i]����0��(i=0-ģ�����Ϊusable��term_set�ṹ����fun����Ԫ�ظ�����
			if��enabled_usable_terms����Ϊ�գ�{
				atom_typing_used_����usable_terms.sizeΪģ�����Ϊusable��fun[i]��atom_typing_usedö���ͣ�
		    else
			    �ж�atom_typing_used_�Ƿ����usable_terms.sizeΪģ�����Ϊusable��fun[i]��atom_typing_usedö���ͣ��ǵĻ�����������ֹ����
			��enabled_usable_terms�����������i��ţ�
			cutoff_��ֵΪģ�����Ϊusable��fun[i]��cutoff
			}
*/
weighted_terms::weighted_terms(const terms* t, const flv& weights) : t(t), weights(weights), cutoff_(0) { // does not own t
	VINA_CHECK(t->distance_additive_terms.num_enabled() == 0); //�ú�������ģ�����Ϊdistance_additive��term_set�ṹ����enabled�д���0��ֵ�ĸ���
	VINA_CHECK(t->         additive_terms.num_enabled() == 0); //�ú�������ģ�����Ϊadditive��term_set�ṹ����enabled�д���0��ֵ�ĸ���
	VINA_CHECK(t->   intermolecular_terms.num_enabled() == 0); //�ú�������ģ�����Ϊintermolecular��term_set�ṹ����enabled�д���0��ֵ�ĸ���      
	VINA_FOR_IN(i, t->usable_terms)                            //for(i=0;i<usable_terms.size;i++)��usable_terms.sizeΪģ�����Ϊusable��term_set�ṹ����fun����Ԫ�ظ���
		if(t->usable_terms.enabled[i]) {
			if(enabled_usable_terms.empty())					//enabled_usable_termsΪunint����      
				atom_typing_used_ = t->usable_terms[i].atom_typing_used;//atom_typing_used_��atom_type::tö�����ͣ�usable_terms[i]=fun[i]
			else
				VINA_CHECK(atom_typing_used_ == t->usable_terms[i].atom_typing_used);

			enabled_usable_terms.push_back(i);
			cutoff_ = (std::max)(cutoff_, t->usable_terms[i].cutoff);
		}
}

/*
	input����t1��t2:unint,   r:double;
	inter����acc��double��enabled_usable_terms��unint������weights��double������
	for(i=0;i<enabled_usable_terms.size;i++)
   ��weight[i]�� t->usable_terms[enabled_usable_terms[i]].eval(t1, t2, r) �ĳ˻������ۼ�
   ������һ������t->usable_terms[enabled_usable_terms[i]].eval(t1, t2, r)�Ƿ���0��������������������������������
*/

fl weighted_terms::eval(sz t1, sz t2, fl r) const { // intentionally not checking for cutoff
	fl acc = 0;
	VINA_FOR_IN(i, enabled_usable_terms)                  //for(i=0;i<enabled_usable_terms.size;i++)
		acc += weights[i] * t->usable_terms[enabled_usable_terms[i]].eval(t1, t2, r);     // t->usable_terms[enabled_usable_terms[i]].eval(t1, t2, r)����0��������
	return acc;
}




/*
	input����m:model�ṹ��,e:double;
	inter����enabled_usable_terms��unint������weights��double������it����������in������Ϊmodel�ṹ���conf_independent_inputs�ṹ�壻
	output����tmp��dopuble
	it������ָ��weights�����е�enabled_usable_terms����Ԫ�ظ�����Ԫ�أ�
	tmp����e + read_iterator(i) * in.******   // read_iterator(i)Ϊ������ָ����һ��Ԫ�� ��******Ӧ����ȡ���ڵ������������ǰһ��everthing.cpp�еĺ���������
	����ж�it�Ƿ�ָ��weights���������һλ���ǵĻ�����tmp
*/

fl weighted_terms::conf_independent(const model& m, fl e) const {
	flv::const_iterator it = weights.begin() + enabled_usable_terms.size();     //ָ��weights�����е�enabled_usable_terms����Ԫ�ظ�����Ԫ��
	conf_independent_inputs in(m); // FIXME quite inefficient, but I think speed is irrelevant here, right?
	fl tmp = t->eval_conf_independent(in, e, it);
	assert(it == weights.end());        //�ж�it�Ƿ�ָ��weights���������һλ
	return tmp;
}
