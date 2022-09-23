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

#include "mutate.h"

/*
* ������count_mutable_entities
* input:conf�����
* output��unsigned int����counter
* counter��СΪ��2��conf�����ĳ�Աligands�Ĵ�С+ligands������ÿ��Ԫ�أ�ligand_conf�ࣩ�ĳ�Աtorsions�����Ĵ�С֮��+flex������ÿ��Ԫ�أ�residue_conf�ࣩ�ĳ�Աtorsions�����Ĵ�С֮��
* ligands������flex��������conf������ݳ�Ա
*/
sz count_mutable_entities(const conf& c) {
	//unsigned int���ͼ���ֵ
	sz counter = 0;

	//��������ligands����������Ԫ��Ϊligand_conf�ṹ��
	//ѭ����������ִ�����²�����1.��������ligands�ĵ�i��Ԫ�أ������ݳ�Աtorsions������doubleԪ�ص��������Ĵ�С
	//2.counter=counter+��1���õ��Ĵ�С+2

	VINA_FOR_IN(i, c.ligands)
		counter += 2 + c.ligands[i].torsions.size();

	//��������flex����������Ԫ��Ϊresidue_conf�ṹ��
    //ѭ����������ִ�����²�����1.��������flex�ĵ�i��Ԫ�أ������ݳ�Աtorsions������doubleԪ�ص��������Ĵ�С
    //2.counter=counter+��1���õ��Ĵ�С

	VINA_FOR_IN(i, c.flex)
		counter += c.flex[i].torsions.size();
	//����counter
	return counter;
}


/*
* ������mutate_conf
* input:conf�����model����󣬸�����amplitude��rng��
* output:��
* ����������conf������һ�����ݳ�Ա��ֵ�����帳ֵ��������ϸ����
* ligands������flex��������conf������ݳ�Ա
*/
// does not set model
void mutate_conf(conf& c, const model& m, fl amplitude, rng& generator) { // ONE OF: 2A for position, similar amp for orientation, randomize torsion
	//����conf�����Ŀɱ�ʵ���������ظ�mutable_entities_num
	sz mutable_entities_num = count_mutable_entities(c);
	//�ɱ�ʵ����Ϊ0�����˳��������������ִ������
	if(mutable_entities_num == 0) return;
	//����һ��0��mutable_entities_num - 1֮������������which_int������������Ͼ��ȷֲ�
	int which_int = random_int(0, int(mutable_entities_num - 1), generator);
	//�ж����������which_int�Ƿ���ڵ���0�����򱨴���ֹ����ִ��
	VINA_CHECK(which_int >= 0);
	//��������which_intת����unsigned int����which
	sz which = sz(which_int);
	//�ж�which�Ƿ�С��mutable_entities_num�����򱨴���ֹ����ִ��
	VINA_CHECK(which < mutable_entities_num);

	//��������ligands����������Ԫ��Ϊligand_conf�ṹ��
	//ѭ����������ִ�����²�����
	//1.���which����0��
	//������ligands�ĵ�i��Ԫ�أ������ݳ�Աrigid��rigid_conf�ࣩ�����ݳ�Աposition��vecʸ����������double�������ݣ�ִ�й�ʽ1����ֱ���˳�����
	//������е�2��

	//2.���which������0��which�Ƚ����Լ�1�����ж�which�Ƿ����0��
	//�������0������һ��������gr������ѭ������i��ֵ��������㹫ʽ��model�ļ�gyration_radius����)
	//�ٱȽ�gr�Ƿ�������б������ļ��������ʶ�����С���㸡�������������
	//��1������һ��������double����Ԫ�ص�ʸ��rotation����ʼ��rotation����ʽ2��,�ٶ�����ligands�ĵ�i��Ԫ�أ������ݳ�Աrigid��rigid_conf�ࣩ�����ݳ�Աorientation����rotationʸ������quaternion_increment���¸���ֵ��������㹫ʽ��quaternion.cpp�ļ���
	//���ֱ���˳�������������е�3��������������2����ͷ��������ף�




	//3.which�Ƚ����Լ�1�����ж�which�Ƿ�С��{����ligands�ĵ�i��Ԫ�أ������ݳ�Աtorsions������doubleԪ�ص��������Ĵ�С},���С�ڣ�
	//������ligands�ĵ�i��Ԫ�أ������ݳ�Աtorsions������doubleԪ�ص���������which��Ԫ�ص���
	//-3.1415926535897931��3.1415926535897931֮���double���������������������Ͼ��ȷֲ�����ֱ���˳�������
	//������which����which-{����ligands�ĵ�i��Ԫ�أ������ݳ�Աtorsions������doubleԪ�ص��������Ĵ�С}
	//ligands������flex��������conf������ݳ�Ա


	VINA_FOR_IN(i, c.ligands) {
		if(which == 0) { c.ligands[i].rigid.position += amplitude * random_inside_sphere(generator); return; }
		--which;
		if(which == 0) { 
			fl gr = m.gyration_radius(i); 
			if(gr > epsilon_fl) { // FIXME? just doing nothing for 0-radius molecules. do some other mutation?
				vec rotation; 
				rotation = amplitude / gr * random_inside_sphere(generator); 
				quaternion_increment(c.ligands[i].rigid.orientation, rotation);
			}
			return; 
		}
		--which;
		if(which < c.ligands[i].torsions.size()) { c.ligands[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= c.ligands[i].torsions.size();
	}

	//��������flex����������Ԫ��Ϊresidue_conf�ṹ��
	//ѭ����������ִ�����²�����1.��������flex�ĵ�i��Ԫ�أ������ݳ�Աtorsions������doubleԪ�ص��������Ĵ�С
	
	//2.���whichС�ڵ�1���ҵ�������������flex�ĵ�i��Ԫ�أ������ݳ�Աtorsions������doubleԪ�ص���������which��Ԫ�ص���
	//-3.1415926535897931��3.1415926535897931֮���double���������������������Ͼ��ȷֲ�����ֱ�ӽ����˳�����
	
	//3.���whichС�ڵ�1���ҵ���������which����which-��1���ҵ�����

	VINA_FOR_IN(i, c.flex) {
		if(which < c.flex[i].torsions.size()) { c.flex[i].torsions[which] = random_fl(-pi, pi, generator); return; }
		which -= c.flex[i].torsions.size();
	}
}
