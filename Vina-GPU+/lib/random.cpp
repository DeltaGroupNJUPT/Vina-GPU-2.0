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

#include <ctime> // for time (for seeding)

#include "random.h"
#include "my_pid.h"



/*input  a��b   ��double
	out��random_fl   ��double
function:	��double���͵�����a��b֮�����ɾ��ȷֲ���������������غ���  double����
			a.��֤a��B�����
			b.����bosst�ĺ���  һ�����������������  ��һ�����ƶ�����   ����������������
			c.��֤�����������double���͵������
*/


fl random_fl(fl a, fl b, rng& generator) { // rng�������� �����
	assert(a < b); // ��֤����
	typedef boost::uniform_real<fl> distr;// ָ�����ȷֲ������䣬���������˵㣬ȡdouble����
	boost::variate_generator<rng&, distr> r(generator, distr(a, b));  // �������Ҫ��������������
	fl tmp = r();
	assert(tmp >= a); //��֤�����������������
	assert(tmp <= b);
	return tmp;
}



/*input  mean��sigma   ��double
	out��random_normal   ��double
function:	��double���͵�����mean��sigma֮��������̬�ֲ���������������غ���  double����
			a.��֤sigma��С��0
			b.����bosst�ĺ���  һ�����������������  ��һ�����ƶ�����Ҫ������   ����������������
			c.����double���͵������
*/

fl random_normal(fl mean, fl sigma, rng& generator) { 
	assert(sigma >= 0); // ��̬�ֲ���֤sigma���ڵ���0
	typedef boost::normal_distribution<fl> distr; //������һ��Ĭ�Ϸ��� double �͸���ֵ����̬�ֲ�
	boost::variate_generator<rng&, distr> r(generator, distr(mean, sigma));// �������Ҫ��������������
	return r();
}


/*input  a��b   ��int
	out��random_int   ��int
function:	��int���͵�����a��b֮�����ɾ��ȷֲ�����������������غ���  int����
			a.��֤a��B�����
			b.����bosst�ĺ���  һ�����������������  ��һ�����ƶ�����   ����������������
			c.��֤�����������int���͵������
*/

//  ����һ����ab֮��������  ��������
int random_int(int a, int b, rng& generator) { // rng�������� �����
	assert(a <= b); // ��֤����
	typedef boost::uniform_int<int> distr;  // ָ�����ȷֲ������䣬���������˵㣬ȡ����ֵ
	boost::variate_generator<rng&, distr> r(generator, distr(a, b));   // �������Ҫ��������������
	int tmp = r();
	assert(tmp >= a);    //��֤�����������������
	assert(tmp <= b);
	return tmp;
}
/*input  a��b   ��unsigned int
	out��random_int   ��unsigned int
function:	��unsigned int���͵�����a��b֮�����ɾ��ȷֲ�����������������غ��� unsigned int����
			a.��֤a��B�����  ����a��b��С��0
			b.����bosst�ĺ���  һ�����������������  ��һ�����ƶ�����   ����������������
			c.��֤�����������int���͵������
*/
sz random_sz(sz a, sz b, rng& generator) { // ����һ����ab֮��������  ����
	assert(a <= b);
	assert(int(a) >= 0);
	assert(int(b) >= 0);   //��֤a��B�����  ����a��b��С��0
	int i = random_int(int(a), int(b), generator);    //���������
	assert(i >= 0);
	assert(i >= int(a));
	assert(i <= int(b));    // ��֤�������������
	return static_cast<sz>(i);
}
/*input  true   �� int
	out��tmp����   ��double
function:	���ź�ture������ʱ�򣬷���һ�����飬�������������������0ΪԲ�ģ��뾶Ϊ1��������  double
			a.�������������  double
			b.�������
			c.��֤������
*/
vec random_inside_sphere(rng& generator) {
	while(true) { // on average, this will have to be run about twice
		fl r1 = random_fl(-1, 1, generator);    //����-1����1�����������
		fl r2 = random_fl(-1, 1, generator);
		fl r3 = random_fl(-1, 1, generator);

		vec tmp(r1, r2, r3);//�����������������
		if(sqr(tmp) < 1)//r1��ƽ��+r2��ƽ��+r3��ƽ���ٿ�����<1  ȷ��������
			return tmp;
	}
}
/*input   corner1��corner2    �� double����  ����
	out��tmp����   ��double
function:	�����������飬����һ��������飬Ҫ��������corner1[i]-corner2[i]֮��   double

*/
vec random_in_box(const vec& corner1, const vec& corner2, rng& generator) { // expects corner1[i] < corner2[i]
	vec tmp;  //��������
	VINA_FOR_IN(i, tmp)     //forѭ��  i��0������ĳ���
		tmp[i] = random_fl(corner1[i], corner2[i], generator);//�����������
	return tmp;
}
/*
	out��auto_seed   ��unsigned int
function:	����ʱ��͵�ǰ���̵ı�ʾ��(PID)�ĳ˻�   int���ͣ�������
*/
int auto_seed() { 
	return my_pid() * int(std::time(NULL));
}
