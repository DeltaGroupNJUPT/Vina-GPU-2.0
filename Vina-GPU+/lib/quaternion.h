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

#ifndef VINA_QUATERNION_H
#define VINA_QUATERNION_H

#include <boost/math/quaternion.hpp>
#include <boost/serialization/split_free.hpp>

#include "common.h"
#include "random.h"

typedef boost::math::quaternion<fl> qt;              //��Ԫ��
bool quaternion_is_normalized(const qt & q);         //����bool���͵�quaternion_is_normalized���������뺯��Ϊ��Ԫ��

// non-intrusive free function split serialization


/*
		save������input����qt:��Ԫ��
		����Ԫ����һ��ʵ�����������洢���ļ��л���ͨ�����緢��ȥ��

		load���������ļ��е���Ԫ����һ��ʵ����������д��qt��Ԫ����
		���л�����save�����������л�����load������
		����const qt qt_identity(1, 2, 3, 4);
		R_component_1()//��Ԫ����һ��Ԫ��1;
		R_component_2()//��Ԫ����һ��Ԫ��2;
		R_component_3()//��Ԫ����һ��Ԫ��3;
		R_component_4()//��Ԫ����һ��Ԫ��4;


*/
namespace boost {
	namespace serialization {
		template<class Archive>
		void save(Archive& ar, const qt& q, const unsigned version) {//����const qt qt_identity(1, 2, 3, 4);
			fl q1 = q.R_component_1();						//q1:��Ԫ����һ��Ԫ��1;
			fl q2 = q.R_component_2();						//q2:��Ԫ����һ��Ԫ��2;
			fl q3 = q.R_component_3();						//q3:��Ԫ����һ��Ԫ��3;
			fl q4 = q.R_component_4();						//q4:��Ԫ����һ��Ԫ��4;

			ar & q1;
			ar & q2;
			ar & q3;
			ar & q4;
		}


		template<typename Archive>
		void load(Archive& ar, qt& q, const unsigned version) {
			fl a, b, c, d;
			ar & a;
			ar & b;
			ar & c;
			ar & d;
			q = qt(a, b, c, d);
		}
	}
}


BOOST_SERIALIZATION_SPLIT_FREE(qt)//BOOST_SERIALIZATION_SPLIT_FREE()������Boost���л���ʹ��save��load����serialize����




/*
	��������
*/
  
bool eq(const qt& a, const qt& b); // elementwise approximate equality - may return false for equivalent rotations
const qt qt_identity(1, 0, 0, 0);   //������һ��ʵ��Ϊ1���鲿��Ϊ0����Ԫ��
qt angle_to_quaternion(const vec& axis, fl angle); // axis is assumed to be a unit vector
qt angle_to_quaternion(const vec& rotation); // rotation == angle * axis
vec quaternion_to_angle(const qt& q);
mat quaternion_to_r3(const qt& q);

/*    ��Ԫ������ֵ��ƽ��

		input����qt:��Ԫ��
		��qt=a+bi+cj+dk;
		�������أ�a^2+b^2+c^2+d^2

*/
inline fl quaternion_norm_sqr(const qt& q) { // equivalent to sqr(boost::math::abs(const qt&))
	return sqr(q.R_component_1()) + sqr(q.R_component_2()) + sqr(q.R_component_3()) + sqr(q.R_component_4());
}



/*	��Ԫ����׼��

	input����q����Ԫ��
	���������Ҫͨ��quaternion_norm_sqr����������õ���Ԫ��ϵ����ƽ����Ӻ��ڽ��俪��������Ԫ������ֵ
	����õ�ֵ���ڼ������ʶ��ĸ������ɱ�ʾ����Сֵ2.22045e-16����
		q=q*��1/a��
	����ж�bool(abs(boost::math::abs(q) - 1)<0.001)����bool(abs(quaternion_norm_sqr(q) - 1)<0.001)��Ϊ��������������ֹ����
*/

inline void quaternion_normalize(qt& q) {
	const fl s = quaternion_norm_sqr(q);
	assert(eq(s, sqr(boost::math::abs(q))));  //����Լ�д��ϵ��ƽ������ϵͳ�ĺ����������Ԫ������ֵƽ���Ĳ�ľ���ֵС��0.001�򷵻�true������false   
    const fl a = std::sqrt(s);                //		����Ԫ������ֵ
	assert(a > epsilon_fl);                   //epsilon_fl�������ʶ��ĸ������ɱ�ʾ����Сֵ���ҵļ������ʶ�����С��������2.22045e-16
	q *= 1/a;                                 //
	assert(quaternion_is_normalized(q));
}


/*��Ԫ�����Ʊ�׼��
		input����q����Ԫ����tolerance��double��1e-6��
		��ϵ����ƽ���ͼ�1�ľ���ֵС��1e-6�����޲��������ֳ�ʼֵ�������ж�ϵ����ƽ���͵Ŀ������Ƿ���ڼ������ʶ��ĸ������ɱ�ʾ����Сֵ2.22045e-16���ǵĻ���
			q=q*��1/a��
	����ж�bool(abs(boost::math::abs(q) - 1)<0.001)����bool(abs(quaternion_norm_sqr(q) - 1)<0.001)��Ϊ��������������ֹ����


*/
inline void quaternion_normalize_approx(qt& q, const fl tolerance = 1e-6) {
	const fl s = quaternion_norm_sqr(q);      //ϵ����ƽ����
	assert(eq(s, sqr(boost::math::abs(q))));  //����Լ�д��ϵ��ƽ������ϵͳ�ĺ����������Ԫ������ֵƽ���Ĳ�ľ���ֵС��0.001�򷵻�true������false   
    if(std::abs(s - 1) < tolerance)            //ϵ����ƽ���ͼ�һ�ľ���ֵС��1e-6
        ; // most likely scenario
    else {
        const fl a = std::sqrt(s);             //����Ԫ������ֵ
        assert(a > epsilon_fl);
        q *= 1/a;
        assert(quaternion_is_normalized(q));
    }
}

/*
��������
*/

qt random_orientation(rng& generator);
void quaternion_increment(qt& q, const vec& rotation);
vec quaternion_difference(const qt& b, const qt& a); // rotation that needs to be applied to convert a to b
void print(const qt& q, std::ostream& out = std::cout); // print as an angle

#endif
