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

#include "quaternion.h"
/*	�ж��Ƿ���Ԫ����׼��
                        input����q����Ԫ����Ҫ������ת��
                        ����ֵΪbool���͵ĺ�����
                        eq(boost::math::abs(q), 1)   ����abs(boost::math::abs(q) - 1)<0.001������1����0��abs(q) Ϊϵ����ƽ����ӿ����ţ�����ʹ�õ���ϵͳ����������Ԫ������ֵ;
                        eq(quaternion_norm_sqr(q), 1)����abs(quaternion_norm_sqr(q) - 1)<0.001������1����0��quaternion_norm_sqr(q)������Ԫ��ϵ����ƽ����ӽ��
                        ���߽������뷵�����ս��
*/
bool quaternion_is_normalized(const qt& q) { // not in the interface, used in assertions
	return eq(quaternion_norm_sqr(q), 1) && eq(boost::math::abs(q), 1);
}


/*	�ж���Ԫ���Ƿ�������
	input����a��b����Ԫ����

	����const qt qt_identity(1, 2, 3, 4);
	R_component_1()//��Ԫ����һ��Ԫ��1;
	R_component_2()//��Ԫ����һ��Ԫ��2;
	R_component_3()//��Ԫ����һ��Ԫ��3;
	R_component_4()//��Ԫ����һ��Ԫ��4;

	eq(a.R_component_1(), b.R_component_1())����������a.R_component_1()-b.R_component_1()�Ĳ�ľ���ֵС��0.001�򷵻�true������false


*/

bool eq(const qt& a, const qt& b) { // elementwise approximate equality - may return false for equivalent rotations
	return eq(a.R_component_1(), b.R_component_1()) && \
		eq(a.R_component_2(), b.R_component_2()) && \
		eq(a.R_component_3(), b.R_component_3()) && \
		eq(a.R_component_4(), b.R_component_4());
}



/*	�Ƕ���Ԫ��
	input����axis��vec�ṹ�壻angle��double��
	�����ж�axis�ṹ����(data[0]^2+data[1]^2+data[2]^2)^(1/2)��1�Ĳ�ľ���ֵ�Ƿ�С��0.001���ǵĻ�������������ֹ���Ƿ��׼����
	normalize_angle(angle)�����Ƕȹ�һ����
	������Ԫ����cos��angle/2����sin��angle/2��*data[0]��sin��angle/2��*data[1]��sin��angle/2��*data[2]��
*/
qt angle_to_quaternion(const vec& axis, fl angle) { // axis is assumed to be a unit vector
													//assert(eq(tvmet::norm2(axis), 1));
	assert(eq(axis.norm(), 1));                    //�ж�axis�ṹ����(data[0]^2+data[1]^2+data[2]^2)^(1/2)��1�Ĳ�ľ���ֵ�Ƿ�С��0.001���ǵĻ�������������ֹ
	normalize_angle(angle); // this is probably only necessary if angles can be very big //�Ƕȹ�һ����ֻ�ڽǶȺܴ��ʱ��ſ�������
	fl c = std::cos(angle / 2);                    //cos��x/2��
	fl s = std::sin(angle / 2);						//sin��x/2��
	return qt(c, s*axis[0], s*axis[1], s*axis[2]);  //������Ԫ����axis[0]=data[0],axis[1]=data[1],axis[2]=data[2]
}



/*   pre_angle_to_quaternion
	input����rotation��vec�ṹ�壻
	inter����angle��double��axis����vec�ṹ�壻
	output������Ԫ����
	angle=(data[0]^2+data[1]^2+data[2]^2)^(1/2)��data������vec�ṹ���Ա����
		if��angle>2.22045e-16��
		�ع�����Ԫ���˷�����axis = ((1/angle)*data[0](1/angle)*data[1],(1/angle)*data[2])��
		Ȼ��ͨ��angle_to_quaternion(axis, angle)�������ؽǶ���Ԫ��
		else
		����ʵ��1
*/
qt angle_to_quaternion(const vec& rotation) {
	//fl angle = tvmet::norm2(rotation); 
	fl angle = rotation.norm();
	if (angle > epsilon_fl) {
		//vec axis; 
		//axis = rotation / angle;	
		vec axis = (1 / angle) * rotation;   //�ع�����Ԫ���˷�����((1/angle)*data[0](1/angle)*data[1],(1/angle)*data[2])
		return angle_to_quaternion(axis, angle);
	}
	return qt_identity;
}



/*
	input����q����Ԫ����
	inter����c��angle��s��double��axis����vec�ṹ�壻

	axis�ṹ���е�data���鴢��q��Ԫ����������ֵ
	�����ж�q�Ƿ��׼����Ȼ��qʵ����ֵ��c��
		if��c��(-1, 1)�ڣ�
			angle=2*acos(c); //������acos
			if��angle>pi��
				angle=angle-2*pi��
			if��abs��sin(angle / 2)��<2.22045e-16��
				���� vec(0, 0, 0)��
			data[0]=data[0]*(angle/sin(angle / 2));
			data[1]=data[1]*(angle/sin(angle / 2));
			data[2]=data[2]*(angle/sin(angle / 2));
			����axis	
	     else
			 ���� vec(0, 0, 0)��
*/

vec quaternion_to_angle(const qt& q) {
	assert(quaternion_is_normalized(q));		//�ж��Ƿ���Ԫ����׼��
	const fl c = q.R_component_1();			    //qʵ������
	if (c > -1 && c < 1) {						// c may in theory be outside [-1, 1] even with approximately normalized q, due to rounding errors
		fl angle = 2 * std::acos(c);		    // acos is in [0, pi]  ������
		if (angle > pi)
			angle -= 2 * pi;					// now angle is in [-pi, pi]
		vec axis(q.R_component_2(), q.R_component_3(), q.R_component_4());  //q����������
		fl s = std::sin(angle / 2);				 // perhaps not very efficient to calculate sin of acos
		if (std::abs(s) < epsilon_fl)
			return zero_vec;
		axis *= (angle / s);                     //�ع�*
		return axis;
	}
	else // when c = -1 or 1, angle/2 = 0 or pi, therefore angle = 0
		return zero_vec;
}



/*
	input����q����Ԫ����
	inter����a��b��c��d��aa��ab��ac��ad��bb��bc��bd��cc��cd��dd��double��
	output����tmp��mat�ṹ�壻
	�����ж�q�Ƿ���Ԫ����׼����Ȼ��a����q��ʵ����b���ڵ�һ������ֵ��c���ڵڶ�������ֵ��d���ڵ���������ֵ��
	aa = a*a;ab = a*b;ac = a*c;ad = a*d;bb = b*b
	bc = b*c;bd = b*d;cc = c*c;cd = c*d;dd = d*d
	�����ж�aa + bb + cc + dd�Ƿ���Ƶ���1���ǵĻ�������������ֹ��
	mat�ṹ���У�
	data[0]=(aa + bb - cc - dd)��data[3]=2 * (-ad + bc)     ��data[6]=2 * (ac + bd)      ;
	data[1]=( 2 * (ad + bc)    ��data[4]=(aa - bb + cc - dd)��data[7]= 2 * (-ab + cd)    ;
	data[2]=2 * (-ac + bd)     ��data[5]=2 * (ab + cd)      ��data[8]=(aa - bb - cc + dd);
*/

mat quaternion_to_r3(const qt& q) {
	assert(quaternion_is_normalized(q));   //�ж�q�Ƿ���Ԫ����׼��

	const fl a = q.R_component_1();
	const fl b = q.R_component_2();
	const fl c = q.R_component_3();
	const fl d = q.R_component_4();

	const fl aa = a*a;
	const fl ab = a*b;
	const fl ac = a*c;
	const fl ad = a*d;
	const fl bb = b*b;
	const fl bc = b*c;
	const fl bd = b*d;
	const fl cc = c*c;
	const fl cd = c*d;
	const fl dd = d*d;

	assert(eq(aa + bb + cc + dd, 1));

	mat tmp;

	// from http://www.boost.org/doc/libs/1_35_0/libs/math/quaternion/TQE.pdf
	tmp(0, 0) = (aa + bb - cc - dd);
	tmp(0, 1) = 2 * (-ad + bc);
	tmp(0, 2) = 2 * (ac + bd);

	tmp(1, 0) = 2 * (ad + bc);
	tmp(1, 1) = (aa - bb + cc - dd);
	tmp(1, 2) = 2 * (-ab + cd);

	tmp(2, 0) = 2 * (-ac + bd);
	tmp(2, 1) = 2 * (ab + cd);
	tmp(2, 2) = (aa - bb - cc + dd);

	return tmp;
}


/*�������
	input����generator�������������
	inter����q����Ԫ����nrm��double��
	random_normal(0, 1, generator)����double���͵�����0��1֮��������̬�ֲ���������������غ���
	q���ĸ�Ԫ��Ϊ�������
	nrmΪq�ľ���ֵ������q��a,b,c,d����qd�ľ���ֵΪ(a^2+b^2+c^2+d^2)^(1/2)
	if��nrm>2.22045e-16)
	q(a/nrm,b/nrm,c/nrm,d/nrm);�����Ʊ�׼������q
	else
	�����������������һ��Ӧ�ò��ᷢ��
*/
qt random_orientation(rng& generator) {
	qt q(random_normal(0, 1, generator),
		random_normal(0, 1, generator),
		random_normal(0, 1, generator),
		random_normal(0, 1, generator));
	fl nrm = boost::math::abs(q);      //�����Ԫ���ľ���ֵ��ϵ��ƽ����ӿ�����
	if (nrm > epsilon_fl) {           //epsilon_fl=2.22045e-16
		q /= nrm;
		assert(quaternion_is_normalized(q));
		return q;
	}
	else
		return random_orientation(generator); // this call should almost never happen�����Ӧ�÷���
}


/*
	input����q����Ԫ����rotation����vec�ṹ��
	�����ж�q�Ƿ��׼����
	�ǵĻ�q����angle_to_quaternion(rotation)���ص���Ԫ����q�ĳ˻�
	��Ԫ���˷���q1��w1,x1,y1,z1��;q2��w2,x2,y2,z2��
	q1 * q2 =
	(w1*w2 - x1*x2 - y1*y2 - z1*z2) +
	(w1*x2 + x1*w2 + y1*z2 - z1*y2) i +
	(w1*y2 - x1*z2 + y1*w2 + z1*x2) j +
	(w1*z2 + x1*y2 - y1*x2 + z1*w2) k

	���q���Ʊ�׼��
*/
void quaternion_increment(qt& q, const vec& rotation) {
	assert(quaternion_is_normalized(q));              //�Ƿ��׼��
	q = angle_to_quaternion(rotation) * q;
	quaternion_normalize_approx(q); // normalization added in 1.1.2
									//quaternion_normalize(q); // normalization added in 1.1.2
}


/*
	input����a��b����Ԫ��;
	inter����tmp����Ԫ����
	�����ж�a��b�Ƿ��׼����Ȼ��b���Ƹ�tmp��
	tmp=tmp/a;
	��Ԫ��������q1/q2����q1��w1,x1,y1,z1��;q2��w2,x2,y2,z2��
	 deno = w2*w2+x2*x2+y2*y2+z2*z2��
 ((w1*w2+x1*x2+y1*y2+z1*z2)/deno,(-w1*x2+x1*w2-y1*z2+z1*y2)/deno,(-w1*y2+x1*z2+y1*w2-z1*x2)/deno,(-w1*z2-x1*y2+y1*x2+z1*w2)/deno)
   ��󷵻�quaternion_to_angle(tmp)��
*/
vec quaternion_difference(const qt& b, const qt& a) { // rotation that needs to be applied to convert a to b
	quaternion_is_normalized(a);
	quaternion_is_normalized(b);
	qt tmp = b;
	tmp /= a; // b = tmp * a    =>   b * inv(a) = tmp 
	return quaternion_to_angle(tmp); // already assert normalization
}


/*
	��ӡ����Ԫ����������qt(a,b,c,d)

*/
void print(const qt& q, std::ostream& out) { // print as an angle
	print(quaternion_to_angle(q), out);
}
