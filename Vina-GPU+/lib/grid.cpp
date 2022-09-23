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

#include "grid.h"
/*input �� �����ṹ����ɵĽṹ����
function�����ṹ���е�������鸳ֵ
output��void    �ı�������ֵ

*/
void grid::init(const grid_dims& gd) {//����Ϊ�����ṹ����ɵĽṹ����
	m_data.resize(gd[0].n+1, gd[1].n+1, gd[2].n+1);   //�����ṹ���е���n+1Ȼ���   ������������������
	m_init = vec(gd[0].begin, gd[1].begin, gd[2].begin);//����grid�е����鸳ֵ     ÿ���ṹ���begin
	m_range = vec(gd[0].span(), gd[1].span(), gd[2].span());//����grid�е����鸳ֵ  ÿ���ṹ���end - begin
	assert(m_range[0] > 0);//ȷ����������>0
	assert(m_range[1] > 0);
	assert(m_range[2] > 0);
	m_dim_fl_minus_1 = vec(m_data.dim0() - 1.0, 
	                       m_data.dim1() - 1.0,
			               m_data.dim2() - 1.0);//��������������һ   ����grid�е����鸳ֵ
	VINA_FOR(i, 3) {
		m_factor[i] = m_dim_fl_minus_1[i] / m_range[i];   //m_dim_fl_minus_1/m_range    ����grid�е����鸳ֵ
		m_factor_inv[i] = 1 / m_factor[i];//1/m_factor     ����grid�е����鸳ֵ
	}
}// 
/*input ������3double��������location  ֻ��     double��  slope��v   ָ������deriv/NULL
output��  double
inter������3double��������s   miss  ����3int��������region  ����3uint��������a
function�����double���͵������Ҹı�deriv��Ӧ��ֵ  ���ļ�    //��Ϊָ����Ч����Ч����  

*/
fl grid::evaluate_aux(const vec& location, fl slope, fl v, vec* deriv) const { // sets *deriv if not NULL
	vec s  = elementwise_product(location - m_init, m_factor); //a0*b0    a1*b1  a2*b2   ��������1   �������е���������������s��ֵ

	vec miss(0, 0, 0);    //��������2
	boost::array<int, 3> region;  //��������3
	boost::array<sz, 3> a;  //��������4

	VINA_FOR(i, 3) {
		if(s[i] < 0) {    //��������1   ������1234��ֵ
			miss[i] = -s[i];
			region[i] = -1;
			a[i] = 0; 
			s[i] = 0;
		}       
		else if(s[i] >= m_dim_fl_minus_1[i]) {//��������1   ������1234��ֵ
			miss[i] = s[i] - m_dim_fl_minus_1[i];
			region[i] = 1;
			assert(m_data.dim(i) >= 2);
			a[i] = m_data.dim(i) -  2; 
			s[i] = 1;
		}
		else {
			region[i] = 0; // now that region is boost::array, it's not initialized
			a[i] = sz(s[i]); //��������1   ������1234��ֵ
			s[i] -= a[i];
		}
		assert(s[i] >= 0);
		assert(s[i] <= 1);
		assert(a[i] >= 0);
		assert(a[i]+1 < m_data.dim(i));//�ж�s[i]��0-1֮��  a[i]>=0��a[0]<m_i a[1]<m_j a[2]<m_k
	}
	const fl penalty = slope * (miss * m_factor_inv); // �����һ����

	const sz x0 = a[0];
	const sz y0 = a[1];
	const sz z0 = a[2];

	const sz x1 = x0+1;
	const sz y1 = y0+1;
	const sz z1 = z0+1;


	const fl f000 = m_data(x0, y0, z0);
	const fl f100 = m_data(x1, y0, z0);
	const fl f010 = m_data(x0, y1, z0);
	const fl f110 = m_data(x1, y1, z0);
	const fl f001 = m_data(x0, y0, z1);
	const fl f101 = m_data(x1, y0, z1);
	const fl f011 = m_data(x0, y1, z1);
	const fl f111 = m_data(x1, y1, z1);

	const fl x = s[0];
	const fl y = s[1];
	const fl z = s[2];

	const fl mx = 1-x;
	const fl my = 1-y;
	const fl mz = 1-z;

	fl f = 
		f000 *  mx * my * mz  +
		f100 *   x * my * mz  +
		f010 *  mx *  y * mz  + 
		f110 *   x *  y * mz  +
		f001 *  mx * my *  z  +
		f101 *   x * my *  z  +
		f011 *  mx *  y *  z  +
		f111 *   x *  y *  z  ;

	if(deriv) { // valid pointer  ��ָ�����Ͳ���
		const fl x_g = 
			f000 * (-1)* my * mz  +
			f100 *   1 * my * mz  +
			f010 * (-1)*  y * mz  + 
			f110 *   1 *  y * mz  +
			f001 * (-1)* my *  z  +
			f101 *   1 * my *  z  +
			f011 * (-1)*  y *  z  +
			f111 *   1 *  y *  z  ;


		const fl y_g = 
			f000 *  mx *(-1)* mz  +
			f100 *   x *(-1)* mz  +
			f010 *  mx *  1 * mz  + 
			f110 *   x *  1 * mz  +
			f001 *  mx *(-1)*  z  +
			f101 *   x *(-1)*  z  +
			f011 *  mx *  1 *  z  +
			f111 *   x *  1 *  z  ;


		const fl z_g =  
			f000 *  mx * my *(-1) +
			f100 *   x * my *(-1) +
			f010 *  mx *  y *(-1) + 
			f110 *   x *  y *(-1) +
			f001 *  mx * my *  1  +
			f101 *   x * my *  1  +
			f011 *  mx *  y *  1  +
			f111 *   x *  y *  1  ;

		vec gradient(x_g, y_g, z_g);
		curl(f, gradient, v);//�ı�f��gradient
		vec gradient_everywhere;

		VINA_FOR(i, 3) {//�ı�deriv��Ӧ��ֵ
			gradient_everywhere[i] = ((region[i] == 0) ? gradient[i] : 0);
			(*deriv)[i] = m_factor[i] * gradient_everywhere[i] + slope * region[i];
		}

		return f + penalty;//�������
	}
	else {  //������ָ�����Ͳ���
		curl(f, v);
		return f + penalty;
	}
} 
