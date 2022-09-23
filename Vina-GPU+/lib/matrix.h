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

#ifndef VINA_MATRIX_H
#define VINA_MATRIX_H

#include <vector>
#include "triangular_matrix_index.h"

// these 4 lines are used 3 times verbatim - defining a temp macro to ease the pain  ���������
#define VINA_MATRIX_DEFINE_OPERATORS \
	const T& operator()(sz i) const { return m_data[i]; } \
	      T& operator()(sz i)       { return m_data[i]; } \
	const T& operator()(sz i, sz j) const { return m_data[index(i, j)]; } \
	      T& operator()(sz i, sz j)       { return m_data[index(i, j)]; } 
/*
����ṹ��matrix 
��������m_data  ����Ϊģ�����T
����uint�ͺ���index   ����  uint  i  j  ���  uint  i + m_i*j
���幹�캯������ֵ


���庯��resize   input  uint  m  n    ģ�����T����  filler_val ֻ��  
function  ���ı�m_data  m_i m_j��ֵ
�ж�
ѭ��ִ��  ij��00 �� m_im_j   �Ѳ���Ϊ��i  j���Ķ���ֵ�� tmp[i+m*j]
��m_data  m_i m_j  ��ֵ


���庯��append  input  ��matrix����ģ�����ΪT�Ĳ���x   ֻ��      ģ�����T���͵Ĳ���filler_val  ���viod
function  �Բ���Ϊ(i+m, j+n)�Ķ���ֵ

����uint���͵ĺ���dim_1dim_2  function ���� m_i  m_j
*/
template<typename T>
class matrix {
	std::vector<T> m_data;//���� m_data
	sz m_i, m_j;
public:
	sz index(sz i, sz j) const {//input   uint  i j   output  uint  i + m_i*j
		assert(j < m_j);
		assert(i < m_i); 
		return i + m_i*j; // column-major
	}
	matrix() : m_i(0), m_j(0) {}
	matrix(sz i, sz j, const T& filler_val) : m_data(i*j, filler_val), m_i(i), m_j(j) {}
	void resize(sz m, sz n, const T& filler_val) { // new sizes should be the same or greater than the old
		if(m == dim_1() && n == dim_2()) return; // �ж�m n   
		VINA_CHECK(m >= dim_1());
		VINA_CHECK(n >= dim_2());
		std::vector<T> tmp(m*n, filler_val);//M*N��Ԫ��  ÿ������filler_val
		VINA_FOR(i, m_i)
			VINA_FOR(j, m_j)
				tmp[i+m*j] = (*this)(i, j);//  �Ѳ���Ϊ��i  j���Ķ���ֵ�� tmp[i+m*j]
		m_data = tmp;			//��ֵm_data 
		m_i = m;
		m_j = n;
	}
	void append(const matrix<T>& x, const T& filler_val) {
		sz m = dim_1();
		sz n = dim_2();
		resize(m + x.dim_1(), n + x.dim_2(), filler_val);//��m_data  m_i m_j  ��ֵ
		VINA_FOR(i, x.dim_1())
			VINA_FOR(j, x.dim_2())(i+m, j+n)
				(*this)(i+m, j+n) = x(i, j);//�Բ���Ϊ(i+m, j+n)�Ķ���ֵ
	}
	VINA_MATRIX_DEFINE_OPERATORS // temp macro defined above
	sz dim_1() const { return m_i; }
	sz dim_2() const { return m_j; }
};
/*������triangular_matrix    ģ�����ΪT
* 
����index  uint����  input  uint i  j  return  i + j*(j+1)/2;ֻ��
����index_permissive  uint����   input  uint  i  j   ֻ��  ����i + j*(j+1)/2;or//j + i*(i+1)/2;
����dim   uint����  ����  uint����   m_dim
*/
template<typename T>
class triangular_matrix {
public:
	std::vector<T> m_data;//����  T����
	sz m_dim;
public:
	sz index(sz i, sz j) const { return triangular_matrix_index(m_dim, i, j); }//i + j*(j+1)/2; 
	sz index_permissive(sz i, sz j) const { return (i < j) ? index(i, j) : index(j, i); }//i + j*(j+1)/2;or//j + i*(i+1)/2;
	triangular_matrix() : m_dim(0) {}
	triangular_matrix(sz n, const T& filler_val) : m_data(n*(n+1)/2, filler_val), m_dim(n) {} 
	VINA_MATRIX_DEFINE_OPERATORS // temp macro defined above     ��
	sz dim() const { return m_dim; }
};
/*������strictly_triangular_matrix  ģ����� T
*
����index    uint����   input  i  j  uint      ֻ��   �ж�ij��ϵȻ�� return   i + j*(j-1)/2

����index_permissive   uint����   input  i  j  uint   ֻ��    return   i + j*(j-1)/2 orreturn   j + i*(i-1)/2

����resize  inout   uint  n    ֻ��T����filler_val  function  ����m_dim��ֵ    ��m_data�������ĳ�*(n-1)/2�����ݵĲ�������Ϊfiller_val

����append 	
input m  strictly_triangular_matrix��     filler_val  T����
out   void
function ��  �����еĲ���  m_dim��ֵ   ����    m_data�������ĳ�*(n-1)/2�����ݵĲ�������Ϊfiller_val
	               Ȼ�� ͨ��������ĵ�ǰ����Ĳ���  

����append  	
input   rectangular  matrix����  ֻ��          triangular     strictly_triangular_matrix����  ֻ��
out    void
function��      1���ж�rectangular��dim_1dim_2��   dim_2==0  ֱ�ӽ���   dim_1=0   ��triangular��ֵ��ǰ����  ����
				2��ִ�к���		append(triangular, filler_val)   *�����еĲ���  m_dim��ֵ   ����    ��ǰm_data�������ĳ�*(n-1)/2�����ݵĲ�������Ϊfiller_val
	                                             ��ͨ������triangular���ĵ�ǰ�����m_data     ����		filler_val=m_data[0]
	            3��Ȼ��  ����rectangular��m_data���ĵ�ǰ�����m_data


*/
template<typename T>
class strictly_triangular_matrix {
	std::vector<T> m_data;  //����   ����ΪT
	sz m_dim;
public:
	sz index(sz i, sz j) const {
		assert(j < m_dim);
		assert(i < j); 
		assert(j >= 1);
		return i + j*(j-1)/2;
	}
	sz index_permissive(sz i, sz j) const { return (i < j) ? index(i, j) : index(j, i); }
	strictly_triangular_matrix() : m_dim(0) {}
	strictly_triangular_matrix(sz n, const T& filler_val) : m_data(n*(n-1)/2, filler_val), m_dim(n) {}
	void resize(sz n, const T& filler_val) {
		if(n == m_dim) return; 
		VINA_CHECK(n > m_dim);
		m_dim = n;
		m_data.resize(n*(n-1)/2, filler_val); // preserves original data//��m_data�������ĳ�*(n-1)/2�����ݵĲ�������Ϊfiller_val
	}


	void append(const strictly_triangular_matrix<T>& m, const T& filler_val) { 
		sz n = dim();//m_dim
		resize(n + m.dim(), filler_val);//��m_dim��ֵ    ��m_data�������ĳ�*(n-1)/2�����ݵĲ�������Ϊfiller_val
		VINA_FOR(i, m.dim())
			VINA_RANGE(j, i+1, m.dim())
				(*this)(i+n, j+n) = m(i, j);//ͨ��������ĵ�ǰ�����m_data  
	}//        m_data[ i+n + m_i*��i+n��]  = m_data[ i + m_i*j]



	void append(const matrix<T>& rectangular, const strictly_triangular_matrix<T>& triangular) {
		VINA_CHECK(dim() == rectangular.dim_1());
		VINA_CHECK(rectangular.dim_2() == triangular.dim());
		if(rectangular.dim_2() == 0) return;//m_i  ==0?
		if(rectangular.dim_1() == 0) {//m_j==0?
			(*this) = triangular;
			return;
		}
		const T& filler_val = rectangular(0, 0); // needed by 'append below'

		sz n = dim();//n=m_dim
		append(triangular, filler_val);
		/*�����еĲ���  m_dim��ֵ   ����    m_data�������ĳ�*(n-1)/2�����ݵĲ�������Ϊfiller_val
	               Ȼ�� ͨ��������ĵ�ǰ�����m_data*/
		VINA_FOR(i, rectangular.dim_1())   //i=0 ;i<m_i;i++
			VINA_FOR(j, rectangular.dim_2())//j=0 ;j<m_j;j++
				(*this)(i, n + j) = rectangular(i, j);//ͨ������rectangular��m_data���ĵ�ǰ�����m_data
		//   m_data[ i + m_i*��J+n��]  = m_data[ i + m_i*j]
	}
	VINA_MATRIX_DEFINE_OPERATORS // temp macro defined above���������
	sz dim() const { return m_dim; }
};

#undef VINA_MATRIX_DEFINE_OPERATORS

#endif
