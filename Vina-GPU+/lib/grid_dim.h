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

#ifndef VINA_GRID_DIM_H
#define VINA_GRID_DIM_H

#include <boost/array.hpp>

#include "common.h"


/*�ṹ�壺
* double�� begin��end    Uint��n
* function1������ double��begin��end    Uint��n  ��ֵ0
* function2������ end- begin��ֵ  double  ֻ��
* function3�� �ж�n>0    bool  ֻ��
* function4: ���룺�����ṹ�� �����bool  �ж������ṹ���n�Ƿ�һ�� ����begin��end�Ĳ�ľ���ֵ�Ƿ�<0.01
* function5�����룺�����ṹ���飨�����ṹ����ɣ� �����bool  �жϽṹ�����е����������е�a[i]��b[i]�Ĳ�ľ���ֵ�Ƿ�<0.001
* function6�����룺�ṹ���飨�����ṹ����ɣ�  std::ostream�������ã�Ĭ��ֵΪstd::cout   ���ܣ���ӡ
* function7�����룺�ṹ���飨�����ṹ����ɣ�  ��� ������3���� double  ����������begin��������
* function8�����룺�ṹ���飨�����ṹ����ɣ�   ��� ������3���� double   ����������end��������
* ����ʽ���л�     private����
*/
struct grid_dim {
	fl begin;
	fl end;
	sz n; // number of intervals == number of sample points - 1
	grid_dim() : begin(0), end(0), n(0) {}
	fl span() const { return end - begin; }
	bool enabled() const { return (n > 0); }
private:                                           //����ʽ���л�
	friend class boost::serialization::access;
	template<class Archive>
	void serialize(Archive& ar, const unsigned version) {
		ar & begin;
		ar & end;
		ar & n;
	}
};

inline bool eq(const grid_dim& a, const grid_dim& b) {
	return a.n == b.n && eq(a.begin, b.begin) && eq(a.end, b.end);//�ж������ṹ���n�Ƿ�һ�� ����begin��end�Ĳ�ľ���ֵ�Ƿ�<0.01
}

typedef boost::array<grid_dim, 3> grid_dims;

inline bool eq(const grid_dims& a, const grid_dims& b) {
	return eq(a[0], b[0]) && eq(a[1], b[1]) && eq(a[2], b[2]);//�жϽṹ�����е����������е�a[i]��b[i]�Ĳ�ľ���ֵ�Ƿ�<0.001
}       

inline void print(const grid_dims& gd, std::ostream& out = std::cout) {
	VINA_FOR_IN(i, gd)
		std::cout << gd[i].n << " [" << gd[i].begin << " .. " << gd[i].end << "]\n";
}//    ��ӡn[b..e]   

inline vec grid_dims_begin(const grid_dims& gd) {
	vec tmp;
	VINA_FOR_IN(i, gd)
		tmp[i] = gd[i].begin;
	return tmp;
}//������begin��������

inline vec grid_dims_end(const grid_dims& gd) {
	vec tmp;
	VINA_FOR_IN(i, gd)
		tmp[i] = gd[i].end;
	return tmp;
}//������end��������

#endif
