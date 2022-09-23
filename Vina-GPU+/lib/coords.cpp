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

#include "coords.h"
/*
	input����a��b��vec�ṹ��������
	inter����acc:double;
	demo:ǰ��a��bԪ�ظ���һ�������Ҳ�Ϊ�գ����򷵻�0
	a{��7��8��9������5��6��7��}��b{��1��2��3������3��4��5��} size=2
	acc=��(7-1)^2+(8-2)^2+(9-3)^2��+��(5-3)^2+(6-4)^2+(7-5)^2��=120
	���� ��120/2��^(1/2)=7.74597
*/
fl rmsd_upper_bound(const vecv& a, const vecv& b) {
	VINA_CHECK(a.size() == b.size());
	fl acc = 0;
	VINA_FOR_IN(i, a)       //for(i=0;i<a.size;i++)
		acc += vec_distance_sqr(a[i], b[i]);
	return (a.size() > 0) ? std::sqrt(acc / a.size()) : 0;
}


/*
	input����a��b.coords��vec�ṹ��������b��output_type�ṹ������
	output����tmp��<unint,1.79769e+308>,������first��second����
	demo:ǰ��a��b.coordsԪ�ظ���һ�������Ҳ�Ϊ�գ����򷵻�res=0
	 a{��7��8��9������5��6��7��}��b.coords{��1��2��3������3��4��5��} size=2
	res=[(��(7-1)^2+(8-2)^2+(9-3)^2��+��(5-3)^2+(6-4)^2+(7-5)^2��)/2]^(1/2)=7.74597
	����tmp<2,7.74597>
*/

std::pair<sz, fl> find_closest(const vecv& a, const output_container& b) {
	std::pair<sz, fl> tmp(b.size(), max_fl);          //max_fl=1.79769e+308
	VINA_FOR_IN(i, b) { 
		fl res = rmsd_upper_bound(a, b[i].coords);
		if(i == 0 || res < tmp.second)
			tmp = std::pair<sz, fl>(i, res);
	}
	return tmp;
}


/*
	input����out��output_type�ṹ��������t��output_type�ṹ�壬min_rmsd��double��max_size��unint
	����find_closest������demo��closest_rmsd<2,7.74597>
	�����out����Ԫ��>2��min_rmsd>7.74597����ʾ��һ���ǳ�С��ֵ������t.e < out[2].e
		out[2]=t;
	�����򣨱�ʾû��С��ֵ��
	���out����Ԫ��<max_size,��out����ĩβ����t��output_type���ṹ��
	�������out������Ϊ�ղ���t.e < outĩβԪ�ص�eֵ����outĩβԪ�ظ�ֵΪt
	�����out������out�ṹ����e������������

*/
void add_to_output_container(output_container& out, const output_type& t, fl min_rmsd, sz max_size) {
	std::pair<sz, fl> closest_rmsd = find_closest(t.coords, out);
	if(closest_rmsd.first < out.size() && closest_rmsd.second < min_rmsd) { // have a very similar one
		if(t.e < out[closest_rmsd.first].e) { // the new one is better, apparently
			out[closest_rmsd.first] = t; // FIXME? slow
		}
	}
	else { // nothing similar
		if(out.size() < max_size)   //new���ٵĿռ��ڶ��ϣ���һ�������ı��������ջ��
			out.push_back(new output_type(t)); // the last one had the worst energy - replacing 
		else
			if(!out.empty() && t.e < out.back().e) // FIXME? - just changed
				out.back() = t; // FIXME? slow
	}
	out.sort();
}
