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

#ifndef VINA_STATISTICS_H
#define VINA_STATISTICS_H

#include <algorithm> // sort
#include <cmath> // sqrt
#include "common.h" // for flv

// simple but inefficient implementations


/*   ���ֵ
	input���� v:double����, ,ֻ����
	inter����acc:double;
	output����mean:double;
	v.size()����v����Ԫ�ظ�����
	v.empty()�����ж�v�����Ƿ�Ϊ�գ�Ϊ�շ���1����Ϊ0��
	        acc=v[i]+acc;i=0-v.size
			if(v����Ϊ��)
			mean=0;
			else 
			mean=acc/v.size;
*/

inline fl mean(const flv& v) {
	fl acc = 0;
	VINA_FOR_IN(i, v)						  //for(uint VINA_MACROS_TMP=v.size,i=0; i<VINA_MACROS_TMP;i++)
		acc += v[i];						  //�������е�ֵ�����ۼ�
	return v.empty() ? 0 : (acc/v.size());    //��vΪ����������0�����򷵻��ۼӵ�ֵ������Ԫ�ظ���
}




/*   ��ƫ��
input���� v:double����, ,ֻ����
inter����acc:double��m:double;
output����deviation:double;
v.size()����v����Ԫ�ظ�����
v.empty()�����ж�v�����Ƿ�Ϊ�գ�Ϊ�շ���1����Ϊ0��
m = mean(v)������v�����ľ�ֵ
             acc=(v[i]-mean(m))^2+acc;i=0-v.size
             if(v����Ϊ��)
             deviation=0;
             else
             deviation=(acc/v.size)^(1/2);
*/

inline fl deviation(const flv& v) {
	fl m = mean(v);
	fl acc = 0;
	VINA_FOR_IN(i, v)								//for(uint VINA_MACROS_TMP=v.size,i=0; i<VINA_MACROS_TMP;i++)
		acc += sqr(v[i] - m);
	return v.empty() ? 0 : std::sqrt(acc/v.size());
}


/*   ��������
input���� a��b:double���� ,ֻ����
inter����acc:double��m:double;
output����rmsd:double;
a.size()����a����Ԫ�ظ�����
a.empty()�����ж�a�����Ƿ�Ϊ�գ�Ϊ�շ���1����Ϊ0��
����ǰҪ���ж�a��b����Ԫ�ظ����Ƿ���ȣ���Ȳ��ܽ��н����������㡣�����˳���
              acc=(a[i]-b[i])^2+acc;i=0-a.size
              if(v����Ϊ��)
              rmsd=0;
              else
              rmsd=(acc/a.size)^(1/2);
*/

inline fl rmsd(const flv& a, const flv& b) {
	VINA_CHECK(a.size() == b.size());                 //���� assert(a.size() == b.size())����a��b����������ͬ����ִ�У�������ֹ����
	fl acc = 0;
	VINA_FOR_IN(i, a)                                //for(uint VINA_MACROS_TMP=v.size,i=0; i<VINA_MACROS_TMP;i++)
		acc += sqr(a[i] - b[i]);
	return a.empty() ? 0 : std::sqrt(acc/a.size());
}



/*  
input���� a��b:double���� ,ֻ����
inter����acc:double��m:double;
output����average_difference:double;
a.size()����a����Ԫ�ظ�����
a.empty()�����ж�a�����Ƿ�Ϊ�գ�Ϊ�շ���1����Ϊ0��
����ǰҪ���ж�a��b����Ԫ�ظ����Ƿ���ȣ���Ȳ��ܽ��н����������㡣�����˳���
           acc=(b[i]-a[i])+acc;i=0-a.size
           if(v����Ϊ��)
           average_difference=0;
           else
           average_difference=(acc/a.size);
*/

inline fl average_difference(const flv& b, const flv& a) { // b - a
	VINA_CHECK(a.size() == b.size());                      //���� assert(a.size() == b.size())����a��b����������ͬ����ִ�У�������ֹ����
	fl acc = 0;
	VINA_FOR_IN(i, a)									   //for(uint VINA_MACROS_TMP=v.size,i=0; i<VINA_MACROS_TMP;i++)
		acc += b[i] - a[i];
	return a.empty() ? 0 : (acc/a.size());
}





/*   
input���� x��y:double����,ֻ����
inter����sum_x,sum_y,sum_x_sq,sum_y_sq,sum_prod��sd_x��sd_y��cov��tmp:double;
output����pearson:double;
x.size()����a����Ԫ�ظ�����
x.empty()�����ж�a�����Ƿ�Ϊ�գ�Ϊ�շ���1����Ϊ0��
epsilon_flΪ�������ʶ��ĸ������ɱ�ʾ����Сֵ��2.22045e-16��
����ǰҪ���ж�a��b����Ԫ�ظ����Ƿ���ȣ���Ȳ��ܽ��н����������㡣�����˳���
������Ϊ��ֱ�ӷ���0��
                     sum_x=sum_x+x[i];i=0-n(a.size)
					 sum_y=sum_y+y[i];
					 sum_x_sq=sum_x_sq+x[i]^2;
					 sum_y_sq=sum_y_sq+y[i]^2;
					 sum_prod=sum_prod+y[i]*x[i];
					 sd_x=(sum_x_sq/n - (sum_x/n)^2)^(1/2);//��׼��
					 sd_y=(sum_y_sq/n - (sum_y/n)^2)^(1/2);//��׼��
					 cov= sum_prod/n - (sum_x/n) * (sum_y/n);//Э����
					 tmp=sd_x * sd_y;
					 if��abs(tmp) < epsilon_fl��
					 pearson= 0;
					 else 
					 pearson= cov / tmp;

*/

inline fl pearson(const flv& x, const flv& y) {
	sz n = x.size();
	VINA_CHECK(n == y.size());                          //���� assert(a.size() == b.size())����a��b����������ͬ����ִ�У�������ֹ����
	if(n == 0) return 0; 
	fl sum_x = 0;
	fl sum_y = 0;
	fl sum_x_sq = 0;
	fl sum_y_sq = 0;
	fl sum_prod = 0;
	
	VINA_FOR(i, n) {                                  //for(uint VINA_MACROS_TMP=n,i=0; i<VINA_MACROS_TMP;i++)
		sum_x += x[i];
		sum_y += y[i];
		sum_x_sq += sqr(x[i]);
		sum_y_sq += sqr(y[i]);
		sum_prod += x[i] * y[i];
	}
	fl sd_x = std::sqrt(sum_x_sq/n - sqr(sum_x/n));  // FIXME the argument is supposed to be >= 0, but ...
	fl sd_y = std::sqrt(sum_y_sq/n - sqr(sum_y/n));
	fl cov = sum_prod/n - (sum_x/n) * (sum_y/n);
	fl tmp = sd_x * sd_y;
	if(std::abs(tmp) < epsilon_fl) return 0;        // epsilon_flΪ�������ʶ��ĸ������ɱ�ʾ����Сֵ��2.22045e-16
	return cov / tmp;
}





/*
x:double;
i:unsigned int;
������spearman_aux�е�ֵ����x��i��
*/
struct spearman_aux {
	fl x;
	sz i;
	spearman_aux(fl x, sz i) : x(x), i(i) {}
};




/*�������Ϊ�����get_rankings�������ṩ����ʽ��������
�������е��βθ�Ϊָ��a��b�����Լ����ڴ�ռ䣬
����const֮�������operator������һ����a.x��b.x��a.i��b.i���޸ĵĲ����ͻᱨ��
���Է�ֹ���ǵ������
	if(a.x<b.x)
	return 1��
	else
	return 0;
*/
inline bool operator<(const spearman_aux& a, const spearman_aux& b) {     //
	return a.x < b.x;
}


/*ʵ�ֶ�����x�����е�ֵ�����������򣬲���¼�����꣬����������tmp������

input���� x:double����,ֻ����
output����tmp:double����,Ϊx.size����ʼ��ȫΪ0������;
inter����to_sort���ṹ������daxiao ��
����to_sort�ṹ��x������������������������x[i]��i=0-x.size���������ṹ���е�x��ֵ��������x;
to_sort�ṹ��i������������������������i��i=0-x.size�������Խṹ���е�i��ֵ;
Ȼ������һ��operator�������غ�������x��������
����ڰ���x����˳��ȡi��ֵ��ʹ��tmp[i]=rank :rank=0-x.size;

*/

inline flv get_rankings(const flv& x) {
	std::vector<spearman_aux> to_sort;                 //to_sortΪ�ṹ������
	VINA_FOR_IN(i, x)								   //for(uint VINA_MACROS_TMP=x.size,i=0; i<VINA_MACROS_TMP;i++)
	to_sort.push_back(spearman_aux(x[i], i));      // push_back����:��һ���µ�Ԫ�ؼӵ�vector������棬λ��Ϊ��ǰ���һ��Ԫ�ص���һ��Ԫ��
	std::sort(to_sort.begin(), to_sort.end());         //��1����һ������first����Ҫ������������ʼ��ַ��
													   //��2���ڶ�������last���ǽ����ĵ�ַ�����һ�����ݵĺ�һ�����ݵĵ�ַ��
													   //��3������������comp������ķ����������Ǵ�����Ҳ���ǽ���
													   //���������������д����Ĭ�ϵ����򷽷��Ǵ�С��������

	flv tmp(x.size(), 0);							   //����һ����СΪx.size����ʼ��ȫΪ0������
	VINA_FOR_IN(rank, to_sort)						   //for(uint VINA_MACROS_TMP=to_sort.size,rank=0; rank<VINA_MACROS_TMP;rank++)
		tmp[to_sort[rank].i] = rank;
	return tmp;
}

/*
input���� x��y:double����,ֻ����
output����pearson��double��
x��y��������󷵻ظ��Ե�����������get_rankings���ã�x1��y1(����Ϊ�����������Լ������)��
��������x1��y1������ѧ���㣨pearson���ã������ؽ����
*/

inline fl spearman(const flv& x, const flv& y) {
	return pearson(get_rankings(x), get_rankings(y));
}

#endif
