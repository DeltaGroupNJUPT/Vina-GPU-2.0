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

#ifndef VINA_PRECALCULATE_H
#define VINA_PRECALCULATE_H

#include "scoring_function.h"
#include "matrix.h"
/*
  �ṹ��precalculate_element���������Ϊunint n��double factor
		fast����double������n������Ԫ�أ���ʼ��ȫΪ0��
		smooth�����ṹ����������������Ա����first��second��������������double���ͣ�n������Ԫ�أ���ʼ��ȫΪ0��
		factor����double���ͣ����������factor��
*/
//���pair
// pair<int, double> p1;  ʹ��Ĭ�Ϲ��캯������һ��������Ҫ����2�����ݵ�ʱ��ѡ��pairʵ����һ���ṹ�壬
//��Ҫ��������Ա������Ϊfirst��second ��Ϊ��ʹ��struct����class�����Կ���ֱ��ʹ��pair�ĳ�Ա����
//p1.first = 1;
//p1.second = 2.5;
//cout << p1.first << p1.second << endl;
//���1  2.5

struct precalculate_element {
	precalculate_element(sz n, fl factor_) : fast(n, 0), smooth(n, pr(0, 0)), factor(factor_) {}

	/*
	input����r2��double
	inter����i��unint
	output����fast[i]��double�����е�iλ����ֵ��
	fast����double����������Ԫ�س�ʼ��ȫΪ0��
	fast.size()����fast����Ԫ�ظ�����
	�������ж�r2 * factor�Ƿ�С��fast������Ԫ�ظ������ǵĻ��ͽ���ִ������ĳ��򣬷�����ֹ����
	Ȼ���factor * r2ǿ��ȡ���õ�i��
	����ж�i �Ƿ�С��fast.size()���ǵĻ��ͽ���ִ������ĳ��򼴷���fast[i]��ֵ��������ֹ����
	*/
	fl eval_fast(fl r2) const {
		assert(r2 * factor < fast.size());//�����������ش�������ֹ����ִ��
		sz i = sz(factor * r2);           //r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
		assert(i < fast.size());          //�����������ش�������ֹ����ִ��
		return fast[i];
	}


	/*
	input����r2��double��
	inter����r2_factored��rem��double�� i1��i2��unint, p1��p2�����캯����
	output����pr(e, dor)������������ֵ����Ϊdouble���ͣ�
	smooth.size()����fast����Ԫ�ظ�����
	smooth�����ṹ����������������Ա����first��second��������������double���ͣ�n������Ԫ�أ���ʼ��ȫΪ0��
	epsilon_fl�����������ʶ��ĸ������ɱ�ʾ����Сֵ���ҵļ������ʶ�����С��������2.22045e-16

	�������ж�r2_factored=r2 * factor�Ƿ�С��smooth������Ԫ�ظ������ǵĻ��ͽ���ִ������ĳ��򣬷�����ֹ����
	Ȼ���r2_factoredȡ���õ�i1����i2=i1+1��
	Ȼ���ж�i1��i2 �Ƿ�С��smooth.size()���ǵĻ��ͽ���ִ������ĳ��򼴼���remΪfactor * r2��С�����ݣ�������ֹ����
	����ж�rem�Ƿ����-epsilon_fl��С��1 + epsilon_fl���ǵĻ����ṹ������smooth[i1]��smooth[i2]��������Ϊp1��p2��������p1.first��p1.second��p2.first��p2.second��������ֹ����
				e   = smooth[sz(factor * r2)].first  + (factor * r2 - sz(factor * r2)) * (smooth[sz(factor * r2)+1].first  - smooth[sz(factor * r2)].first );
				dor = smooth[sz(factor * r2)].second  + (factor * r2 - sz(factor * r2)) * (smooth[sz(factor * r2)+1].second   - smooth[sz(factor * r2)].second  );

	*/

	pr eval_deriv(fl r2) const {
		fl r2_factored = factor * r2;         
		assert(r2_factored + 1 < smooth.size());              //�����������ش�������ֹ����ִ��
		sz i1 = sz(r2_factored);
		sz i2 = i1 + 1;                                      // r2 is expected < cutoff_sqr, and cutoff_sqr * factor + 1 < n, so no overflow
		assert(i1 < smooth.size());
		assert(i2 < smooth.size());
		fl rem = r2_factored - i1;
		assert(rem >= -epsilon_fl);                          //epsilon_fl�������ʶ��ĸ������ɱ�ʾ����Сֵ���ҵļ������ʶ�����С��������2.22045e-16
		assert(rem < 1 + epsilon_fl);
		const pr& p1 = smooth[i1];	                         //	�����&���Ǹ�ֵ��ַ������һ�����ã�Ҳ����˵�Ǳ���			
		const pr& p2 = smooth[i2];                           //
		fl e   = p1.first  + rem * (p2.first  - p1.first);      
		fl dor = p1.second + rem * (p2.second - p1.second); 
		return pr(e, dor);
	}




/*
input����rs��double������
inter����n��unint��dor��delta��r��f1��f2��double
smooth.size()����smooth����Ԫ�ظ�����
smooth�����ṹ����������������Ա����first��second��������������double���ͣ�����Ԫ�س�ʼ��ȫΪ0��
fast����double����������Ԫ�س�ʼ��ȫΪ0��
fast.size()����fast����Ԫ�ظ�����
���ȼ���smooth����Ԫ�ظ�������ֵ��n��
Ȼ���ж�����rs����Ԫ�ظ�����fast����Ԫ�ظ����Ƿ�Ϊn����Ļ�����ֹ�����ǵĻ������������
���ͨ��һ��0<=i<nѭ������dor��fast����Ԫ�أ�����dorΪsmooth[i].second���޸�dor�͵��޸�smooth[i].second��


		����dor�����Ƚ������ĵ�һ��Ԫ�غ����һ��Ԫ�ض���ֵΪ0������

						dor=��smooth[i+1].first - smooth[i-1].first��/��rs[i+1] - rs[i-1]��*rs[i]����

	    ����fast���������ȵ�i == n-1ʱf2=0������f2=smooth[i+1].first�����ԣ�

				     	fast[i]=(smooth[i].first+f2)/2��
*/

	void init_from_smooth_fst(const flv& rs) {
		sz n = smooth.size();
		VINA_CHECK(rs.size() == n);														//assert����
		VINA_CHECK(fast.size() == n);												    //assert����
		VINA_FOR(i, n) {																//for(uint VINA_MACROS_TMP=n,i=0; i<VINA_MACROS_TMP;i++)ѭ��
			// calculate dor's
			fl& dor = smooth[i].second;                                                 //�����&���Ǹ�ֵ��ַ������һ�����ã�Ҳ����˵�Ǳ������޸�dor�͵��޸�smooth[i].second	
			if(i == 0 || i == n-1)                                                      //�����ĵ�һ��Ԫ�غ����һ��Ԫ�ض�Ϊ0
				dor = 0;
			else {
				fl delta = rs[i+1] - rs[i-1];                                           //���������һ��Ԫ�ص�������ֵ�Ĳ�ֵ��������������Ԫ�ظ�ֵ�� r
				fl r = rs[i];
				dor = (smooth[i+1].first - smooth[i-1].first) / (delta * r);		   
			}
			// calculate fast's from smooth.first's
			fl f1 = smooth[i].first;
			fl f2 = (i+1 >= n) ? 0 : smooth[i+1].first;
			fast[i] = (f2 + f1) / 2;
		}
	}


/*			����smooth.first����Сֵ������
	input�����ޣ�
	inter����i��unint��
	output����tmp��unint��
	smooth.size()����smooth����Ԫ�ظ�����
	smooth�����ṹ����������������Ա����first��second��������������double���ͣ�����Ԫ�س�ʼ��ȫΪ0��
	
	�ȶ�tmp����ֵ0��
	���ͨ��һ��0<=i_inv<smooth.sizeѭ����tmp���и�ֵ������������£�

								    i = smooth.size() - i_inv - 1; 

					��i_inv == 0 ���� smooth[i].first < smooth[tmp].firstʱ��

									tmp = i;
*/

	sz min_smooth_fst() const {
		sz tmp = 0;                                                 // returned if smooth.empty()
		VINA_FOR_IN(i_inv, smooth) {                               //for(uint VINA_MACROS_TMP=smooth.size,i_inv=0; i_inv<VINA_MACROS_TMP;i_inv++)
			sz i = smooth.size() - i_inv - 1;                      // i_inv < smooth.size()  => i_inv + 1 <= smooth.size()
			if(i_inv == 0 || smooth[i].first < smooth[tmp].first)
				tmp = i;
		}
		return tmp;
	}


	/*
	input����rs��double������left��right��double��
	inter����tmp����smoothͬһԪ��������double������ min_index��unint�� optimal_r��r��double��

	smooth.size()����smooth����Ԫ�ظ�����
	smooth�����ṹ����������������Ա����first��second��������������double���ͣ�����Ԫ�س�ʼ��ȫΪ0��

	�����ǵ�������min_smooth_fst()����smooth.first����Сֵ��������ֵ��min_index��
	���ж������double����rs��Ԫ�ظ����Ƿ�����������ֵ���Լ��Ƿ�smooth��Ԫ�ظ�����ȣ��ǵĻ���ѵ������򣬷�����ֹ����
	Ȼ��ͨ��0<=i<smooth.sizeѭ����tmp[i]��ֵ��
					    //����rs[i]�޸Ĳ����޸������rs[i]����Ϊ��ͨ��r������
						if (rs[i] < rs[min_index] - left ) 
								rs[i] = r+left;
						else if(rs[i] > rs[min_index]+ right)
								rs[i] = rs[i]-right;
						else  
						        rs[i] = rs[min_index];

						if(rs[i] < 0) 
							rs[i] = 0;
						if(rs[i] > rs.back())
						    rs[i] = rs.back();                   //rs.back() ����rs��ĩһ��Ԫ��

						tmp[i] = eval_deriv(rs[i]^2).first;      //eval_deriv()��������pr(e, dor)������.firstӦ��ȡֵe��ֵ

	���ͨ��0<=i<smooth.sizeѭ����smooth[i].first]��ֵ��

						smooth[i].first = tmp[i];

	*/

	void widen_smooth_fst(const flv& rs, fl left, fl right) {
		flv tmp(smooth.size(), 0);								// the new smooth[].first
		sz min_index = min_smooth_fst();						// smooth.first����Сֵ������
		VINA_CHECK(min_index < rs.size());					    //assert // won't hold for n == 0
		VINA_CHECK(rs.size() == smooth.size());					//assert
		fl optimal_r   = rs[min_index];							//
		VINA_FOR_IN(i, smooth) {                                //for(uint VINA_MACROS_TMP=smooth.size,i=0; i<VINA_MACROS_TMP;i++)
			fl r = rs[i];
			if     (r < optimal_r - left ) r += left;
			else if(r > optimal_r + right) r -= right;
			else                           r = optimal_r;

			if(r < 0) r = 0;
			if(r > rs.back()) r = rs.back();                   //back() ������ĩһ��Ԫ��

			tmp[i] = eval_deriv(sqr(r)).first;                 //eval_deriv()��������pr(e, dor)
		}
		VINA_FOR_IN(i, smooth)									//for(uint VINA_MACROS_TMP=smooth.size,i=0; i<VINA_MACROS_TMP;i++)                         
			smooth[i].first = tmp[i];
	}



    /*���������Ҫ����������widen_smooth_fst������init_from_smooth_fst������
	input����rs��double������left��right��double��
	*/

	void widen(const flv& rs, fl left, fl right) {
		widen_smooth_fst(rs, left, right);
		init_from_smooth_fst(rs);
	}


	flv fast;                           //double����
	prv smooth;                         //�ṹ������
	fl factor;                          //double
};



/*
      input����sf��scoring_function�ṹ��;   v:�������ʶ���double���ɱ�ʾ�����ֵ���ҵļ������ʶ������double����00DB17EE,
	  factor_:double���ͳ�ʼֵ32�� 
	  data����ģ�����Ϊprecalculate_element�ṹ���triangular_matrix�ṹ�壻
	  m_cutoff_sqr=weighted_terms�ṹ�����ع���cutoff���������ķ���ֵcutoff_��ƽ����
	  n=unint��m_cutoff_sqr*factor_��+3��factor=factor_(��ʼֵ)��
	  data�ṹ��������m_data={ԭ��������unint*��ԭ��������unint-1��/2�����루n,factor_��Ϊprecalculate_element�ṹ��}��
	  m_dim=ԭ��������unint��
	  m_atom_typing_used��ԭ��������ö�����ͣ�
		�ʼҪ��֤factor���ڼ������ʶ��ĸ������ɱ�ʾ����Сֵ��2.22045e-16��������nҪ����sz(m_cutoff_sqr*factor) + 1��m_cutoff_sqr*factor + 1 
	  

*/
struct precalculate {
	    precalculate(const scoring_function& sf, fl v = max_fl, fl factor_ = 32) : // sf should not be discontinuous, even near cutoff, for the sake of the derivatives
		m_cutoff_sqr(sqr(sf.cutoff())),     //����weighted_terms�ṹ�����ع���cutoff���������ķ���ֵcutoff_��ƽ��
		n(sz(factor_ * m_cutoff_sqr) + 3),  // sz(factor * r^2) + 1 <= sz(factor * cutoff_sqr) + 2 <= n-1 < n  // see assert below
		factor(factor_),

		data(num_atom_types(sf.atom_typing_used()), precalculate_element(n, factor_)),//num_atom_types������������weighted_terms�ṹ����atom_typing_used_����ʾԭ��������unint
			//strictly_triangular_matrix(sz n, const T& filler_val) : m_data(n*(n - 1) / 2, filler_val), m_dim(n) {}
			
			m_atom_typing_used(sf.atom_typing_used()) {                                   //

		VINA_CHECK(factor > epsilon_fl);              //factor���ڼ������ʶ��ĸ������ɱ�ʾ����Сֵ���ҵļ������ʶ�����С��������2.22045e-16
		VINA_CHECK(sz(m_cutoff_sqr*factor) + 1 < n); // cutoff_sqr * factor is the largest float we may end up converting into sz, then 1 can be added to the result
		VINA_CHECK(m_cutoff_sqr*factor + 1 < n);

		flv rs = calculate_rs();                                                   //rs[i]=��i/factor��^(1/2);(0-n)
														                          
		VINA_FOR(t1, data.dim())                                                  //for(t1=0;t1<m_dim;t1++)  m_dim=ԭ��������unint
			VINA_RANGE(t2, t1, data.dim()) {                                      //for(t2=t1;i<m_dim;t2++) Ƕ��ѭ�� m_dim����
				precalculate_element& p = data(t1, t2);	                         //���� data��t1, t2���ṹ��                        
				// init smooth[].first					                          
				VINA_FOR_IN(i, p.smooth)                                          //for(i=0;i<smooth.size;i++)
					p.smooth[i].first = (std::min)(v, sf.eval(t1, t2, rs[i]));   // sf.eval(t1, t2, rs[i]���س��ۼӽ����
																				  //��ȡ�������ȡ�����doubleֵ�ͳ��ۼӽ������Сֵ

				// init the rest
				p.init_from_smooth_fst(rs);                                       //��ʼ��
			}
	}

/*	
	input����type_pair_index��unint��r2��double��
    m_cutoff_sqr=weighted_terms�ṹ�����ع���cutoff���������ķ���ֵcutoff_��ƽ����
	��֤�����r2<= m_cutoff_sqr��
	��󷵻�fast[sz(factor * r2)},	precalculate_element�ṹ����fast����double����������Ԫ�س�ʼ��ȫΪ0��
*/
	fl eval_fast(sz type_pair_index, fl r2) const {
		assert(r2 <= m_cutoff_sqr);
		return data(type_pair_index).eval_fast(r2);
	}

	/*
		input����type_pair_index��unint��r2��double��
		m_cutoff_sqr=weighted_terms�ṹ�����ع���cutoff���������ķ���ֵcutoff_��ƽ����
		��֤�����r2<= m_cutoff_sqr��
		��󷵻�(e, dor)��
		e   = smooth[sz(factor * r2)].first  + (factor * r2 - sz(factor * r2)) * (smooth[sz(factor * r2)+1].first  - smooth[sz(factor * r2)].first );
		dor = smooth[sz(factor * r2)].second  + (factor * r2 - sz(factor * r2)) * (smooth[sz(factor * r2)+1].second   - smooth[sz(factor * r2)].second  );

	*/

	pr eval_deriv(sz type_pair_index, fl r2) const {
		assert(r2 <= m_cutoff_sqr);
		return data(type_pair_index).eval_deriv(r2);
	}


	sz index_permissive(sz t1, sz t2) const { return data.index_permissive(t1, t2); }//���� (t1 < t2) ? index(t1, t2) : index(t2, t1)
	atom_type::t atom_typing_used() const { return m_atom_typing_used; }             //���� ԭ��������ö������
	fl cutoff_sqr() const { return m_cutoff_sqr; }                                   //����m_cutoff_sqr=weighted_terms�ṹ�����ع���cutoff���������ķ���ֵ
																				    //cutoff_��ƽ����

	/*
	
	
	*/
	void widen(fl left, fl right) {
		flv rs = calculate_rs();                           //rs[i]=��i/factor��^(1/2);(0-n)   
		VINA_FOR(t1, data.dim())                          //for(t1=0;t1<m_dim;t1++)��m_dim=ԭ��������unint
			VINA_RANGE(t2, t1, data.dim())                //for(t2=t1;i<m_dim;t2++) Ƕ��ѭ�� m_dim����
				data(t1, t2).widen(rs, left, right);     //����������widen_smooth_fst������init_from_smooth_fst����
	}


public:

	/* 
	��������double����tmp��tmp�е�ֵΪ��i/factor��^(1/2)
	*/
	flv calculate_rs() const {
		flv tmp(n, 0);
		VINA_FOR(i, n)											//for(uint VINA_MACROS_TMP=n,i=0; i<VINA_MACROS_TMP;i++)
			tmp[i] = std::sqrt(i / factor);                     //tmp[i]=(i / factor)^(1/2)
		return tmp;
	}



	fl m_cutoff_sqr;                                          //double
	sz n;													  //unint 
	fl factor;												  //double
	atom_type::t m_atom_typing_used;                          //atom_type::tö����

	triangular_matrix<precalculate_element> data;             //ģ�����Ϊprecalculate_element�ṹ���triangular_matrix�ṹ��

};

#endif
