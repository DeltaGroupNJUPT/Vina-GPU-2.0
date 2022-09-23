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

#ifndef VINA_RECENT_HISTORY_H
#define VINA_RECENT_HISTORY_H

// Used by manifold, which itself is only used in "design"

#include "common.h"

struct recent_history {
	/*
	�ṹ�崫������ֵ�� initial_x_estimate�� initial_error_estimate�� lifetime��double;
	x_estimate��derror_estimate_sqr��weight��Ϊdouble���ͣ�
	               x_estimate=initial_x_estimate��
	               derror_estimate_sqr=initial_error_estimate^2;
	               weight=1/��max��1.5��lifetime����;
		
	*/
	recent_history(fl initial_x_estimate, fl initial_error_estimate, fl lifetime)
		: x_estimate(initial_x_estimate), error_estimate_sqr(sqr(initial_error_estimate)), weight(1/(std::max)(fl(1.5), lifetime)) {}
	//(std::max)(fl(1.5), lifetime)���ڶ�lifetime��1.5���Ƚ�



	/*
		input����x:double��
		inter����this_error_sqr��double��
			this_error_sqr=��x - x_estimate��^2;
			error_estimate_sqr=weight * this_error_sqr + (1 - weight) * error_estimate_sqr
	        x_estimate         = weight * x              + (1 - weight) * x_estimate;
	*/
	void add(fl x) {																	    
		fl this_error_sqr = sqr(x - x_estimate);											
		error_estimate_sqr = weight * this_error_sqr + (1 - weight) * error_estimate_sqr;
		x_estimate         = weight * x              + (1 - weight) * x_estimate;
	}

	/*	

	input����x:double��
	output����0/1��
	��������x�Ⱥ�x_estimate�����е���add�������Ǹú����������x_estimate�����Ƚϣ���x>x_estimate �򷵻�1��true����
	��������ж�(x - x_estimate)^2 < 4 * error_estimate_sqr(�����е���add�������Ǹú����������error_estimate_sqr��),
	��ȷ����1��true��������0��false��
	*/
	bool possibly_smaller_than(fl x) const {
		if(x_estimate < x) return true;
		return sqr(x - x_estimate) < 4 * error_estimate_sqr;
	}
private:
	fl x_estimate;
	fl error_estimate_sqr;
	fl weight;
};

#endif
