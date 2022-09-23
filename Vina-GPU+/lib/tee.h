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

#ifndef VINA_TEE_H
#define VINA_TEE_H

#include <iostream>
#include "file.h"


/*


*/
struct tee {
	ofile* of;
	tee() : of(NULL) {}                        //�������ļ�

/*
     ����ofile�ṹ��of�����ڴ���һ��������name��ֵ���ļ�
*/
	void init(const path& name) {
		of = new ofile(name);                 
	}

	/*
	      ������һ���麯������ɾ��of�ṹ��
	*/
	virtual ~tee() { delete of; } 


	/*
	  flush() �ǰѻ�����������ǿ�����, ��Ҫ����IO�У�����ջ��������ݣ�
	  һ���ڶ�д��(stream)��ʱ���������ȱ��������ڴ��У��ٰ�����д���ļ��У�
	  �������ݶ����ʱ�򲻴�����������Ѿ�д���ˣ���Ϊ����һ�����п��ܻ������ڴ�����������С�
	  ��ʱ������������close()�����ر��˶�д������ô�ⲿ�����ݾͻᶪʧ������Ӧ���ڹرն�д��֮ǰ��flush()��
	  ��������һ���ļ������ڻ�����������д���������ļ���
	*/

	void flush() {
		std::cout << std::flush;                            
		if(of)                      
			(*of) << std::flush;
	}



	/*
	endl() ���в�ˢ�»���������flush�������
	*/

	void endl() {                 
		std::cout << std::endl;
		if(of)
			(*of) << std::endl;
	}

	/*
	ͨ���ú���ʵ��of�ṹ�����ļ���a�����Ŀ��Ʒ��Ÿ�ʽ��������Ž����������ӣ�
	https://blog.csdn.net/yedawei_1/article/details/105538439
	*/

	void setf(std::ios::fmtflags a) {   //��ʽ������Ϣ��ö������fmtflags ��Ӱ�쵽��ν����������еĸ�ʽ��
		                                //�������������еĸ�ʽ������������16���ƻ���10���Ʊ�ʾ���������ǿ�ѧ���������Ƕ�����ʽ�ȣ�
										//a�����������������ʽ
		std::cout.setf(a);
		if(of)
			of->setf(a);
	}
	/*
	ͨ���ú���ʵ��of�ṹ�����ļ���a��b��Ͽ��Ʒ��Ÿ�ʽ��������Ž����������ӣ�
	https://blog.csdn.net/yedawei_1/article/details/105538439
	*/
	void setf(std::ios::fmtflags a, std::ios::fmtflags b) {
		std::cout.setf(a, b);
		if(of)
			of->setf(a, b);
	}
};



/*�ع�tee�ṹ��
	input ����out��tee�ṹ�壬x���������ͣ�
	���x�������������ļ��ͽ�xд���ļ�����󷵻�out�ṹ��



*/
template<typename T>
tee& operator<<(tee& out, const T& x) {
	std::cout << x;
	if(out.of)
		(*out.of) << x;
	return out;
}

#endif
