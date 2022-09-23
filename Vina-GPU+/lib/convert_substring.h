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

#ifndef VINA_CONVERT_SUBSTRING_H
#define VINA_CONVERT_SUBSTRING_H

#include <cctype> // for isspace
#include <boost/lexical_cast.hpp>
#include "common.h"


struct bad_conversion {};

/*  ת�����ַ���
	input����str��string��i��j��unint
	output����tmp����������
*/
template<typename T>
T convert_substring(const std::string& str, sz i, sz j) { // indexes are 1-based, the substring should be non-null
	if(i < 1 || i > j+1 || j > str.size()) throw bad_conversion();    //�׳��쳣

	// omit leading whitespace
	while(i <= j && std::isspace(str[i-1]))      //  isspace�ж��ַ����ǲ��ǿո��ַ����ǵĻ�Ϊtrue����Ϊfalse
		++i;                                     
	/*
	 catch ���ܹ����� try ���׳����κ����͵��쳣����������쳣������������ʹ��ʡ�Ժ� ...��
	 ����������Ҫ������lexical_cast��������ת������ص�����
	 https://www.cnblogs.com/lidabo/archive/2012/12/06/2804252.html
	*/
	T tmp;
	try {
		tmp = boost::lexical_cast<T>(str.substr(i-1, j-i+1));   //lexical_cast��������ת����
	}
	catch(...) {
		throw bad_conversion();
	}
	return tmp;
}

/*�ж����ַ����Ƿ�Ϊ�գ��ǵĻ�Ϊtrue������Ϊfalse
		input����str��string��i��j��unint
		input����bool
*/
inline bool substring_is_blank(const std::string& str, sz i, sz j) { // indexes are 1-based, the substring should be non-null
	if(i < 1 || i > j+1 || j > str.size()) throw bad_conversion();  //�׳��쳣
	VINA_RANGE(k, i-1, j)											//for(k=i-1;k<j;k++)
		if(!std::isspace(str[k]))                                  // isspace�ж��ַ����ǲ��ǿո��ַ����ǵĻ�Ϊtrue����Ϊfalse
			return false;
	return true;
}


/*  ת�����ַ��������ַ���ת��unint
	input����str��string��i��j��unint
	output����tmp����������
*/
// when this was written, lexical cast to unsigned didn't work very well with "-123", etc.
template<>
inline unsigned convert_substring<unsigned>(const std::string& str, sz i, sz j) { // indexes are 1-based, the substring should be non-null
	int tmp = convert_substring<int>(str, i, j); //ģ���β�Ϊint�ͣ���strת����int��
	if(tmp < 0) throw bad_conversion();          //�׳��쳣
	return static_cast<unsigned>(tmp);           //static_cast<newType>(data)����tmp��intת��unsigend int
}

#endif
