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

#ifndef VINA_PARSE_ERROR_H
#define VINA_PARSE_ERROR_H

#include "common.h"
/*
* �ṹ��parse_error
*���ݳ�Ա����:boost::filesystem��path��file;�޷�����������line;string��reason
*���캯��parse_error�ڴ���parse_error����ʱ�����б��ṩpath,line,reasonʱִ��Ĭ�ϳ�ʼ������֮���γ�ʼ����Աpath,line,reason
*path��Ĺ��캯���ɽ���char*���ͺ�string���͵Ĳ������죬Ҳ������һ���ַ���������Χ,·���ķָ����constexpr preferred_separator����,UNIX����б��(/),WINDOWS�Ƿ�б��(\),C++����Ҫת��;
*��ʼ��path�����ʱ�������б����ļ��ڵ��Եĵ�ַ 
*/

struct parse_error {
	path file;
	unsigned line;
	std::string reason;
	parse_error(const path& file_, unsigned line_, const std::string& reason_ = "") : file(file_), line(line_), reason(reason_) {}
private:
	parse_error() {}
};

#endif
