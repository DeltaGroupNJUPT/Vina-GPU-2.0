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


/*
*���������WIN32��������ͷ�ļ�Windows.h���������ͷ�ļ�unistd.h
*/

#ifdef WINDOWS
#include <Windows.h>
#else
#include <unistd.h>
#endif

/*
*����ͷ�ļ�my_pid.h
*/

#include "my_pid.h"

//���������WIN32������GetCurrentProcessId()���ص�ǰ���̵Ľ���ID����������getpid()���ص�ǰ���̵Ľ���ID
int my_pid() {
#ifdef WINDOWS
	return GetCurrentProcessId();
#else
	return getpid();
#endif
}

