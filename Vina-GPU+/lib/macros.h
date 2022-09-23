/*

   Copyright (c) 2004-2010, Dr. Oleg Trott <ot14@columbia.edu>

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

*/

#ifndef VINA_MACROS_H
#define VINA_MACROS_H
/*   i��0ѭ����V�ĳ��ȴ�  
input  V ������  unit����   i���������±�  UINT����       �����������   
size_t���ͱ�ʾC���κζ������ܴﵽ����󳤶ȣ������޷���������
*/
#define VINA_FOR_IN(i, v) for(std::size_t VINA_MACROS_TMP = (v).size(), (i) = 0; (i) < VINA_MACROS_TMP; ++(i))


//  i��0ѭ����n��        input  i ��nUINT����        �����������    size_t���ͱ�ʾC���κζ������ܴﵽ����󳤶ȣ������޷���������
#define VINA_FOR(i, n)    for(std::size_t VINA_MACROS_TMP = (n),        (i) = 0; (i) < VINA_MACROS_TMP; ++(i))



//  i��0ѭ����n��            input  i��n UINT����         �����������
#define VINA_U_FOR(i, n)  for(unsigned    VINA_MACROS_TMP = (n),        (i) = 0; (i) < VINA_MACROS_TMP; ++(i))



//   i��aѭ����b��      input  i ��a ��b  UINT����       �����������   size_t���ͱ�ʾC���κζ������ܴﵽ����󳤶ȣ������޷���������
#define VINA_RANGE(i, a, b)   for(std::size_t VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))


//   i��aѭ����b��       input  i ��a ��b  UINT����       �����������
#define VINA_U_RANGE(i, a, b) for(unsigned    VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))



//   i��aѭ����b��       input  i ��a  ��  b  int����    �����������
#define VINA_I_RANGE(i, a, b) for(int         VINA_MACROS_TMP = (b), (i) = (a); (i) < VINA_MACROS_TMP; ++(i))



//  ��������������Ԫ��ֵ   input  i ��v UINT����  t:  vector<int>  ���������ܸı�    �����������
#define VINA_LOOP_CONST(t, i, v) for(t::const_iterator (i) = (v).begin(); (i) != (v).end(); ++(i))
//  ��������������Ԫ��ֵ   input  i ��v UINT����  t:  vector<int>  ���������Ըı� 
#define       VINA_LOOP(t, i, v) for(t::iterator       (i) = (v).begin(); (i) != (v).end(); ++(i))




//ִ�д�ӡһ����� x=        
#define VINA_SHOW(x)       do { std::cout << #x << " = " << (x) << std::endl; } while(false)
//ִ�д�ӡһ����� x=
#define VINA_SHOW_FAST(x)  do { std::cout << #x << " = " << (x) <<      '\n'; } while(false)
//ִ�д�ӡһ����� x=        ����Ѹ�����������Ϣ���Ǳ�׼����
#define VINA_ESHOW(x)      do { std::cerr << #x << " = " << (x) <<      '\n'; } while(false)

#endif 
