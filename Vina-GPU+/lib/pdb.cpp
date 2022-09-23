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

#include "pdb.h"
#include "parse_error.h"
#include "file.h"
#include "convert_substring.h"


/*
* �ṹ��pdb��Ա����check����
* input:˫���ȸ�����min_distance����С����
* �ú������ã��ṹ��pdb�����ݳ�Ա����atoms�У������������ṹ���ʸ������С������˫���ȸ�����min_distance��ƽ���ģ����������
* output����
*/
void pdb::check(fl min_distance) const {
	VINA_FOR_IN(i, atoms) { //ѭ��atoms����������pdb_atom�ṹ��,�ֱ��Ӧ�����±�Ϊi
		const pdb_atom& a = atoms[i];//a����atoms�������i��pdb_atom�ṹ��
		VINA_RANGE(j, i+1, atoms.size()) {//ѭ��atoms���������дӵ�i+1��pdb_atom�ṹ��ʣ�µ�pdb_atom�ṹ�壨����i+1)���ֱ��Ӧ�����±�Ϊj
			const pdb_atom& b = atoms[j];//b����atoms�������j��pdb_atom�ṹ��
			fl d2 = vec_distance_sqr(a.coords, b.coords);//����ṹ��a,b�ж�Ӧcoordsʸ����ƽ�����룬���ظ�double��������d2
			
			//���d2С�������˫���ȸ�����min_distance��ƽ��
			//������ѣ�
			//The distance between �ṹ��a��id:�ṹ��a��name:�ṹ��a��element and �ṹ��b��id:�ṹ��b��name:�ṹ��b��element is d2�Ŀ��������У�
			if(d2 < sqr(min_distance)) {
				std::cout << "The distance between " 
					<< a.id << ":" << a.name << ":" << a.element
					<< " and " 
					<< b.id << ":" << b.name << ":" << b.element
					<< " is " << std::sqrt(d2) << '\n';
			}
		}
	}
}


/*
* ������string_to_pdb_atom
* input:string����str
* output���ṹ��pdb_atom
*/

pdb_atom string_to_pdb_atom(const std::string& str) {
	//���string����str���ַ���С��66���׳��쳣��bad_conversion()����ֹ��ǰ�����ִ��
	if(str.size() < 66) throw bad_conversion(); // b-factor is in 61-66
	pdb_atom tmp;//����һ���ṹ��pdb_atom�Ķ���tmp
	
	tmp.id           = convert_substring<unsigned>   (str,  7, 11);     //��string����str��( 7-1)��(11-1)�±���ַ�������ǰ���ո�ת����unsigned       �������ݸ������tmp�ĳ�Աid          
	tmp.name         = convert_substring<std::string>(str, 13, 16);     //��string����str��(13-1)��(16-1)�±���ַ�������ǰ���ո�ת����std::string    �������ݸ������tmp�ĳ�Աname        
	tmp.residue_id   = convert_substring<int>        (str, 23, 26);     //��string����str��(23-1)��(26-1)�±���ַ�������ǰ���ո�ת����int            �������ݸ������tmp�ĳ�Աresidue_id  
	tmp.residue_name = convert_substring<std::string>(str, 18, 20);     //��string����str��(18-1)��(20-1)�±���ַ�������ǰ���ո�ת����string         �������ݸ������tmp�ĳ�Աresidue_name
	tmp.coords[0]    = convert_substring<fl>         (str, 31, 38);     //��string����str��(31-1)��(38-1)�±���ַ�������ǰ���ո�ת����double         �������ݸ������tmp�ĳ�Աcoords[0]   
	tmp.coords[1]    = convert_substring<fl>         (str, 39, 46);     //��string����str��(39-1)��(46-1)�±���ַ�������ǰ���ո�ת����double         �������ݸ������tmp�ĳ�Աcoords[1]   
	tmp.coords[2]    = convert_substring<fl>         (str, 47, 54);     //��string����str��(47-1)��(54-1)�±���ַ�������ǰ���ո�ת����double         �������ݸ������tmp�ĳ�Աcoords[2]   
	tmp.b_factor     = convert_substring<fl>         (str, 61, 66);     //��string����str��(61-1)��(66-1)�±���ַ�������ǰ���ո�ת����double         �������ݸ������tmp�ĳ�Աb_factor    
	tmp.element      = convert_substring<std::string>(str, 77, 78);     //��string����str��(77-1)��(78-1)�±���ַ�������ǰ���ո�ת����string         �������ݸ������tmp�ĳ�Աelement     
	//���ؽṹ�����tmp
	return tmp;
}

/*
* ������parse_pdb
* input:path�࣬�������ʱ�����������ļ�·��
* output����ȡ�ļ����ݲ��������еĽṹ��pdb
*/
pdb parse_pdb(const path& name) {
	//����һ���ṹ��ifile�Ķ���in
	ifile in(name);

	//����һ���ṹ��pdbΪtmp
	pdb tmp;

	//����string����str
	std::string str;
	
	//�����޷�������count����ʼֵΪ0
	unsigned count = 0;

	//��������in�����ַ����浽string����str�У�ÿ�ζ���һ�����У�ֱ���������ļ�������־
	while(std::getline(in, str)) {
	//count���ڱ����к�
		++count;
	//�жϵ�ǰ�ַ����Ƿ�������һ�����������ַ�����ATOM  ������"HETATM"��ͷ������ִ��if�����Ĳ��������������ȡ��һ��
		if(starts_with(str, "ATOM  ") || starts_with(str, "HETATM")) {
			try {
				tmp.atoms.push_back(string_to_pdb_atom(str));
				//��pdb_atom�ṹ�����string_to_pdb_atom����tmp�����ݳ�Ա������atoms
			}
			catch(...) { // bad_conversion, presumably; but can lexical_cast throw its own errors?
				throw parse_error(name, count, "ATOM syntax incorrect");
				//�ڽ�strת����string_to_pdb_atom���ݳ�Ա�Ĺ��̷����쳣���׳��쳣��parse_error()����ֹ��ǰ�����ִ��
			}
		}
	}
	return tmp;
	//���ض����ļ����ݺ��pdb����tmp
}
