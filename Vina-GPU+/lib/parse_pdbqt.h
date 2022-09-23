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

#ifndef VINA_PARSE_PDBQT_H
#define VINA_PARSE_PDBQT_H

#include "model.h"
/*��������parse_receptor_pdbqt
* input��path�ೣ������rigid��path�ೣ������flex
* output: model��
*/
model parse_receptor_pdbqt(const path& rigid, const path& flex); // can throw parse_error
/*������������parse_receptor_pdbqt
* input��path�ೣ������rigid
* output: model��
*/
model parse_receptor_pdbqt(const path& rigid); // can throw parse_error
/*��������parse_ligand_pdbqt
* input��path�ೣ������name
* output: model��
*/
model parse_ligand_pdbqt  (const path& name); // can throw parse_error

#endif
