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

#ifndef VINA_TREE_H
#define VINA_TREE_H

#include "conf.h"
#include "atom.h"



/*
	orientation_m����mat�ṹ�壻orientation_q������Ԫ��
	origin=origin_����vec�ṹ����������������Ԫ�������鲿��
	��ʼ����
	orientation_q������Ԫ����(1, 0, 0, 0)��
	orientation_m����mat�ṹ�壨�ԽǾ���
	data[0]=1 ��data[3]=0 ��data[6]=0 ;
	data[1]=0 ��data[4]=1 ��data[7]=0 ;
	data[2]=0 ��data[5]=0 ��data[8]=1 ;
*/

struct frame {
	frame(const vec& origin_) : origin(origin_), orientation_q(qt_identity), orientation_m(quaternion_to_r3(qt_identity)) {}
	/*
		input����local_coords��vec�ṹ��;
		output����tmp��vec�ṹ��;
		demo��
		origin��data0[0],data0[1],data0[2]��
		orientation_m��data1[0],data1[3],data1[6]
				       data1[1],data1[4],data1[7]
					   data1[2],data1[5],data1[8]��
		local_coords��data2[0],data2[1],data2[2]��

		tmp(data0[0]+data1[0]*data2[0] + data1[3]*data2[1] + data1[6]*data2[2],
		    data0[1]+data1[1]*data2[0] + data1[4]*data2[1] + data1[7]*data2[2],
		    data0[2]+data1[2]*data2[0] + data1[5]*data2[1] + data1[8]*data2[2])
	*/
	vec local_to_lab(const vec& local_coords) const {
		vec tmp;
		tmp = origin + orientation_m*local_coords; 
		return tmp;
	}

	/*
		input����local_direction��vec�ṹ��;
		output����tmp��vec�ṹ��;
		demo:
		orientation_m��data1[0],data1[3],data1[6]
					   data1[1],data1[4],data1[7]
					   data1[2],data1[5],data1[8]��
		local_direction��data2[0],data2[1],data2[2]��

		tmp(data1[0]*data2[0] + data1[3]*data2[1] + data1[6]*data2[2],
			data1[1]*data2[0] + data1[4]*data2[1] + data1[7]*data2[2],
			data1[2]*data2[0] + data1[5]*data2[1] + data1[8]*data2[2])
	
	*/
	vec local_to_lab_direction(const vec& local_direction) const {
		vec tmp;
		tmp = orientation_m * local_direction;
		return tmp;
	}


	const qt& orientation() const { return orientation_q; }//����orientation_q��Ԫ��
	const vec& get_origin() const { return origin; }       //����origin�ṹ��
	const mat& get_orientation_m() const { return orientation_m; }
protected:
	vec origin;
	/*	������Ԫ��orientation_q=q��orientation_m�ṹ��
		input����q����Ԫ��
	*/
	void set_orientation(const qt& q) { // does not normalize the orientation
		orientation_q = q;
		orientation_m = quaternion_to_r3(orientation_q);
	}
	mat orientation_m;
	qt  orientation_q;
};

/*
	begin��end��unint
*/
struct atom_range {
    sz begin;
    sz end;
	atom_range(sz begin_, sz end_) : begin(begin_), end(end_) {}
	template<typename F>
	/*	ת��
		input����f���κ�����
		diff��unint
	*/
	void transform(const F& f) {
		sz diff = end - begin;
		begin = f(begin);
		end   = begin + diff;
	}
};


/*	ԭ�ӽṹ���̳�atom_range��frame�ṹ��

    frame�ṹ���У�
	origin=origin_����vec�ṹ����������������Ԫ�������鲿��
	��ʼ����
	orientation_q������Ԫ����(1, 0, 0, 0)��
	orientation_m����mat�ṹ�壨�ԽǾ���
	data[0]=1 ��data[3]=0 ��data[6]=0 ;
	data[1]=0 ��data[4]=1 ��data[7]=0 ;
	data[2]=0 ��data[5]=0 ��data[8]=1 ;
	
	atom_range�ṹ���У�
	begin=begin_��end=end_��

*/
struct atom_frame : public frame, public atom_range {
	atom_frame(const vec& origin_, sz begin_, sz end_) : frame(origin_), atom_range(begin_, end_) {}
	/*
		input����atoms��atom�ṹ��������coords��vec�ṹ������
		demo��
		origin��data0[0],data0[1],data0[2]��
		orientation_m��data1[0],data1[3],data1[6]
					   data1[1],data1[4],data1[7]
					   data1[2],data1[5],data1[8]��
		coords ��data2[0],data2[1],data2[2]��

		tmp(data0[0]+data1[0]*data2[0] + data1[3]*data2[1] + data1[6]*data2[2],
		    data0[1]+data1[1]*data2[0] + data1[4]*data2[1] + data1[7]*data2[2],
		    data0[2]+data1[2]*data2[0] + data1[5]*data2[1] + data1[8]*data2[2])
			����coords�ṹ��һ��Ԫ�ص����ݣ��У�begin��end����
	*/
	void set_coords(const atomv& atoms, vecv& coords) const {
		VINA_RANGE(i, begin, end)
			coords[i] = local_to_lab(atoms[i].coords);
	}

	/*
		input����coords��forces��vec�ṹ��������
		output����tmp��<vec,vec>�ṹ��������
		�ȳ�ʼ��tmp�ṹ������<(0,0,0),(0,0,0)>��tmp��������Ա����first��second
		demo��
			forces{��0��0��0������1��1��1������2��2��2��}
			coords{��1��1��1������1��1��1������2��2��2��}
			origin��1,1,1��
			coords[i] - origin={(1-1,1-1,1-1),(1-1,1-1,1-1),(2-1,2-1,2-1)}={��0��0��0������0��0��0������1��1��1��}
			tmp<(0+0+1+2=3,0+0+1+2=3,0+0+1+2=3),(0+0+0+0,0+0+0+0,0+0+0+0)>=<(3,3,3),(0,0,0)>
	
	*/
	vecp sum_force_and_torque(const vecv& coords, const vecv& forces) const {
		vecp tmp;
		tmp.first.assign(0);
		tmp.second.assign(0);
		VINA_RANGE(i, begin, end) {
			tmp.first  += forces[i]; 
			tmp.second += cross_product(coords[i] - origin, forces[i]);
		}
		return tmp;
	}
};

/*  ���壬�̳�atom_frame�ṹ��
	origin_��vec�ṹ�壬begin_��end_��unint��
	frame�ṹ���е�origin=origin_��atom_range�ṹ����begin=begin_��end=end_
*/
struct rigid_body : public atom_frame {
	rigid_body(const vec& origin_, sz begin_, sz end_) : atom_frame(origin_, begin_, end_) {}

	/*��������
	input����atoms��atom�ṹ��������coords��vec�ṹ��������c��rigid_conf�ṹ��
		��atom_frame.frame.origin=rigid_conf.position
		��������Ԫ��orientation_q=rigid_conf.orientation��orientation_m�ṹ��
		������coords�ṹ������
		     demo��
		     origin��data0[0],data0[1],data0[2]��
		     orientation_m��data1[0],data1[3],data1[6]
		     			   data1[1],data1[4],data1[7]
		     			   data1[2],data1[5],data1[8]��
		     coords ��data2[0],data2[1],data2[2]��
		     
		     tmp(data0[0]+data1[0]*data2[0] + data1[3]*data2[1] + data1[6]*data2[2],
		    data0[1]+data1[1]*data2[0] + data1[4]*data2[1] + data1[7]*data2[2],
		    data0[2]+data1[2]*data2[0] + data1[5]*data2[1] + data1[8]*data2[2])
			����coords�ṹ��һ��Ԫ�ص����ݣ��У�begin��end����
	*/
	void set_conf(const atomv& atoms, vecv& coords, const rigid_conf& c) {
		origin = c.position;
		set_orientation(c.orientation);
		set_coords(atoms, coords);
	}

	void count_torsions(sz& s) const {} // do nothing�������û��ʲô����

	/*����������
		input����force_torque��<vec, vec>�ṹ��������������Ա����first��second��c:rigid_change�ṹ��
	*/
	void set_derivative(const vecp& force_torque, rigid_change& c) const {
		c.position     = force_torque.first;
		c.orientation  = force_torque.second;
	}
};
/*����ṹ���̳�atom_frame�ṹ��
	origin_��axis_root��diff��vec�ṹ�壬begin_��end_��unint��nrm��double��
	��frame�ṹ���е�origin=origin_��atom_range�ṹ����begin=begin_��end=end_
	��diff=origin�ṹ���ȥaxis_root�ṹ�壨��-���ع���
	��nrm=diff�ṹ����data�����Ԫ��ƽ���ͣ�������(�ռ�����ģ)
	��axis��vec �ṹ��
	demo��
		diff��2��3��4����1/nrm=2;
		 axis=��4,6,8)
*/
struct axis_frame : public atom_frame {
	axis_frame(const vec& origin_, sz begin_, sz end_, const vec& axis_root) : atom_frame(origin_, begin_, end_) {
		vec diff; diff = origin - axis_root;
		fl nrm = diff.norm();
		VINA_CHECK(nrm >= epsilon_fl);     //Ҫ����2.22045e-16
		axis = (1/nrm) * diff;
	}

	/* ����������
		input����force_torque��<vec,vec>�ṹ��������c��double
		axis��vec �ṹ��
		demo��
		force_torque.second��1��2��3����axis��2��3��4��
		c=1*2+2*3+3*4=20
	*/
	void set_derivative(const vecp& force_torque, fl& c) const {
		c = force_torque.second * axis;
	}
	vec get_axis() { return axis; } // ����axis

protected:
	vec axis;
};
/*   �̳�axis_frame�ṹ��
	origin_��axis_root��relative_axis��relative_origin��vec�ṹ�壬begin_��end_��unint��parent��frame�ṹ�壻
	��frame�ṹ���е�origin=origin_��atom_range�ṹ����begin=begin_��end=end_
	��diff=origin�ṹ���ȥaxis_root�ṹ�壨��-���ع���
	��axis_frame�ṹ����nrm=diff�ṹ����data�����Ԫ��ƽ���ͣ�������(�ռ�����ģ)
	��axis_frame�ṹ����axis��vec �ṹ��
		demo��
			diff��2��3��4����1/nrm=2;
			 axis=��4,6,8)
	���ж�frame�ṹ����orientation_q�Ƿ����(1, 0, 0, 0)
	��relative_axis=axis
	��relative_origin��
	  demo��origin��4��5��6����frame�ṹ����origin��1��1��1��
	  relative_origin=��3��4��5��
*/
struct segment : public axis_frame {
	segment(const vec& origin_, sz begin_, sz end_, const vec& axis_root, const frame& parent) : axis_frame(origin_, begin_, end_, axis_root) {
		VINA_CHECK(eq(parent.orientation(), qt_identity)); // the only initial parent orientation this c'tor supports
		relative_axis = axis;
		relative_origin = origin - parent.get_origin();
	}

	/*������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
		input����parent��frame�ṹ�壬atoms��coords��atom�ṹ��������c��������
		torsion��double
		tmp���Ƕ���Ԫ��*frame�ṹ����orientation_q��Ԫ��
	
	*/
	void set_conf(const frame& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
		const fl torsion = *c;  //ָ���Ԫ�ص�����
		++c;                    //ָ����һ��Ԫ��
		origin = parent.local_to_lab(relative_origin);          //�������
		axis = parent.local_to_lab_direction(relative_axis);    //������
		qt tmp = angle_to_quaternion(axis, torsion) * parent.orientation();//�Ƕ���Ԫ��*frame�ṹ����orientation_q��Ԫ��
		quaternion_normalize_approx(tmp); // normalization added in 1.1.2  �ж���Ԫ���Ƿ���Ʊ�׼��
		//quaternion_normalize(tmp); // normalization added in 1.1.2
		set_orientation(tmp);
		set_coords(atoms, coords);
	}
	/*Ť�ؼ���
	 input����s:unint
	*/
	void count_torsions(sz& s) const {
		++s;
	}
public:
	vec relative_axis;
	vec relative_origin;
};


/*�̳�axis_frame�ṹ��

*/
struct first_segment : public axis_frame {
	first_segment(const segment& s) : axis_frame(s) {}
	first_segment(const vec& origin_, sz begin_, sz end_, const vec& axis_root) : axis_frame(origin_, begin_, end_, axis_root) {}
	/*���ýǶ���Ԫ��orientation_q=q��orientation_m�ṹ��
		input����atoms��atom�ṹ�壬coords��vec�ṹ�壬torsion��double
	*/
	void set_conf(const atomv& atoms, vecv& coords, fl torsion) {
		set_orientation(angle_to_quaternion(axis, torsion));//�Ƕ���Ԫ��
  	}

	/*Ť�ؼ���
	input����s:unint
	*/
	void count_torsions(sz& s) const {
		++s;
	}
};

/*
 ������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
*/
template<typename T> // T == branch
void branches_set_conf(std::vector<T>& b, const frame& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
	VINA_FOR_IN(i, b)                            //for(i=0;i<v.size;i++)
		b[i].set_conf(parent, atoms, coords, c);
}


/*��Ҫ����b�����͵��ú���
	input����b������������origin��vec�ṹ�壬forces��out��coords��vec�ṹ������
	d��������
*/
template<typename T> // T == branch
void branches_derivative(const std::vector<T>& b, const vec& origin, const vecv& coords, const vecv& forces, vecp& out, flv::iterator& d) { // adds to out
	VINA_FOR_IN(i, b) {											//for(i=0;i<v.size;i++)
		vecp force_torque = b[i].derivative(coords, forces, d);
		out.first  += force_torque.first;
		vec r; r = b[i].node.get_origin() - origin;
		out.second += cross_product(r, force_torque.first) + force_torque.second;
	}
}


/*
    node���������ͱ��� ��children��tree�ṹ������
	node=node_
*/
template<typename T> // T == segment
struct tree {
	T node;
	std::vector< tree<T> > children;
	tree(const T& node_) : node(node_) {}
	/*	���ã�����node����ȥ��������ڵĺ���
		input����parent��frame�ṹ�壬atoms��atom�ṹ��������coords��vec�ṹ������
		c��������
	*/
	void set_conf(const frame& parent, const atomv& atoms, vecv& coords, flv::const_iterator& c) {
		node.set_conf(parent, atoms, coords, c);
		branches_set_conf(children, node, atoms, coords, c);
	}

	/*	����node����ȥ��������ڵĺ���
		input����coords��forces��vec�ṹ��������p��������
		output����force_torque��<vec, vec>�ṹ������
	*/
	vecp derivative(const vecv& coords, const vecv& forces, flv::iterator& p) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		fl& d = *p; // reference   ȡ��ǰ������ָ���Ԫ��
		++p;        //������ָ����һ��Ԫ��
		branches_derivative(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, d);
		return force_torque;
	}
};

typedef tree<segment> branch;       //�β�Ϊsegment�ṹ���tree�ṹ��
typedef std::vector<branch> branches;  //branch�ṹ������

template<typename Node> // Node == first_segment || rigid_body
struct heterotree {
	Node node;          //�������ͱ���
	branches children;  //branch�ṹ������
	heterotree(const Node& node_) : node(node_) {}
	/*
	    ������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
		input����atoms��atom�ṹ��������coords��vec�ṹ��������c��ligand_conf�ṹ��
	*/
	void set_conf(const atomv& atoms, vecv& coords, const ligand_conf& c) {
		node.set_conf(atoms, coords, c.rigid);
		flv::const_iterator p = c.torsions.begin();   //p������ָ��torsions�����ĵ�һ��Ԫ��
		branches_set_conf(children, node, atoms, coords, p);
		assert(p == c.torsions.end());               //p������ָ��torsions���������һ��Ԫ��
	}

	/* ������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
	 input����atoms��atom�ṹ��������coords��vec�ṹ��������c��residue_conf�ṹ��
	*/
	void set_conf(const atomv& atoms, vecv& coords, const residue_conf& c) {
		flv::const_iterator p = c.torsions.begin();//p������ָ��torsions�����ĵ�һ��Ԫ��
		node.set_conf(atoms, coords, *p);
		++p;
		branches_set_conf(children, node, atoms, coords, p);
		assert(p == c.torsions.end());               //p������ָ��torsions���������һ��Ԫ��
	}

	/*
	input����coords��forces��vec�ṹ��������c��ligand_change�ṹ��
	*/
	void derivative(const vecv& coords, const vecv& forces, ligand_change& c) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		flv::iterator p = c.torsions.begin();             //p������ָ��torsions�����ĵ�һ��Ԫ��
		branches_derivative(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, c.rigid);
		assert(p == c.torsions.end());                    //p������ָ��torsions���������һ��Ԫ��
	}

	/*��Ҫ����node���͵�������ڵĺ���
	input����atoms��atom�ṹ��������coords��vec�ṹ��������c��residue_conf�ṹ��
	*/
	void derivative(const vecv& coords, const vecv& forces, residue_change& c) const {
		vecp force_torque = node.sum_force_and_torque(coords, forces);
		flv::iterator p = c.torsions.begin();
		fl& d = *p; // reference
		++p;
		branches_derivative(children, node.get_origin(), coords, forces, force_torque, p);
		node.set_derivative(force_torque, d);
		assert(p == c.torsions.end());
	}
};

/*Ť�ؼ���
	input����t���������ͣ�s��unint
*/
template<typename T> // T = main_branch, branch, flexible_body
void count_torsions(const T& t, sz& s) {
	t.node.count_torsions(s);
	VINA_FOR_IN(i, t.children)
		count_torsions(t.children[i], s);
}

typedef heterotree<rigid_body> flexible_body;   //�β�Ϊrigid_body�ṹ���heterotree�ṹ��
typedef heterotree<first_segment> main_branch; //�β�Ϊfirst_segment�ṹ���heterotree�ṹ��

/*
	vector_mutable�ṹ��̳�ϵͳvector��

*/
template<typename T> // T == flexible_body || main_branch

struct vector_mutable : public std::vector<T> {
	template<typename C>
	/*������ꡢ����Լ�������Ԫ��orientation_q=q��orientation_m�ṹ���coords
		input���� atoms��atom�ṹ��������coords��vec�ṹ��������c��c���͵�����
	*/
	void set_conf(const atomv& atoms, vecv& coords, const std::vector<C>& c) { // C == ligand_conf || residue_conf
		VINA_FOR_IN(i, (*this))                                //for(i=0;i<v.size;i++)
			(*this)[i].set_conf(atoms, coords, c[i]);
	}
	/*Ť�ؼ���
		output����tmp��unint����
	*/
	szv count_torsions() const {
		szv tmp(this->size(), 0);
		VINA_FOR_IN(i, (*this))
			::count_torsions((*this)[i], tmp[i]);
		return tmp;
	}

	/*��Ҫ����c���͵�������ڵĺ���
		input����coords��forces��vec�ṹ��������c��������������
	
	*/
	template<typename C>
	void derivative(const vecv& coords, const vecv& forces, std::vector<C>& c) const { // C == ligand_change || residue_change
		VINA_FOR_IN(i, (*this))
			(*this)[i].derivative(coords, forces, c[i]);
	}
};

/* ��Ҫ����t��f���͵��ú���         
	input����t��f���������ͱ�����
*/
template<typename T, typename F> // tree or heterotree - like structure
void transform_ranges(T& t, const F& f) {
	t.node.transform(f);
	VINA_FOR_IN(i, t.children)
		transform_ranges(t.children[i], f);
}

#endif
