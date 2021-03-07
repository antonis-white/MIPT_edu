#include<iostream>
#include <math.h>
#include <algorithm>
#include <utility>
#include <ctime>

using std::cerr;

///------------------------DISCLAIMER--------------------------------------------//
//This is implicit AVL-tree, which was made for specific problem

//all comparisons {<, >, <=, ==, ...} of obj_type must be defined

//I'm not sure, that this code will work with other than numeric types, but I don't care

///Note: obj_type() must return neutral element
///-------------------------------------------------------------------------------//
template <typename obj_type = long long>
class impl_AVL_tree {
private:
	static int const V = 130047;  //maximal number of vertices in the system of AVL-trees
	static int const H = 50; //maximal height of tree (and maybe extra for insurance)

	struct vertex {
		//parameters that distinguish assignment and addition
		bool assign = 0;
		//parameter that shows that segment (subtree) must be reversed
		bool rev = 0;
		//height
		short h = 1;
		//number of vertices in tree
		int sz = 1;
		//length of the longest non-decreasing suffix
		int nd_suf_len = 1;
		//length of the longest non-increasing suffix
		int ni_suf_len = 1;
		//length of the longest non-decreasing prefix
		int nd_pref_len = 1;
		//length of the longest non-increasing prefix
		int ni_pref_len = 1;	
		//left child
		vertex* l = nullptr; 
		//right child
		vertex* r = nullptr;
		//object
		obj_type val;
		//sum in subtree (left associative)
		obj_type sum;
		//value to assign/add
		obj_type push = 0;
		//value of leftist element
		obj_type leftist;
		//value of rightist element
		obj_type rightist;
		vertex() = default;
	};
	using node = vertex*;
	
	//return new node with given key
	static node new_node(const obj_type& obj) {
		static vertex* vertices;
		static int used = 0;
		if (!used) 
			vertices = new vertex[V];
		node v = &vertices[used++];
		v->val = obj;     		
		v->sum = obj;
		v->leftist = obj;
		v->rightist = obj;
		return v;		
	}
	
	//methods returning parameters of node (maybe empty)
	
	static int sz(node v) { return v ? v->sz : 0; }

	static short height(node v) { return v ? v->h : 0; }

	static obj_type sum(node v) { return v ? v->sum : 0; }

	static int nd_suf_len(node v) { return v ? v->nd_suf_len : 0; }

	static int ni_suf_len(node v) { return v ? v->ni_suf_len : 0; }

	static int nd_pref_len(node v) { return v ? v->nd_pref_len : 0; }

	static int ni_pref_len(node v) { return v ? v->ni_pref_len : 0; } 
	
	//reverses subtree
	static void assign(node v, const obj_type& val) { 
		if (v) {
			v->assign = 1;
			v->nd_suf_len = v->sz;
			v->ni_suf_len = v->sz;
			v->nd_pref_len = v->sz;
			v->ni_pref_len = v->sz;
			v->val = val;
			v->sum = val * (v->sz);
			v->push = val;
			v->leftist = val;
			v->rightist = val;			
		}
        }	

        static void add(node v, const obj_type& val) {
        	if (v) {        		
        		v->val += val;      
        		v->sum += val * (v->sz);
        		v->push += val;
        		v->leftist += val;
        		v->rightist += val;	
        	}
        }

       	static void reverse(node v) {
       		if (v) {
       			std::swap(v->nd_suf_len, v->ni_pref_len);
			std::swap(v->ni_suf_len, v->nd_pref_len);
			std::swap(v->leftist, v->rightist);
			v->rev ^= 1;
		}			
       	}
	//pushes updates
	static void push(node v) {
		if (!v)	
			return;
		if (v->rev) {
			std::swap(v->l, v->r);
			reverse(v->l);
			reverse(v->r);
			v->rev = 0;		
		}
		if (v->assign) {
			assign(v->l, v->push);
			assign(v->r, v->push);
			v->push = 0;
			v->assign = 0;
			return;
		}
		if (v->push != 0) {
			add(v->l, v->push);
			add(v->r, v->push);
			v->push = 0;	
		}
	}

	///recalculation of interesting prefix/suffixes lengths (assuming v is non-empty)

	static void upd_nd_suf_len(node v) {
		v->nd_suf_len = nd_suf_len(v->r);
		if (v->nd_suf_len != sz(v->r))
			return;
		//checking that element in current vertex continues suffix 
		if (v->r && v->val > v->r->leftist)
			return;
		++(v->nd_suf_len);
		//adding left subtree
		if (v->l && v->l->rightist <= v->val)
			v->nd_suf_len += v->l->nd_suf_len;	
	}

	static void upd_ni_suf_len(node v) {
		v->ni_suf_len = ni_suf_len(v->r);
		if (v->ni_suf_len != sz(v->r))
			return;
		//checking that element in current vertex continues suffix 
		if (v->r && v->val < v->r->leftist)
			return;
		++(v->ni_suf_len);
		//adding left subtree
		if (v->l && v->l->rightist >= v->val)
			v->ni_suf_len += v->l->ni_suf_len;	
	}

	static void upd_nd_pref_len(node v) {
		v->nd_pref_len = nd_pref_len(v->l);
		if (v->nd_pref_len != sz(v->l))
			return;
		//checking that element in current vertex continues prefix 
		if (v->l && v->l->rightist > v->val) 
			return;
		++(v->nd_pref_len);
		//adding right subtree
		if (v->r && v->val <= v->r->leftist)
			v->nd_pref_len += v->r->nd_pref_len;	
	}

	static void upd_ni_pref_len(node v) {
		v->ni_pref_len = ni_pref_len(v->l);
		if (v->ni_pref_len != sz(v->l))
			return;
		//checking that element in current vertex continues prefix 
		if (v->l && v->l->rightist < v->val) 
			return;
		++(v->ni_pref_len);
		//adding right subtree
		if (v->r && v->val >= v->r->leftist)
			v->ni_pref_len += v->r->ni_pref_len;	
	}

	static void upd_fun(node v) {
		upd_nd_suf_len(v);
		upd_ni_suf_len(v);
		upd_nd_pref_len(v);
		upd_ni_pref_len(v);
		v->sum = sum(v->l) + v->val + sum(v->r);
		v->leftist = (v->l ? v->l->leftist : v->val);
		v->rightist = (v->r ? v->r->rightist : v->val);	
	}

	//recalc quantities of tree with root V
	static void upd_node(node v) {
		if (!v)
			return;
		v->h = std::max(height(v->l), height(v->r)) + 1;
		v->sz = sz(v->l) + sz(v->r) + 1;
		upd_fun(v);
	}

	//output in stderr all objects in subtree with root V (in increasing order) 
	static void print_array(node v, std::ostream& out) {
		if (!v) 
			return;
		push(v);	
		print_array(v->l, out);
		out << v->val << ' ';
		print_array(v->r, out);
	}

	//returns true if and only if there is node with given key in tree
	static node search(node v, int key) {
		if (!v)
			return nullptr;
	       	push(v);
		int cur_key = sz(v->l) + 1;
		if (cur_key == key)
			return v;
	       	if (key < cur_key)
	       		return search(v->l, key);
	       	return search(v->r, key - cur_key);	
	}
	                
	//hang tree with root V as right child of TO
	static void hang_right(node v, node to) {
		push(to);
		to->r = v;
		upd_node(to);	
	}

	//hang tree with root V as left child of TO
	static void hang_left(node v, node to) {
		push(to);
		to->l = v;
		upd_node(to);	
	}

	//cut left subtree of node v
	static node cut_left(node v) {
		if (!v)
			return nullptr;
		push(v);
		node res = v->l;
		v->l = nullptr;
		upd_node(v);
		return res;
	}

	//cut left subtree of node v
	static node cut_right(node v) {
		if (!v)
			return nullptr;
		push(v);
		node res = v->r;
		v->r = nullptr;
		upd_node(v);
		return res;
	}

	//do small left rotation in tree
	static void small_left_rotation(node& root) {
		//new root
		node new_root = root->r;   
		push(new_root);
		//hanging middle subtree to the previous root
		hang_right(new_root->l, root);
		//hanging previous root to the new one   
		hang_left(root, new_root);
		//changing root
		root = new_root;		
	}

	//do small right rotation in tree  
	static void small_right_rotation(node& root) {
		//new root
		node new_root = root->l;
		push(new_root);                                                                   	
		//hang middle subtree to the previous root 
		hang_left(new_root->r, root);
		//hang previous root to the new one     
		hang_right(root, new_root);         
		//changing root
		root = new_root;	
	}

	//do big left rotate in tree
	static void big_left_rotation(node& root) {
		push(root->r);
		small_right_rotation(root->r);
		small_left_rotation(root);				
	}

	//do big right rotation in tree 
	static void big_right_rotation(node& root) {
		push(root->l);
		small_left_rotation(root->l);
		small_right_rotation(root);			
	}       	

	//rebalance subtree (balance factor of V lies in [-2; 2]) 
	static void rebalance(node& root) {
		if (!root)
			return;	  
		push(root);
		//balance factor
		short dif = height(root->r) - height(root->l);
		if (abs(dif) != 2)
			return;
		//right-heavy subtree       	
		if (dif == 2) {
			push(root->r);
			if (height(root->r->l) <= height(root->r->r))
				small_left_rotation(root);
			else
				big_left_rotation(root);			
		} else {
			push(root->l);
			if (height(root->l->r) <= height(root->l->l))
				small_right_rotation(root);
			else
				big_right_rotation(root);			 	
		}	
	}

	//insert vertex V in tree after key-indexed node
	//don't store identical keys in tree in this realisation
	static void insert(node& root, int key, node v) {
		//position for new vertex found
		if (!root) {
			root = v;
			return;
		}
		push(root);
		int cur_key = sz(root->l) + 1;
		if (key < cur_key)   //V should be in left subtree	
			insert(root->l, key, v);
		else               //V should be in right subtree left subtree
			insert(root->r, key - cur_key, v);
		upd_node(root);
		rebalance(root);
	}       	

	//extract node with maximal key from tree (nullptr from empty tree)
	static node extract_last(node& root) {
		node ans = nullptr;
		//empty tree
		if (!root)
			return ans;
		push(root);	
		if (root->r) {   //there is bigger key in right subtree
			ans = extract_last(root->r);
			upd_node(root);
			rebalance(root);
		} else {       //found maximal key
			ans = root;
			//extract node and hang its left subtree back
			root = cut_left(ans);	
		}	
		return ans;
	}

	//merge 2 trees and one node, where keys are stricly increasing throw this system of trees
	//returns combined tree
	static node big_merge(node l_root, node mid, node r_root) {
		//diffrence of heights
		short dif = height(r_root) - height(l_root);
		//it's possible to use mid as the root if combined tree
		if (abs(dif) < 2) {
			hang_left(l_root, mid);
			hang_right(r_root, mid);
			return mid;
		}
		push(l_root);
		push(r_root);
		if (dif > 0) { ///right tree is heavier
			//make root of right tree as root of combined tree
			hang_left(big_merge(l_root, mid, cut_left(r_root)), r_root);
			rebalance(r_root);
			return r_root;
		}
		///left tree is heavier
		//make root of left tree as root of combined tree
		hang_right(big_merge(cut_right(l_root), mid, r_root), l_root);
		rebalance(l_root);
		return l_root;              
	} 

	//concatenate two trees into one
	static void merge(node& root, node l_root, node r_root) {
		if (!l_root || !r_root) {
			root = (r_root ? r_root : l_root);  
			return;
		}
		node mid = extract_last(l_root);
		root = big_merge(l_root, mid, r_root);
	}

	struct fast_stack {
		std::pair<node, node> st[H]; //contains pairs {tree, vertex}
		int it = -1;
		bool left = 1; //orientation of merge
		fast_stack(bool left) : left(left) {} 
		//add pair {tree, vertex} on top of stack 
		void add(std::pair<node, node> q) {
			///instantly merges trees with same height
			while (it != -1) {
				auto [lst, mid] = st[it];
				if (height(lst) != height(q.first))
					break;
				--it;
				if (left)
					q.first = big_merge(lst, mid, q.first);	
				else
					q.first = big_merge(q.first, mid, lst);
			}
			++it;
			st[it] = q;
		}
		//combines all trees (last tree is given)
		node combine(node res) {	    
			while (it != -1) {
				if (left)
					res = big_merge(st[it].first, st[it].second, res);
				else
					res = big_merge(res, st[it].second, st[it].first);
				--it;		
			}
			return res;
		}
	};

	//divides one tree into two, where left one has first key elements of given tree
	static void split(node root, node& l_root, node& r_root, int key) {
		///fast stacks to merge subtrees of resulted trees
		static fast_stack l_stack(1), r_stack(0);
		while (root) {
			push(root);
			int cur_key = 1 + sz(root->l);
			if (cur_key == key)
				break;
			node l = cut_left(root);
			node r = cut_right(root);
			if (cur_key < key) {			
				l_stack.add({l, root});	
				root = r;
				key -= cur_key;
			} else {
				r_stack.add({r, root});
				root = l;				
			}											
		}
		l_root = l_stack.combine(cut_left(root));
		r_root = r_stack.combine(cut_right(root));
		if (root)
			merge(l_root, l_root, root);
	}       	

	//erase node with given key from tree (if there is one)
	static void erase(node& root, int key) {
		node l_part, r_part;
		split(root, l_part, r_part, key);
		//there isn't node with given key
		if (!l_part) 
			return;
		//extract vertex with given key
		extract_last(l_part);
		merge(root, l_part, r_part);		
	}

	//does next permutation in subtree
	static void next_permutation(node& root) {	
		if (root->ni_suf_len == root->sz) {
			reverse(root);
			return;
		}      
		node suf_part;
		node pivot;
		split(root, root, suf_part, sz(root) - root->ni_suf_len - 1);
		split(suf_part, pivot, suf_part, 1);
		//divided into: root, pivot, suf_part
		reverse(suf_part);
		//count number of element in suffix with <= value relative to the value of pivot element
		//NOTE: valus in suf_part are non_decreasing 
		int less_eq = 0;
		node v = suf_part;
		while (v) {
			push(v);
			if (v->val > pivot->val)
				v = v->l;
		       	else
		       		less_eq += 1 + sz(v->l), v = v->r;
		}
		node suf_l_part;
		node suf_pivot;
		node suf_r_part;
		split(suf_part, suf_l_part, suf_r_part, less_eq);
		split(suf_r_part, suf_pivot, suf_r_part, 1);
		//divided into: root, pivot, suf_l_part, suf_pivot, suf_r_part
		//swap pivot and suf_pivot and merge all back
		std::swap(pivot, suf_pivot);
		suf_part = big_merge(suf_l_part, suf_pivot, suf_r_part);
		root = big_merge(root, pivot, suf_part); 
	}

	//does next permutation in subtree
	static void prev_permutation(node& root) {	
		if (root->nd_suf_len == root->sz) {
			reverse(root);
			return;
		}
		node suf_part;
		node pivot;
		split(root, root, suf_part, sz(root) - root->nd_suf_len - 1);
		split(suf_part, pivot, suf_part, 1);
		//divided into: root, pivot, suf_part
		reverse(suf_part);
		//count number of element in suffix with greater value relative to the value of pivot element
		//NOTE: valus in suf_part are non_increasing 
		int gre_eq = 0;
		node v = suf_part;
		while (v) {
			push(v);
			if (v->val < pivot->val)
				v = v->l;
		       	else
		       		gre_eq += 1 + sz(v->l), v = v->r;
		}
		node suf_l_part;
		node suf_pivot;
		node suf_r_part;
		split(suf_part, suf_l_part, suf_r_part, gre_eq);
		split(suf_r_part, suf_pivot, suf_r_part, 1);
		//divided into: root, pivot, suf_l_part, suf_pivot, suf_r_part
		//swap pivot and suf_pivot and merge all back
		std::swap(pivot, suf_pivot);
		suf_part = big_merge(suf_l_part, suf_pivot, suf_r_part);
		root = big_merge(root, pivot, suf_part); 
	}
///end of implementation
private:
	node root = nullptr;
	impl_AVL_tree(node root) : root(root) {}
public:
	impl_AVL_tree() = default;
	impl_AVL_tree(const impl_AVL_tree& t) : root(t.root) {}
	//insert vertex with given key
	//don't store identical keys in tree in this realisation
	void insert(int key, const obj_type& obj) { insert(root, key, new_node(obj)); }
	//erase node with given key from tree (if there is one)
	void erase(int key) { erase(root, key + 1); }
	//returns true if and only if there is node with given key in tree
	obj_type find_by_order(int key) {
		node v = search(root, key + 1);
		return v ? v->val : obj_type();
	}
	//print all keys in ostream
	friend std::ostream& operator<<(std::ostream& out, impl_AVL_tree t) { print_array(t.root, out); return out; }
	template <typename k>
	friend void merge(impl_AVL_tree<k>& t, impl_AVL_tree<k> l, impl_AVL_tree<k> r);
	template <typename k>
	friend void split(impl_AVL_tree<k> t, impl_AVL_tree<k>& l, impl_AVL_tree<k>& r, int key);
	void clear() { root = nullptr; };
	int size() { return sz(root); }
	bool empty() { return !root; };
	void push_front(const obj_type& obj) { insert(0, obj); }
	void push_back(const obj_type& obj) { insert(size(), obj); }
	void print(std::ostream& out) const { print_array(root, out); }
	//return minimum on segment
	obj_type sum_on_seg(int ql, int qr) {
		++ql, ++qr; //change of indexation
		node l_part;
		node r_part;
		split(root, l_part, r_part, ql - 1);
		node seg;
		split(r_part, seg, r_part, qr - ql + 1);
		obj_type res = seg->sum;
		merge(r_part, seg, r_part);
		merge(root, l_part, r_part);
		return res;
	}

	void add_on_seg(int ql, int qr, const obj_type& val) {
		++ql, ++qr; //change of indexation
		node l_part;
		node r_part;
		split(root, l_part, r_part, ql - 1);
		node seg;
		split(r_part, seg, r_part, qr - ql + 1);
		add(seg, val);
		merge(r_part, seg, r_part);
		merge(root, l_part, r_part);				
	}

	void assign_on_seg(int ql, int qr, const obj_type& val) {
		++ql, ++qr; //change of indexation
		node l_part;
		node r_part;
		split(root, l_part, r_part, ql - 1);
		node seg;
		split(r_part, seg, r_part, qr - ql + 1);
		assign(seg, val);
		merge(r_part, seg, r_part);
		merge(root, l_part, r_part);				
	}

	void next_perm_on_seg(int ql, int qr) {
		++ql, ++qr; //change of indexation
		node l_part;
		node r_part;
		split(root, l_part, r_part, ql - 1);
		node seg;
		split(r_part, seg, r_part, qr - ql + 1);
		next_permutation(seg);
		merge(r_part, seg, r_part);
		merge(root, l_part, r_part);				
	}

	void prev_perm_on_seg(int ql, int qr) {
		++ql, ++qr; //change of indexation
		node l_part;
		node r_part;
		split(root, l_part, r_part, ql - 1);
		node seg;
		split(r_part, seg, r_part, qr - ql + 1);
		prev_permutation(seg);
		merge(r_part, seg, r_part);
		merge(root, l_part, r_part);				
	}
};

//concatenate two trees into one
template <typename obj_type> 
void merge(impl_AVL_tree<obj_type>& t, impl_AVL_tree<obj_type> l, impl_AVL_tree<obj_type> r) { 
	impl_AVL_tree<obj_type>::merge(t.root, l.root, r.root); 
}

//divides one tree into two, where left one has first key elements of given tree
template <typename obj_type> 
void split(impl_AVL_tree<obj_type> t, impl_AVL_tree<obj_type>& l, impl_AVL_tree<obj_type>& r, int key) {	
	impl_AVL_tree<obj_type>::split(t.root, l.root, r.root, key); 
} 

int main(){
	std::ios_base::sync_with_stdio(false);
	std::cin.tie(NULL);
	#ifdef LOCAL
	freopen("in", "r", stdin);
	freopen("out", "w", stdout);
	#endif
	int n, q;
	std::cin >> n;
	impl_AVL_tree tree;
	while (n--) {
		int a;
		std::cin >> a;
		tree.push_back(a);
	}

	std::cin >> q;
	while (q--) {
		int t, l, r, x, i;
		std::cin >> t;

		switch (t) {
		
		case 2: 
			std::cin >> x >> i;
			tree.insert(i, x);
			break;
		case 3: 
			std::cin >> i;
			tree.erase(i);
			break;
		case 4:	
			std::cin >> x >> l >> r;
			tree.assign_on_seg(l, r, x);
			break;
		case 5:
			std::cin >> x >> l >> r;
			tree.add_on_seg(l, r, x);
			break;
		case 6:
			std::cin >> l >> r;
			tree.next_perm_on_seg(l, r);
			break;
		case 7:
			std::cin >> l >> r;
			tree.prev_perm_on_seg(l, r);
			break;
		default:
			std::cin >> l >> r;
			std::cout << tree.sum_on_seg(l, r) << '\n';
			break;
		}
	}
	tree.print(std::cout);
	return 0;
}       		
