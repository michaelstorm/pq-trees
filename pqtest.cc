// This file is part of the PQ Tree library.
//
// The PQ Tree library is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by the 
// Free Software Foundation, either version 3 of the License, or (at your
// option) any later version.
//
// The PQ Tree Library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License along 
// with the PQ Tree Library.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <set>
#include "pqnode.h"
#include "pqtree.h"

template <typename T>
string ReadableType(typename PQNode<T>::PQNode_types type) {
  if (type == PQNode<T>::leaf) {
    return "leaf";
  } else if (type == PQNode<T>::pnode) {
    return "P-Node";
  } else if (type == PQNode<T>::qnode) {
    return "Q-Node";
  }
  return "unknown";
}

void testbed() {
  set<int> S;
  for(int i=1;i<9;i++)
    S.insert(i);
  PQTree<int> P(S);

  cout << "PQ Tree with 8 elements and no reductions" << endl;  
  cout << P.Print() << endl;
  
  cout << "Reducing by set {4, 3}" << endl;
  S.clear();
  S.insert(4);
  S.insert(3);
  cout << P.Reduce(S) << endl;
  cout << P.Print() << endl;
  cout << "Permutations: " << P.PermutationCount() << endl;
  
  cout << "Reducing by set {4, 3, 6}" << endl;
  S.clear();
  S.insert(6);
  S.insert(4);
  S.insert(3);
  cout << P.Reduce(S) << endl;
  cout << P.Print() << endl;
  cout << "Permutations: " << P.PermutationCount() << endl;
  
  cout << "Reducing by set {4, 3, 5}" << endl;
  S.clear();
  S.insert(4);
  S.insert(3);
  S.insert(5);
  cout << P.Reduce(S) << endl;
  cout << P.Print() << endl;
  cout << "Permutations: " << P.PermutationCount() << endl;
    
  cout << "Reducing by set {4, 5}" << endl;
  S.clear();
  S.insert(4);
  S.insert(5);
  cout << P.Reduce(S) << endl;
  cout << P.Print() << endl;
  cout << "Permutations: " << P.PermutationCount() << endl;
  
  cout << "Reducing by set {2, 6}" << endl;
  S.clear();
  S.insert(2);
  S.insert(6);
  cout << P.Reduce(S) << endl;
  cout << P.Print() << endl;
  cout << "Permutations: " << P.PermutationCount() << endl;
  
  cout << "Reducing by set {1, 2}" << endl;
  S.clear();
  S.insert(1);
  S.insert(2);
  cout << P.Reduce(S) << endl;
  cout << P.Print() << endl;
  cout << "Permutations: " << P.PermutationCount() << endl;
  
  cout << "Reducing by set {4, 5}" << endl;
  S.clear();
  S.insert(4);
  S.insert(5);
  cout << P.Reduce(S) << endl;
  cout << P.Print() << endl;
  cout << "Permutations: " << P.PermutationCount() << endl;

  // Lets actually explore the tree manually
  cout << endl;
  PQNode<int>* root = P.Root();
  cout << "Root Type: " << ReadableType<int>(root->Type()) << endl;
  vector<PQNode<int>*> children;
  root->Children(&children);
  for (int i = 0; i < children.size(); ++i) {
    PQNode<int>* child = children[i];
    cout << "Child " << i+1 << " Type: " << ReadableType<int>(child->Type());
    if (child->Type() == PQNode<int>::leaf) {
      cout << " Value: " << child->LeafValue() << endl;
    } else {
      cout << endl;
      vector<PQNode<int>*> grandchildren;
      child->Children(&grandchildren);
      for (int j = 0; j < grandchildren.size(); ++j) {
	PQNode<int>* grandchild = grandchildren[j];
        cout << "GrandChild " << j+1 << " Type: " 
	     << ReadableType<int>(grandchild->Type());
        if (grandchild->Type() == PQNode<int>::leaf)
	  cout << " Value: " << grandchild->LeafValue();
	cout << endl;
      }
    }
  }
  cout << endl;
 
  // Now, we perform a reduction that will fail. 
  cout << "Reducing by set {5, 3} - will fail" << endl;
  S.clear();
  S.insert(5);
  S.insert(3);
  cout << P.Reduce(S) << endl;
  cout << P.Print() << endl;
}
  

int main(int argc, char **argv) {
  testbed();
}


