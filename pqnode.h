// Internal class used by pqtree to represent individual nodes in the pqtree.

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

#ifndef PQNODE_H
#define PQNODE_H

#include <assert.h>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <cstdio>
using namespace std;

template <typename T> class PQTree;
template <typename T> class QNodeChildrenIterator;

template <typename T>
class PQNode {
 // PQNodes are not exposed by pqtrees, they are internally used only.
 friend class PQTree<T>;
  friend class QNodeChildrenIterator<T>;

 public:
  // Enum types we use throughout.
  enum PQNode_types  {leaf, pnode, qnode};
  enum PQNode_marks  {unmarked, queued, blocked, unblocked};
  enum PQNode_labels {empty, full, partial};

  // Returns the type of the current node, an enum of type PQNode_types.
  PQNode_types Type();

  // Returns the value of the leaf node.  Fails assertion if not leaf node.
  T LeafValue();

  // Returns all of this Node's children if it has any.
  // Return Value is the |children| argument.
  void Children(vector<PQNode<T>*> *children);
  
  int PermutationCount();

 private:
  /***** Used by P Nodes only *****/
  
  // A doubly-linked of links which form the children of a p-node, the order
  // of the list is arbitrary.
  list<PQNode<T>*> circular_link_;
  
  // A count of the number of children used by a node.
  int ChildCount();

  // Returns the first |circular_link_| child with a given label or NULL.
  PQNode<T>* CircularChildWithLabel(PQNode_labels label);

  // Moves the full children of this node to children of |new_node|.
  void MoveFullChildren(PQNode<T>* new_node);
  
  // Replaces the circular_link pointer of |old_child| with |new_child|.
  void ReplaceCircularLink(PQNode<T>* old_child, PQNode<T>* new_child);
  
  /***** Used by Q Nodes only *****/
  
  // A set containing the two endmost children of a Q-node
  PQNode<T> *endmost_children_[2];
  PQNode<T> *pseudo_neighbors_[2];

  // Boolean indicating whether or not this is a pseudonode/child.
  bool pseudonode_, pseudochild_;

  // Returns the first endmost child with a given label or NULL.
  PQNode<T>* EndmostChildWithLabel(PQNode_labels label);
  
  // Returns the first immediate sibling with a given label or NULL.
  PQNode<T>* ImmediateSiblingWithLabel(PQNode_labels label);

  // Returns the first immediate sibling without a given label or NULL.
  PQNode<T>* ImmediateSiblingWithoutLabel(PQNode_labels label);

  // Adds an immediate sibling to this node.
  void AddImmediateSibling(PQNode<T>* sibling);
  
  // Adds an immediate sibling to this node.
  void RemoveImmediateSibling(PQNode<T>* sibling);

  // Nulls otu the immediate siblings of this node.
  void ClearImmediateSiblings();

  // Returns the number of immediate siblings this node has (0, 1, or 2).
  int ImmediateSiblingCount();

  // Replaces the |endmost_children_| pointer to |old_child| with |new_child|.
  void ReplaceEndmostChild(PQNode<T>* old_child, PQNode<T>* new_child);
  
  // Replaces the immediate sibling of |old_child| with |new_child|.
  void ReplaceImmediateSibling(PQNode<T>* old_child, PQNode<T>* new_child);

  // Replaces the partial child |old_child| with |new_child|.
  void ReplacePartialChild(PQNode<T>* old_child, PQNode<T>* new_child);

  // Forces a Q-Node to "forget" it's pointers to it's endmost children.  
  // Useful if you want to delete a Q-Node but not it's children.
  void ForgetChildren();

  // Returns true if all of the full and partial children of this node are
  // consecutive, with the partial children on the outside.
  bool ConsecutiveFullPartialChildren();
  
  /***** Used by Both Node types *****/
  
  // A set containing all the children of a node currently known to be full.
  set<PQNode<T>*> full_children_;
  
  // A set containing all the children of a node currently known to be partial.
  set<PQNode<T>*> partial_children_;
  
  // Only children of Q nodes have more than 0 immediate siblings.  Stores the
  // siblings to either side of this node in its parent's children.  One or both
  // siblings may be NULL.  
  PQNode<T> *immediate_siblings_[2];

  // Label is an indication of whether the node is empty, full, or partial
  enum PQNode_labels label_;
  
  // Mark is a designation used during the first pass of the reduction
  // algorithm.  Every node is initially unmarked.  It is marked
  // queued when it is placed onto the queue during the bubbling up.
  // It is marked either blocked or unblocked when it is processed.  Blocked
  // nodes can become unblocked if their siblings become unblocked.
  enum PQNode_marks mark_;

  // type is a designation telling whether the node is a leaf, P, or Q.
  enum PQNode_types type_;

  // the immediate ancestor of a node.  This field is always
  // valid for children of P-nodes and for endmost children of Q-nodes
  PQNode<T>* parent_;

  // A count of the number of pertinent children currently possessed by a node
  int pertinent_child_count;

  // A count of the number of pertinent leaves currently possessed by a node
  int pertinent_leaf_count;

  // The value of the PQNode<T> if it is a leaf.
  int leaf_value_;

  // Makes a deep copy of a node, sets this to be it's parent and returns copy.
  PQNode<T>* CopyAsChild(const PQNode<T>& to_copy);

  // Makes a deep copy of copy.
  void Copy(const PQNode<T>& copy);
  
  // deep assignment operator
  PQNode<T>& operator=(const PQNode<T>& to_copy);
  
  // Return the next child in the immediate_siblings chain given a last pointer
  // if last pointer is null, will return the first sibling.  Behavior similar
  // to an iterator.
  PQNode<T>* QNextChild(PQNode<T> *last);

  void ReplaceChild(PQNode<T>* old_child, PQNode<T>* new_child);
  
  // removes this node from a q-parent and puts toInsert in it's place
  void SwapQ(PQNode<T> *toInsert);
  
  // deep copy constructor
  PQNode<T>(const PQNode<T>& to_copy);
    
  // Constructor for a leaf PQNode<T>.
  PQNode<T>(T value);
  
  // Constructor for non-leaf PQNode<T>.
  PQNode<T>();
  
  // Deep destructor.
  ~PQNode<T>();
  
  // Label's this node as full, updating the parent if needed.
  void LabelAsFull();

  // Walks the tree to build a map from values to leaf pointers.
  void FindLeaves(map<T, PQNode<T>*> &leafAddress);
  
  // Walks the tree to find it's Frontier, returns one possible ordering.
  void FindFrontier(list<T> &ordering);
  
  // Resets a bunch of temporary variables after the reduce walks
  void Reset();
  
  // Walks the tree and prints it's structure to |out|.  P-nodes are
  // represented as ( ), Q-nodes as [ ], and leafs by their integer id. Used
  // primarily for debugging purposes
  void Print(string *out);
};

// Q-Nodes have an unusual structure that makes iterating over their children
// slightly tricky.  This class makes the iterating much simpler.
//
// A Q-Node has two children whose ordering is defined strictly, but can be
// reversed, ie: (1,2,3,4) or (4,3,2,1) but not (2,1,4,3).  The Q-Node itself
// has only unordered pointers to the 1 or 2 children nodes on the ends of this
// list.  Each child contains pointers to the 1 or two siblings on either side
// of itself, but knows nothing of the ordering.
//
// Usage:
//   for (QNodeChildrenIterator it(candidate_node); !it.IsDone(); it.Next()) {
//     Process(it.Current()); 
//   }
template <typename T>
class QNodeChildrenIterator {
 public:
  // Creates an iterator of the children of |parent| optionally forcing the
  // iteration to start on the |begin_side|.
  QNodeChildrenIterator(PQNode<T>* parent, PQNode<T>* begin_side=NULL);

  // Returns a pointer to the current PQNode<T> child in the list.
  PQNode<T>* Current();

  // Next advances the current position in the child list.  If |IsDone()|,
  // operation has no effect.
  void Next();

  // Resets the iterator to the beginning of the list of children, optionally
  // forcing the iteration to start on the |begin_side|.
  void Reset(PQNode<T>* begin_side=NULL);

  // Indicate whether or not all children have been looped through.  The initial
  // list order (forward or reverse) is consistent during the life of the
  // object.
  bool IsDone();

 private:
  // Next() helper method to deal with pseudonodes.
  void NextPseudoNodeSibling();
  PQNode<T>* parent_;
  PQNode<T>* current_;
  PQNode<T>* next_;
  PQNode<T>* prev_;
};

template <typename T>
typename PQNode<T>::PQNode_types PQNode<T>::Type() {
  return type_;
}

template <typename T>
T PQNode<T>::LeafValue() {
  assert(type_ == leaf);
  return leaf_value_;
}

template <typename T>
int PQNode<T>::PermutationCount() {
  int count = 1;
  if (type_ != leaf) {
    for (typename list<PQNode<T>*>::const_iterator j = circular_link_.begin();
         j != circular_link_.end(); ++j)
      count *= (*j)->PermutationCount();
  }
  
  if (type_ == pnode) {
    for (int i = circular_link_.size(); i > 0; i--)
      count *= i;
  } else if (type_ == qnode) {
    int q_children = 0;
    for(QNodeChildrenIterator<T> qit(this); !qit.IsDone(); qit.Next()) {
      q_children++;
      if (q_children > 1) {
        count *= 2;
        break;
      }
    }
  }
  return count;
}

template <typename T>
void PQNode<T>::Children(vector<PQNode<T>*> *children) {
  assert(children->empty());
  if (type_ == pnode) {
    for (typename list<PQNode<T>*>::const_iterator i = circular_link_.begin();
	 i != circular_link_.end(); ++i)
      children->push_back(*i);
  } else if (type_ == qnode) {
    for(QNodeChildrenIterator<T> qit(this); !qit.IsDone(); qit.Next())
      children->push_back(qit.Current());
  }
} 

template <typename T>
int PQNode<T>::ChildCount() {
  return circular_link_.size();
}

template <typename T>
PQNode<T>* PQNode<T>::CopyAsChild(const PQNode<T>& to_copy) {
  PQNode<T>* temp = new PQNode<T>(to_copy);
  temp->parent_ = this;
  return temp;
}

template <typename T>
PQNode<T>::PQNode(const PQNode<T>& to_copy) {
  Copy(to_copy);
}

template <typename T>
void PQNode<T>::Copy(const PQNode<T>& to_copy) {
  // Copy the easy stuff
  leaf_value_            = to_copy.leaf_value_;
  pertinent_leaf_count  = to_copy.pertinent_leaf_count;
  pertinent_child_count = to_copy.pertinent_child_count;
  type_                  = to_copy.type_;
  mark_                  = to_copy.mark_;
  label_                 = to_copy.label_;
  pseudonode_            = to_copy.pseudonode_;
  pseudochild_           = to_copy.pseudochild_;

  // Make sure that these are unset initially
  parent_ = NULL;
  partial_children_.clear();
  full_children_.clear();
  circular_link_.clear();
  ClearImmediateSiblings();
  ForgetChildren();

  // Copy the nodes in circular link for pnodes.
  // If it is not a pnode, it will be empty, so this will be a no-op.
  for (typename list<PQNode<T>*>::const_iterator i = to_copy.circular_link_.begin();
      i != to_copy.circular_link_.end(); i++)
    circular_link_.push_back(CopyAsChild(**i));

  // Copy the sibling chain for qnodes
  if (type_ == qnode) {
    PQNode<T> *current, *last;
    // Pointers to nodes we are going to copy
    PQNode<T> *curCopy, *lastCopy, *nextCopy;  
    endmost_children_[0] = CopyAsChild(*to_copy.endmost_children_[0]);
    curCopy = to_copy.endmost_children_[0];
    lastCopy = NULL;
    last = endmost_children_[0];
    
    // Get all the intermediate children
    nextCopy = curCopy->QNextChild(lastCopy);
    while (nextCopy != NULL) {
      lastCopy = curCopy;
      curCopy  = nextCopy;
      current  = CopyAsChild(*curCopy);
      current->AddImmediateSibling(last);
      last->AddImmediateSibling(current);
      last = current;
      nextCopy = curCopy->QNextChild(lastCopy);
    }

    // Now set our last endmost_children_ pointer to our last child
    endmost_children_[1] = current;
  }
}

template <typename T>
PQNode<T>& PQNode<T>::operator=(const PQNode<T>& to_copy) {
  // Check for self copy
  if (&to_copy == this)
    return *this;
  Copy(to_copy);
  // Return ourself for chaining.
  return *this;
}

template <typename T>
void PQNode<T>::LabelAsFull() {
  label_ = full;
  if (parent_)
    parent_->full_children_.insert(this);
}

// Return the next child in the immediate_siblings_ chain given a last pointer.
// If last pointer is null, will return the first sibling.
template <typename T>
PQNode<T>* PQNode<T>::QNextChild(PQNode<T> *last) {
  if (immediate_siblings_[0] == last) {
      return immediate_siblings_[1];
  } else {
    if (!last && 2 == ImmediateSiblingCount()) // occurs at edge of pseudonode.
        return immediate_siblings_[1];
    return immediate_siblings_[0];
  }
}

template <typename T>
void PQNode<T>::ReplaceChild(PQNode<T>* old_child, PQNode<T>* new_child) {
  if (type_ == pnode) {
    ReplaceCircularLink(old_child, new_child);
  } else { // qnode
    for (int i = 0; i < 2 && old_child->immediate_siblings_[i]; ++i) {
      PQNode<T> *sibling = old_child->immediate_siblings_[i];
      sibling->ReplaceImmediateSibling(old_child, new_child);
    }
    ReplaceEndmostChild(old_child, new_child);
  }
  new_child->parent_ = old_child->parent_;
  if (new_child->label_ == partial)
    new_child->parent_->partial_children_.insert(new_child);
  if (new_child->label_ == full)
    new_child->parent_->full_children_.insert(new_child);
}
  
// Removes this node from a q-parent and puts toInsert in it's place
template <typename T>
void PQNode<T>::SwapQ(PQNode<T> *toInsert) {
  toInsert->pseudochild_ = pseudochild_;
  toInsert->ClearImmediateSiblings();
  for (int i = 0; i < 2; ++i) {
    if (parent_->endmost_children_[i] == this)
      parent_->endmost_children_[i] = toInsert;
    if (immediate_siblings_[i])
      immediate_siblings_[i]->ReplaceImmediateSibling(this, toInsert);
  }
  ClearImmediateSiblings();
  parent_ = NULL;
}
  
template <typename T>
PQNode<T>::PQNode(T value) {
  leaf_value_             = value;
  type_                  = leaf;
  parent_                = NULL;
  label_                 = empty;
  mark_                  = unmarked;
  pertinent_child_count  = 0;
  pertinent_leaf_count   = 0;
  ClearImmediateSiblings();
  ForgetChildren();
}
  
template <typename T>
PQNode<T>::PQNode() {
  pseudonode_ = false;
  parent_ = NULL;
  label_ = empty;
  mark_ = unmarked;
  pertinent_child_count = 0;
  pertinent_leaf_count = 0;
  ClearImmediateSiblings();
  ForgetChildren();
}

template <typename T>
PQNode<T>::~PQNode<T>() {
  if (type_ == qnode) {
    PQNode<T> *last     = NULL;
    PQNode<T> *current  = endmost_children_[0];
    while(current) {
      PQNode<T> *next = current->QNextChild(last);
      delete last;
      last    = current;
      current = next;
    }
    delete last;
  } else if (type_ == pnode) {  
    for (typename list<PQNode<T>*>::iterator i = circular_link_.begin();
         i != circular_link_.end(); i++)
      delete *i;
    circular_link_.clear();
  }
}

template <typename T>
PQNode<T>* PQNode<T>::CircularChildWithLabel(PQNode_labels label) {
  for (typename list<PQNode<T>*>::iterator i = circular_link_.begin();
       i != circular_link_.end(); i++) {
    if ((*i)->label_ == label)
      return *i;
  }
  return NULL;
}


template <typename T>
PQNode<T>* PQNode<T>::EndmostChildWithLabel(PQNode_labels label) {
  for (int i = 0; i < 2; ++i)
    if (endmost_children_[i] && endmost_children_[i]->label_ == label)
      return endmost_children_[i];
  return NULL;
}

template <typename T>
PQNode<T>* PQNode<T>::ImmediateSiblingWithLabel(PQNode_labels label) {
  for (int i = 0; i < 2 && immediate_siblings_[i]; ++i)
    if (immediate_siblings_[i]->label_ == label)
      return immediate_siblings_[i];
  return NULL;
}

template <typename T>
PQNode<T>* PQNode<T>::ImmediateSiblingWithoutLabel(PQNode_labels label) {
  for (int i = 0; i < 2 && immediate_siblings_[i]; ++i)
    if (immediate_siblings_[i]->label_ != label)
      return immediate_siblings_[i];
  return NULL;
}

template <typename T>
void PQNode<T>::AddImmediateSibling(PQNode<T> *sibling) {
  int null_idx = ImmediateSiblingCount();
  assert(null_idx < 2);
  immediate_siblings_[null_idx] = sibling;
}

template <typename T>
void PQNode<T>::RemoveImmediateSibling(PQNode<T> *sibling) {
  if (immediate_siblings_[0] == sibling) {
    immediate_siblings_[0] = immediate_siblings_[1];
    immediate_siblings_[1] = NULL;
  } else if (immediate_siblings_[1] == sibling) {
    immediate_siblings_[1] = NULL;
  } else {
    assert(false);
  }
}

template <typename T>
void PQNode<T>::ClearImmediateSiblings() {
  for (int i = 0; i < 2; ++i)
    immediate_siblings_[i] = NULL;
}

template <typename T>
int PQNode<T>::ImmediateSiblingCount() {
  int count = 0;
  for (int i = 0; i < 2 && immediate_siblings_[i]; ++i)
    count ++;
  return count;
}

template <typename T>
void PQNode<T>::ReplaceEndmostChild(PQNode<T>* old_child, PQNode<T>* new_child) {
  for (int i = 0; i < 2; ++i) {
    if (endmost_children_[i] == old_child) {
      endmost_children_[i] = new_child;
      return;
    }
  }
}

template <typename T>
void PQNode<T>::ReplaceImmediateSibling(PQNode<T>* old_child, PQNode<T>* new_child) {
  for (int i = 0; i < 2 && immediate_siblings_[i]; ++i)
    if (immediate_siblings_[i] == old_child)
      immediate_siblings_[i] = new_child;
  new_child->immediate_siblings_[new_child->ImmediateSiblingCount()] = this;
}

template <typename T>
void PQNode<T>::ReplacePartialChild(PQNode<T>* old_child, PQNode<T>* new_child) {
  new_child->parent_ = this;
  partial_children_.insert(new_child);
  partial_children_.erase(old_child);
  if (type_ == pnode) {
    circular_link_.remove(old_child);
    circular_link_.push_back(new_child);
  } else {
    old_child->SwapQ(new_child);
  }
}

template <typename T>
void PQNode<T>::ForgetChildren() {
  for (int i = 0; i < 2; ++i)
    endmost_children_[i] = NULL;
}  

template <typename T>
bool PQNode<T>::ConsecutiveFullPartialChildren() {
  // Trivial Case:
  if (full_children_.size() + partial_children_.size() <= 1)
    return true;
  // The strategy here is to count the number of each label of the siblings of
  // all of the full and partial children and see if the counts are correct.
  map<PQNode_labels, int> counts;
  for(typename set<PQNode<T>*>::iterator it = full_children_.begin();
      it != full_children_.end(); ++it) {
    for (int i = 0; i < 2 && (*it)->immediate_siblings_[i]; ++i)
      counts[(*it)->immediate_siblings_[i]->label_]++;
  }
  for(typename set<PQNode<T>*>::iterator it = partial_children_.begin();
      it != partial_children_.end(); ++it) {
    for (int i = 0; i < 2 && (*it)->immediate_siblings_[i]; ++i)
      counts[(*it)->immediate_siblings_[i]->label_]++;
  }
  if (counts[partial] != partial_children_.size())
    return false;
  // Depending on how many partials there are, most full children will get
  // counted twice.
  if (counts[full] != (full_children_.size() * 2) - (2 - counts[partial]))
    return false;
  return true;
}

template <typename T>
void PQNode<T>::MoveFullChildren(PQNode<T>* new_node) {
  for (typename set<PQNode<T>*>::iterator i = full_children_.begin();
       i != full_children_.end(); ++i) {
    circular_link_.remove(*i);
    new_node->circular_link_.push_back(*i);
    (*i)->parent_ = new_node;
  }
}

template <typename T>
void PQNode<T>::ReplaceCircularLink(PQNode<T>* old_child, PQNode<T>* new_child) {
  circular_link_.remove(old_child);
  circular_link_.push_back(new_child);
}

// FindLeaves, FindFrontier, Reset, and Print are very similar recursive
// functions.  Each is a depth-first walk of the entire tree looking for data
// at the leaves.  
// TODO: Could probably be implemented better using function pointers.
template <typename T>
void PQNode<T>::FindLeaves(map<T, PQNode<T>*> &leafAddress) {
  if (type_ == leaf) {
    leafAddress[leaf_value_] = this;
  } else if (type_ == pnode) {
    // Recurse by asking each child in circular_link_ to find it's leaves.
    for (typename list<PQNode<T>*>::iterator i = circular_link_.begin();
         i!=circular_link_.end(); i++)
      (*i)->FindLeaves(leafAddress);
  } else if (type_ == qnode) {
    // Recurse by asking each child in my child list to find it's leaves.
    PQNode<T> *last    = NULL;
    PQNode<T> *current = endmost_children_[0];
    while (current) {
      current->FindLeaves(leafAddress);
      PQNode<T> *next = current->QNextChild(last);
      last    = current;
      current = next;
    }
  }
}
  
template <typename T>
void PQNode<T>::FindFrontier(list<T> &ordering) {
  if (type_ == leaf) {
    ordering.push_back(leaf_value_);
  } else if (type_ == pnode) {
    for(typename list<PQNode<T>*>::iterator i = circular_link_.begin();
        i != circular_link_.end();i++)
      (*i)->FindFrontier(ordering);
  } else if (type_ == qnode) {
    PQNode<T> *last    = NULL;
    PQNode<T> *current = endmost_children_[0];
    while (current) {
      current->FindFrontier(ordering);
      PQNode<T> *next = current->QNextChild(last);
      last    = current;
      current = next;
    }
  }
}

// Resets a bunch of temporary variables after the reduce walks
template <typename T>
void PQNode<T>::Reset() {
  if (type_ == pnode) {
    for (typename list<PQNode<T>*>::iterator i = circular_link_.begin();
         i != circular_link_.end(); i++)
      (*i)->Reset();
  } else if(type_ == qnode) {
    PQNode<T> *last    = NULL;
    PQNode<T> *current = endmost_children_[0];
    while (current) {
      current->Reset();
      PQNode<T> *next = current->QNextChild(last);
      last    = current;
      current = next;
    }
  }
  
  full_children_.clear();
  partial_children_.clear();
  label_                 = empty;
  mark_                  = unmarked;
  pertinent_child_count = 0;
  pertinent_leaf_count  = 0;
  pseudochild_           = false;
  pseudonode_            = false;
}

// Walks the tree from the top and prints the tree structure to the string out.
// Used primarily for debugging purposes.
template <typename T>
void PQNode<T>::Print(string *out) {
  if (type_ == leaf) {
    char value_str[10];
    sprintf(value_str, "%d", leaf_value_);
    *out += value_str;
  } else if (type_ == pnode) {
    *out += "(";
    for (typename list<PQNode<T>*>::iterator i = circular_link_.begin();
         i != circular_link_.end(); i++) {
      (*i)->Print(out);
      // Add a space if there are more elements remaining.
      if (++i != circular_link_.end())
        *out += " ";
      --i;
    }
    *out += ")";
  } else if (type_ == qnode) {
    *out += "[";
    PQNode<T> *last     = NULL;
    PQNode<T> *current  = endmost_children_[0];
    while (current) {
      current->Print(out);
      PQNode<T> *next = current->QNextChild(last);
      last     = current;
      current  = next;
      if (current)
        *out += " ";
    }
    *out += "]";
  }
}

/***** QNodeChildrenIterator<T> class *****/

template <typename T>
QNodeChildrenIterator<T>::QNodeChildrenIterator(PQNode<T>* parent,
					     PQNode<T>* begin_side) {
  parent_ = parent;
  Reset(begin_side);
}

template <typename T>
void QNodeChildrenIterator<T>::Reset(PQNode<T>* begin_side) {
  current_ = parent_->endmost_children_[0];
  if (begin_side)
    current_ = begin_side;
  prev_ = NULL;
  next_ = current_->immediate_siblings_[0];
}

template <typename T>
PQNode<T>* QNodeChildrenIterator<T>::Current() {
  return current_;
}

template <typename T>
void QNodeChildrenIterator<T>::NextPseudoNodeSibling() {
  // This should only be called from the first Next() call after Reset() and
  // only if the first subnode has two immediate siblings.  We want to advance
  // our iterator to the non-empty sibling of |current_|
  prev_ = current_;
  current_ = current_->ImmediateSiblingWithLabel(PQNode<T>::full);
  if (!current_)
    current_ = current_->ImmediateSiblingWithLabel(PQNode<T>::partial);
}

template <typename T>
void QNodeChildrenIterator<T>::Next() {
  // If the first child has 2 immediate siblings, then we are on
  // the edge of a pseudonode.
  if (IsDone())
    return;
  if (prev_ == NULL && current_->ImmediateSiblingCount() == 2) {
    NextPseudoNodeSibling();
  } else {
    prev_ = current_;
    current_ = next_;
  }

  if (current_) {
    next_ = current_->immediate_siblings_[0];
    if (next_ == prev_)
      next_ = current_->immediate_siblings_[1];
  }
}

template <typename T>
bool QNodeChildrenIterator<T>::IsDone() {
  return current_ == NULL;
}

#endif
