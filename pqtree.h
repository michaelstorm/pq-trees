/*
PQ-Tree class based on the paper "Testing for the Consecutive Onces Property,
Interval Graphs, and Graph Planarity Using PQ-Tree Algorithms" by Kellog S. 
Booth and George S. Lueker in the Journal of Computer and System Sciences 13,
335-379 (1976)
*/

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


#include <assert.h>
#include <list>
#include <set>
#include <vector>
#include <queue>
#include <map>
#include <iostream>
#include "pqnode.h"
#include "set_methods.h"

using namespace std;

#ifndef PQTREE_H
#define PQTREE_H

template <typename T>
class PQTree {
  private:
  
  // Root node of the PQTree
  PQNode<T> *root_;
  
  // The number of blocks of blocked nodes during the 1st pass
  int block_count_;

  // The number of blocked nodes during the 1st pass
  int blocked_nodes_;

  // A variable (0 or 1) which is a count of the number of virtual nodes which
  // are imagined to be in the queue during the bubbling up.
  int off_the_top_;

  // Keeps track of all reductions performed on this tree in order
  list<set<T> > reductions_;

  // Keeps a pointer to the leaf containing a particular value this map actually
  // increases the time complexity of the algorithm.  To fix, you can create an
  // array of items so that each item hashes to its leaf address in constant
  // time, but this is a tradeoff to conserve space.
  map<T, PQNode<T>*> leaf_address_;

  // A reference to a pseudonode that cannot be reached through the root
  // of the tree.  The pseudonode is a temporary node designed to handle
  // a special case in the first bubbling up pass it only exists during the
  // scope of the reduce operation
  PQNode<T>* pseudonode_;

  // true if a non-safe reduce has failed, tree is useless.
  bool invalid_;

  // Loops through the consecutive blocked siblings of an unblocked node
  // recursively unblocking the siblings.
  // Args:
  //   sibling: next sibling node to unblock
  //   parent: Node to set as the parent of all unblocked siblings
  //   last: the last node that was unblocked (used in the recursion)
  int UnblockSiblings(PQNode<T>* sibling, PQNode<T>* parent, PQNode<T>* last=NULL);

  // All of the templates for matching a reduce are below.  The template has a
  // letter describing which type of node it refers to and a number indicating
  // the index of the template for that letter.  These are the same indices in
  // the Booth & Lueker paper.  The return value indicates whether or not the
  // pattern accurately matches the template
  bool TemplateL1(PQNode<T>* candidate_node);
  bool TemplateQ1(PQNode<T>* candidate_node);
  bool TemplateQ2(PQNode<T>* candidate_node);
  bool TemplateQ3(PQNode<T>* candidate_node);
  bool TemplateP1(PQNode<T>* candidate_node, bool is_reduction_root);
  bool TemplateP2(PQNode<T>* candidate_node);
  bool TemplateP3(PQNode<T>* candidate_node);
  bool TemplateP4(PQNode<T>* candidate_node);
  bool TemplateP5(PQNode<T>* candidate_node);
  bool TemplateP6(PQNode<T>* candidate_node);
  
  // This procedure is the first pass of the Booth&Leuker PQTree algorithm
  // It processes the pertinent subtree of the PQ-Tree to determine the mark
  // of every node in that subtree
  // the pseudonode, if created, is returned so that it can be deleted at
  // the end of the reduce step
  bool Bubble(set<T> S);

  bool ReduceStep(set<T> S);
        
 public:  
  // Default constructor - constructs a tree using a set
  // Only reductions using elements of that set will succeed
  PQTree(set<T> S);
  PQTree(const PQTree& to_copy);
  ~PQTree();
 
  // Returns the root PQNode<T> used for exploring the tree.
  PQNode<T>* Root();
 
  // Mostly for debugging purposes, Prints the tree to standard out
  string Print();
  
  // Cleans up pointer mess caused by having a pseudonode
  void CleanPseudo();

  // Reduces the tree but protects if from becoming invalid if the reduction
  // fails, takes more time.
  bool SafeReduce(set<T>);
  bool SafeReduceAll(list<set<T> >);

  //reduces the tree - tree can become invalid, making all further
  //reductions fail
  bool Reduce(set<T> S);
  bool ReduceAll(list<set<T> > L);
  
  // Returns 1 possible frontier, or ordering preserving the reductions
  list<T> Frontier();
  
  // Assignment operator
  PQTree& operator=(const PQTree& to_copy);  

  // Copies to_copy.  Used for the copy constructor and assignment operator.
  void CopyFrom(const PQTree& to_copy);
  
  
  // Returns a frontier not including leaves that were not a part of any
  // reduction
  list<T> ReducedFrontier();
  
  // Returns the reductions that have been performed on this tree.
  list<set<T> > GetReductions();
  
  // Returns the set of all elements on which a reduction was performed.
  set<T> GetContained();
};

template <typename T>
PQTree<T>::PQTree(const PQTree<T>& to_copy) {
  CopyFrom(to_copy);
}

template <typename T>
PQTree<T>& PQTree<T>::operator=(const PQTree<T>& to_copy) {
  if (&to_copy != this)
    CopyFrom(to_copy);
  return *this;
}

template <typename T>
void PQTree<T>::CopyFrom(const PQTree<T>& to_copy) {
  root_                 = new PQNode<T>(*to_copy.root_);
  block_count_   = to_copy.block_count_;
  blocked_nodes_ = to_copy.blocked_nodes_;
  invalid_         = to_copy.invalid_;
  off_the_top_   = to_copy.off_the_top_;
  pseudonode_         = NULL;
  reductions_         = to_copy.reductions_;

  leaf_address_.clear();
  root_->FindLeaves(leaf_address_);
}

template <typename T>
int PQTree<T>::UnblockSiblings(PQNode<T>* candidate_node,
                            PQNode<T>* parent,
                            PQNode<T>* last) {
  int unblocked_count = 0;
  if (candidate_node->mark_ == PQNode<T>::blocked) {
    // Unblock the current sibling and set it's parent to |parent|.
    unblocked_count++;
    candidate_node->mark_ = PQNode<T>::unblocked;
    candidate_node->parent_ = parent;

    // Unblock adjacent siblings recursively
    for (int i = 0; i < 2 && candidate_node->immediate_siblings_[i]; ++i) {
      PQNode<T> *sibling = candidate_node->immediate_siblings_[i];
      if (sibling != last)
        unblocked_count += UnblockSiblings(sibling, parent, candidate_node);
    }
  }
  return unblocked_count;
}

// All the different templates for matching a reduce are below.  The template
// has a letter describing which type of node it refers to and a number
// indicating the index of the template for that letter.  These are the same
// indices found in the Booth & Lueker paper.
//
// Each template method is called in order within the reduce phase.  The
// template method determines if the Node passed to it is of the form that the
// template method can reduce.  If so, the method processes the Node and returns
// true.  If not, the method makes no changes and returns false.  Methods
// should be attempted in order until one of them returns true and then the
// remainder of the methods should be skipped.  The template ordering is:
// L1, P1, P2, P3, P4, P5, P6, Q1, Q2, Q3

template <typename T>
bool PQTree<T>::TemplateL1(PQNode<T>* candidate_node) {
  // L1's pattern is simple: the node is a leaf node.
  if (candidate_node->type_ != PQNode<T>::leaf)
    return false;

  candidate_node->LabelAsFull();
  return true;
}

template <typename T>
bool PQTree<T>::TemplateQ1(PQNode<T>* candidate_node) {
  // Q1's Pattern is a Q-Node that has only full children.
  if (candidate_node->type_ != PQNode<T>::qnode)
    return false;
  for (QNodeChildrenIterator<T> it(candidate_node); !it.IsDone(); it.Next()) {
    if (it.Current()->label_ != PQNode<T>::full)
      return false;
  }

  candidate_node->LabelAsFull();
  return true;
}

template <typename T>
bool PQTree<T>::TemplateQ2(PQNode<T>* candidate_node) {
  // Q2's pattern is a Q-Node that either:
  // 1) contains consecutive full children with one end of the consecutive 
  //    ordering being one of |candidate_node|'s |endmost_children|, also full.
  // 2) contains a single partial child that is one of |candidate_node|'s
  //    |endmost_children|.
  // 3) One of |candidate_node|'s |endmost_childdren| is full, consecutively
  //    followed by 0 or more full children, followed by one partial child.
  if (candidate_node->type_ != PQNode<T>::qnode ||
      candidate_node->pseudonode_ ||
      candidate_node->partial_children_.size() > 1 ||
      !candidate_node->ConsecutiveFullPartialChildren())
    return false;

  bool has_partial = candidate_node->partial_children_.size() > 0;
  bool has_full = candidate_node->full_children_.size() > 0;
  
  if (has_full && !candidate_node->EndmostChildWithLabel(PQNode<T>::full))
    return false;
  if (!has_full && !candidate_node->EndmostChildWithLabel(PQNode<T>::partial))
    return false;
  
  // If there is a partial child, merge it's children into the candidate_node.
  if (has_partial) {
    PQNode<T>* to_merge = *candidate_node->partial_children_.begin();
    for (int i = 0; i < 2; ++i) {
      PQNode<T>* child = to_merge->endmost_children_[i];
      PQNode<T>* sibling = to_merge->ImmediateSiblingWithLabel(child->label_);
      if (sibling) {
        sibling->ReplaceImmediateSibling(to_merge, child);
      } else {
        candidate_node->ReplaceEndmostChild(to_merge, child);
        child->parent_ = candidate_node;
      }
    }
    to_merge->ForgetChildren();
    delete to_merge;
  }

  candidate_node->label_ = PQNode<T>::partial;
  if (candidate_node->parent_)
    candidate_node->parent_->partial_children_.insert(candidate_node);
  return true;
}

template <typename T>
bool PQTree<T>::TemplateQ3(PQNode<T>* candidate_node) {
  // Q3's pattern is a Q-Node that contains 0-2 partial children.  It can
  // contain any number of empty and full children, but any full children must
  // be consecutive and sandwiched between any partial children.  Unlike Q2, 
  // the consecutive full and partial children need not be endmost children.
  if (candidate_node->type_ != PQNode<T>::qnode ||
      candidate_node->partial_children_.size() > 2 ||
      !candidate_node->ConsecutiveFullPartialChildren())
    return false;

  // Merge each of the partial children into |candidate_node|'s children
  for (typename set<PQNode<T>*>::iterator j = candidate_node->partial_children_.begin();
       j != candidate_node->partial_children_.end(); j++) {
    PQNode<T>* to_merge = *j;
    for (int i = 0; i < 2; ++i) {
      PQNode<T>* sibling = to_merge->immediate_siblings_[i];
      if (sibling) {
	PQNode<T>* child = to_merge->EndmostChildWithLabel(sibling->label_);
	if (!child)
	  child = to_merge->EndmostChildWithLabel(PQNode<T>::full);
        sibling->ReplaceImmediateSibling(to_merge, child);
      } else {
	PQNode<T>* empty_child = to_merge->EndmostChildWithLabel(PQNode<T>::empty);
        empty_child->parent_ = candidate_node;
        candidate_node->ReplaceEndmostChild(to_merge, empty_child);
      }
    }

    to_merge->ForgetChildren();
    delete to_merge;
  }
  return true;
}

// A note here.  An error in the Booth and Leuker Algorithm fails to consider
// the case where a P-node is full, is the pertinent root, and is not an endmost
// child of a q-node.  In this case, we need to know that the P-node is a
// pertinent root and not try to update its parent whose pointer is possibly
// invalid.
template <typename T>
bool PQTree<T>::TemplateP1(PQNode<T>* candidate_node, bool is_reduction_root) {
    // P1's pattern is a P-Node with all full children.
    if (candidate_node->type_ != PQNode<T>::pnode ||
        candidate_node->full_children_.size() != candidate_node->ChildCount())
      return false;

    candidate_node->label_ = PQNode<T>::full;
    if (!is_reduction_root)
      candidate_node->parent_->full_children_.insert(candidate_node);
    return true;
}

template <typename T>
bool PQTree<T>::TemplateP2(PQNode<T>* candidate_node) {
  // P2's pattern is a P-Node at the root of the perinent subtree containing
  // both empty and full children.
  if (candidate_node->type_ != PQNode<T>::pnode ||
      !candidate_node->partial_children_.empty())
    return false;

  // Move candidate_node's full children into their own P-node
  if (candidate_node->full_children_.size() >= 2) {
    PQNode<T>* new_pnode = new PQNode<T>;
    new_pnode->type_ = PQNode<T>::pnode;
    new_pnode->parent_ = candidate_node;
    candidate_node->MoveFullChildren(new_pnode);
    candidate_node->circular_link_.push_back(new_pnode);
  }
  // Mark the root partial
  candidate_node->label_ = PQNode<T>::partial;

  return true;
}
template <typename T>
bool PQTree<T>::TemplateP3(PQNode<T>* candidate_node) {
  // P3's pattern is a P-Node not at the root of the perinent subtree
  // containing both empty and full children.
  if (candidate_node->type_ != PQNode<T>::pnode ||
      !candidate_node->partial_children_.empty())
    return false;

  // P3's replacement is to create a Q-node that places all of the full
  // elements in a single P-Node child and all of the empty elements in a
  // single Q-Node child.  This new Q-Node is called a pseudonode as it isn't
  // properly formed (Q-Nodes should have at least 3 children) and will not
  // survive in it's current form to the end of the reduction.
  PQNode<T>* new_qnode = new PQNode<T>;
  new_qnode->type_ = PQNode<T>::qnode;
  new_qnode->label_ = PQNode<T>::partial;
  candidate_node->parent_->ReplacePartialChild(candidate_node, new_qnode);

  // Set up a |full_child| of |new_qnode| containing all of |candidate_node|'s
  // full children.
  PQNode<T>* full_child;
  if (candidate_node->full_children_.size() == 1) {
    full_child = *candidate_node->full_children_.begin();
    candidate_node->circular_link_.remove(full_child);
  } else {
    full_child = new PQNode<T>;
    full_child->type_ = PQNode<T>::pnode;
    full_child->label_ = PQNode<T>::full;
    candidate_node->MoveFullChildren(full_child);
  }
  full_child->parent_ = new_qnode;
  full_child->label_ = PQNode<T>::full;
  new_qnode->endmost_children_[0] = full_child;
  new_qnode->full_children_.insert(full_child);

  // Set up a |empty_child| of |new_qnode| containing all of |candidate_node|'s
  // empty children.
  PQNode<T>* empty_child;
  if (candidate_node->circular_link_.size() == 1) {
    empty_child = *candidate_node->circular_link_.begin();
    candidate_node->circular_link_.clear();
    delete candidate_node;
  } else {
    empty_child = candidate_node;
  }
  empty_child->parent_ = new_qnode;
  empty_child->label_ = PQNode<T>::empty;
  new_qnode->endmost_children_[1] = empty_child;

  // Update the immediate siblings links (erasing the old ones if present)
  empty_child->immediate_siblings_[0] = full_child;
  full_child->immediate_siblings_[0] = empty_child;

  return true;
}

template <typename T>
bool PQTree<T>::TemplateP4(PQNode<T>* candidate_node) {
  // P4's pattern is a P-Node at the root of the perinent subtree containing
  // one partial child and any number of empty/full children.
  if (candidate_node->type_ != PQNode<T>::pnode ||
      candidate_node->partial_children_.size() != 1)
    return false;

  PQNode<T>* partial_qnode = *candidate_node->partial_children_.begin();
  assert(partial_qnode->type_ == PQNode<T>::qnode);
  PQNode<T>* empty_child = partial_qnode->EndmostChildWithLabel(PQNode<T>::empty);
  PQNode<T>* full_child = partial_qnode->EndmostChildWithLabel(PQNode<T>::full);

  if (!empty_child || !full_child)
    return false;

  // Move the full children of |candidate_node| to children of |partial_qnode|.
  // TODO: Lots of redundancy between lines 271-287, 326-342, and 345-359.
  if (!candidate_node->full_children_.empty()) {
    PQNode<T> *full_children_root;
    if (candidate_node->full_children_.size() == 1) {
      full_children_root = *candidate_node->full_children_.begin();
      candidate_node->circular_link_.remove(full_children_root);
    } else {
      full_children_root = new PQNode<T>;
      full_children_root->type_ = PQNode<T>::pnode;
      full_children_root->label_ = PQNode<T>::full;
      candidate_node->MoveFullChildren(full_children_root);
    }
    full_children_root->parent_ = partial_qnode;
    partial_qnode->ReplaceEndmostChild(full_child, full_children_root);
    partial_qnode->full_children_.insert(full_children_root);
    full_child->AddImmediateSibling(full_children_root);
    full_children_root->AddImmediateSibling(full_child);
  }

  // If |candidate_node| now only has one child, get rid of |candidate_node|.
  if (candidate_node->circular_link_.size() == 1) {
    if (candidate_node->parent_) {
      candidate_node->parent_->ReplaceChild(candidate_node, partial_qnode);
    } else {
      root_ = partial_qnode;
      partial_qnode->parent_ = NULL;
    }
    candidate_node->circular_link_.clear();
    delete candidate_node;
  }
  return true;
}

template <typename T>
bool PQTree<T>::TemplateP5(PQNode<T>* candidate_node) {
  // P4's pattern is a P-Node not at the root of the perinent subtree
  // containing one partial child and any number of empty/full children.
  if (candidate_node->type_ != PQNode<T>::pnode ||
      candidate_node->partial_children_.size() != 1)
    return false;

  // |partial_qnode| will become the pertinent subtree root after replacement.
  PQNode<T>* partial_qnode = *candidate_node->partial_children_.begin();
  assert(partial_qnode->type_ == PQNode<T>::qnode);
  PQNode<T>* empty_child = partial_qnode->EndmostChildWithLabel(PQNode<T>::empty);
  PQNode<T>* full_child = partial_qnode->EndmostChildWithLabel(PQNode<T>::full);
  PQNode<T>* empty_sibling = candidate_node->CircularChildWithLabel(PQNode<T>::empty);

  if (!empty_child || !full_child)
    return false;

  // Move partial_qnode from candidate_node's child to it's parent's child
  candidate_node->parent_->ReplaceChild(candidate_node, partial_qnode);
  partial_qnode->pertinent_leaf_count = candidate_node->pertinent_leaf_count;
  candidate_node->circular_link_.remove(partial_qnode);

  // Move the full children of |candidate_node| to children of |partial_qnode|.
  if (!candidate_node->full_children_.empty()) {
    PQNode<T> *full_children_root;
    if (candidate_node->full_children_.size() == 1) {
      full_children_root = *candidate_node->full_children_.begin();
      candidate_node->circular_link_.remove(full_children_root);
    } else {
      full_children_root = new PQNode<T>;
      full_children_root->type_ = PQNode<T>::pnode;
      full_children_root->label_ = PQNode<T>::full;
      candidate_node->MoveFullChildren(full_children_root);
    }

    full_children_root->parent_ = partial_qnode;
    full_child->AddImmediateSibling(full_children_root);
    full_children_root->AddImmediateSibling(full_child);
    partial_qnode->ReplaceEndmostChild(full_child, full_children_root);
  }

  // If candidate_node still has some empty children, insert them
  if (candidate_node->ChildCount()) {
    PQNode<T> *empty_children_root;
    if (candidate_node->ChildCount() == 1) {
      empty_children_root = empty_sibling;
    } else {
      empty_children_root = candidate_node;
      empty_children_root->label_ = PQNode<T>::empty;
      empty_children_root->ClearImmediateSiblings();
    }

    empty_children_root->parent_ = partial_qnode;
    empty_child->AddImmediateSibling(empty_children_root);
    empty_children_root->AddImmediateSibling(empty_child);
    partial_qnode->ReplaceEndmostChild(empty_child, empty_children_root);
  }
  if (candidate_node->ChildCount() < 2) {
    // We want to delete candidate_node, but not it's children.
    candidate_node->circular_link_.clear();
    delete candidate_node;
  }

  return true;
}

template <typename T>
bool PQTree<T>::TemplateP6(PQNode<T>* candidate_node)
{
  if (candidate_node->type_ != PQNode<T>::pnode ||
      candidate_node->partial_children_.size() != 2)
    return false;

  // TODO: Convert these to an array so we don't have 2 of everything.
  PQNode<T>* partial_qnode1 = *candidate_node->partial_children_.begin();
  PQNode<T>* partial_qnode2 = *(++(candidate_node->partial_children_.begin()));

  PQNode<T>* empty_child1 = partial_qnode1->EndmostChildWithLabel(PQNode<T>::empty);
  PQNode<T>* full_child1 = partial_qnode1->EndmostChildWithLabel(PQNode<T>::full);
  if (!empty_child1 || !full_child1)
    return false;

  PQNode<T>* empty_child2 = partial_qnode2->EndmostChildWithLabel(PQNode<T>::empty);
  PQNode<T>* full_child2 = partial_qnode2->EndmostChildWithLabel(PQNode<T>::full);
  if (!empty_child2 || !full_child2)
    return false;

  // Move the full children of candidate_node to be children of partial_qnode1
  if (!candidate_node->full_children_.empty()) {
    PQNode<T> *full_children_root = NULL;
    if (candidate_node->full_children_.size() == 1) {
      full_children_root = *candidate_node->full_children_.begin();
      candidate_node->circular_link_.remove(full_children_root);
    } else {
      // create full_children_root to be a new p-node
      full_children_root = new PQNode<T>;
      full_children_root->type_ = PQNode<T>::pnode;
      full_children_root->label_ = PQNode<T>::full;
      candidate_node->MoveFullChildren(full_children_root); 
    }
    full_children_root->parent_ = partial_qnode1;
    full_child2->parent_ = partial_qnode1;

    full_child1->AddImmediateSibling(full_children_root);
    full_child2->AddImmediateSibling(full_children_root);
    full_children_root->AddImmediateSibling(full_child1);
    full_children_root->AddImmediateSibling(full_child2);
  } else  {
    full_child1->AddImmediateSibling(full_child2);
    full_child2->AddImmediateSibling(full_child1);
  }
  partial_qnode1->ReplaceEndmostChild(full_child1, empty_child2);
  empty_child2->parent_ = partial_qnode1;

  // We dont need |partial_qnode2| any more
  candidate_node->circular_link_.remove(partial_qnode2);
  partial_qnode2->ForgetChildren();
  delete partial_qnode2;

  // If |candidate_node| now only has one child, get rid of |candidate_node|.
  if (candidate_node->circular_link_.size() == 1) {
    partial_qnode1->parent_ = candidate_node->parent_;
    partial_qnode1->pertinent_leaf_count = candidate_node->pertinent_leaf_count;
    partial_qnode1->label_ = PQNode<T>::partial;

    if (candidate_node->parent_) {
      candidate_node->parent_->partial_children_.insert(partial_qnode1);
      if (candidate_node->parent_->type_ == PQNode<T>::pnode) {
        candidate_node->parent_->ReplaceCircularLink(candidate_node,
                                                     partial_qnode1);
      } else {
        for (int i = 0; i < 2 && candidate_node->immediate_siblings_[i]; ++i) {
	  	  PQNode<T>* sibling = candidate_node->immediate_siblings_[i];
          sibling->ReplaceImmediateSibling(candidate_node, partial_qnode1);
        }
        candidate_node->parent_->ReplaceEndmostChild(candidate_node,
                                                     partial_qnode1);
      }
    }
    else {
      root_ = partial_qnode1;
      partial_qnode1->parent_ = NULL;

      // Delete candidate_node, but not it's children.
      candidate_node->circular_link_.clear();
      delete candidate_node;
    }
  }
  return true;
}

// This procedure is the first pass of the Booth & Leuker PQTree<T> algorithm.
// It processes the pertinent subtree of the PQ-Tree to determine the mark
// of every node in that subtree.  The pseudonode, if created, is returned so
// that it can be deleted at the end of the reduce step.
template <typename T>
bool PQTree<T>::Bubble(set<T> reduction_set) {
  queue<PQNode<T>*> q;
  block_count_ = 0;
  blocked_nodes_ = 0;
  off_the_top_ = 0;

  // Stores nodes that aren't unblocked
  set<PQNode<T>*> blocked_list;

  // Insert the set's leaves into the queue
  for (typename set<T>::iterator i = reduction_set.begin();
       i != reduction_set.end(); ++i) {
    PQNode<T> *temp = leaf_address_[*i];
    assert (temp);
    q.push(temp);
  }

  while (q.size() + block_count_ + off_the_top_ > 1) {
    if (q.empty())
      return false;

    PQNode<T>* candidate_node = q.front();
    q.pop();
    candidate_node->mark_ = PQNode<T>::blocked;

    // Get the set of blocked and PQNode<T>::unblocked siblings
    set<PQNode<T>*> unblocked_siblings;
    set<PQNode<T>*> blocked_siblings;
    for (int i = 0; i < 2 && candidate_node->immediate_siblings_[i]; ++i) {
      PQNode<T>* sibling = candidate_node->immediate_siblings_[i];
      if (sibling->mark_ == PQNode<T>::blocked) {
        blocked_siblings.insert(sibling);
      } else if (sibling->mark_ == PQNode<T>::unblocked) {
        unblocked_siblings.insert(sibling);
      }
    }

    // We can unblock |candidate_node| if any of there conditions is met:
    //  - 1 or more of its immediate siblings is unblocked.
    //  - It has 1 immediate sibling meaning it is a corner child of a q node.
    //  - It has 0 immediate siblings meaning it is a p node.
    if (!unblocked_siblings.empty()) {
      candidate_node->parent_ = (*unblocked_siblings.begin())->parent_;
      candidate_node->mark_ = PQNode<T>::unblocked;
    } else if (candidate_node->ImmediateSiblingCount() < 2) {
      candidate_node->mark_ = PQNode<T>::unblocked;
    }

    // If |candidate_node| is unblocked, we can process it.
    if (candidate_node->mark_ == PQNode<T>::unblocked) {
      int list_size = 0;
      if (!blocked_siblings.empty()) {
        candidate_node->mark_ = PQNode<T>::blocked;
        list_size = UnblockSiblings(candidate_node, candidate_node->parent_);
        candidate_node->parent_->pertinent_child_count += list_size - 1;
      }

      if (!candidate_node->parent_) {
        off_the_top_ = 1;
      } else {
        candidate_node->parent_->pertinent_child_count++;
        if (candidate_node->parent_->mark_ == PQNode<T>::unmarked) {
          q.push(candidate_node->parent_);
          candidate_node->parent_->mark_ = PQNode<T>::queued;
        }
      }
      block_count_ -= blocked_siblings.size();
      blocked_nodes_ -= list_size; 
    } else {
      block_count_ += 1 - blocked_siblings.size();
      blocked_nodes_ += 1;
      blocked_list.insert(candidate_node);
    }
  }

  if (block_count_ > 1 || (off_the_top_ == 1 && block_count_ != 0))
    return false;

  int correct_blocked_count = 0;
  for (typename set<PQNode<T>*>::iterator blocked_node = blocked_list.begin();
      blocked_node != blocked_list.end(); ++blocked_node) {
      if ((*blocked_node)->mark_ == PQNode<T>::blocked)
        correct_blocked_count++;
  }

  // In this case, we have a block that is contained within a Q-node.  We must
  // assign a psuedonode to handle it.
  if (block_count_ == 1 && correct_blocked_count > 1) {
    pseudonode_ = new PQNode<T>;
    pseudonode_->type_ = PQNode<T>::qnode;
    pseudonode_->pseudonode_ = true;
    pseudonode_->pertinent_child_count = 0;

    // Find the blocked nodes and which of those are endmost children.
    int side = 0;
    for (typename set<PQNode<T>*>::iterator i = blocked_list.begin();
        i != blocked_list.end(); ++i) {
      PQNode<T>* blocked = *i;
      if (blocked->mark_ == PQNode<T>::blocked) {
        pseudonode_->pertinent_child_count++;
        pseudonode_->pertinent_leaf_count += blocked->pertinent_leaf_count;
        int count = 0;  // count number of immediate siblings.
        for (int j = 0; j < 2 && blocked->immediate_siblings_[j]; ++j) {
          PQNode<T>* sibling = blocked->immediate_siblings_[j];
          if (sibling->mark_ == PQNode<T>::blocked) {
            ++count;
          } else {
            blocked->RemoveImmediateSibling(sibling);
            sibling->RemoveImmediateSibling(blocked);
            pseudonode_->pseudo_neighbors_[side] = sibling;
          }
        }
        blocked->parent_ = pseudonode_;
        blocked->pseudochild_ = true;
        if (count < 2)
          pseudonode_->endmost_children_[side++] = blocked;
      }
    }
  }
  return true;
}

template <typename T>
bool PQTree<T>::ReduceStep(set<T> reduction_set) {
  // Build a queue with all the pertinent leaves in it
  queue<PQNode<T>*> q;
  for (typename set<T>::iterator i = reduction_set.begin();
       i != reduction_set.end(); i++) {
    PQNode<T>* candidate_node = leaf_address_[*i];
    if (candidate_node == NULL)
      return false;
    candidate_node->pertinent_leaf_count = 1;
    q.push(candidate_node);
  }

  while (!q.empty()) {
    // Remove candidate_node from the front of the queue
    PQNode<T>* candidate_node = q.front();
    q.pop();

    // We test against different templates depending on whether |candidate_node|
    // is the root of the pertinent subtree.
    if (candidate_node->pertinent_leaf_count < reduction_set.size()) {
      PQNode<T>* candidate_parent = candidate_node->parent_;
      candidate_parent->pertinent_leaf_count += 
	  candidate_node->pertinent_leaf_count;
      candidate_parent->pertinent_child_count--;
      // Push |candidate_parent| onto the queue if it no longer has any
      // pertinent children.
      if (candidate_parent->pertinent_child_count == 0)
        q.push(candidate_parent);

      // Test against each template in turn until one of them returns true.
      if      (TemplateL1(candidate_node)) {}
      else if (TemplateP1(candidate_node, /*is_reduction_root=*/ false)) {}
      else if (TemplateP3(candidate_node)) {}
      else if (TemplateP5(candidate_node)) {}
      else if (TemplateQ1(candidate_node)) {}
      else if (TemplateQ2(candidate_node)) {}
      else {
        CleanPseudo();
        return false;
      }
    } else {  // candidate_node is the root of the reduction subtree
      // Test against each template in turn until one of them returns true.
      if      (TemplateL1(candidate_node)) {}
      else if (TemplateP1(candidate_node, /*is_reduction_root=*/ true)) {}
      else if (TemplateP2(candidate_node)) {}
      else if (TemplateP4(candidate_node)) {}
      else if (TemplateP6(candidate_node)) {}
      else if (TemplateQ1(candidate_node)) {}
      else if (TemplateQ2(candidate_node)) {}
      else if (TemplateQ3(candidate_node)) {}
      else {
        CleanPseudo();
        return false;
      }
    }
  }
  CleanPseudo();
  return true;
}

template <typename T>
void PQTree<T>::CleanPseudo() {
  if (pseudonode_) {
    for (int i = 0; i < 2; i++) {
      pseudonode_->endmost_children_[i]->AddImmediateSibling(
          pseudonode_->pseudo_neighbors_[i]);
      pseudonode_->pseudo_neighbors_[i]->AddImmediateSibling(
          pseudonode_->endmost_children_[i]);
    }

    pseudonode_->ForgetChildren();
    delete pseudonode_;
    pseudonode_ = NULL;
  }
}

// Basic constructor from an initial set.
template <typename T>
PQTree<T>::PQTree(set<T> reduction_set) {
  // Set up the root node as a P-Node initially.
  root_ = new PQNode<T>;
  root_->type_ = PQNode<T>::pnode;
  invalid_ = false;
  pseudonode_ = NULL;
  block_count_ = 0;
  blocked_nodes_ = 0;
  off_the_top_ = 0;
  for (typename set<T>::iterator i = reduction_set.begin();
       i != reduction_set.end(); i++) {
    PQNode<T> *new_node;
    new_node = new PQNode<T>(*i);
    leaf_address_[*i] = new_node;
    new_node->parent_ = root_;
    new_node->type_ = PQNode<T>::leaf;
    root_->circular_link_.push_back(new_node);
  }
}

template <typename T>
PQNode<T>* PQTree<T>::Root() {
  return root_;
}

template <typename T>
string PQTree<T>::Print() {
  string out;
  root_->Print(&out);
  return out;
}

//reduces the tree but protects if from becoming invalid
//if the reduction fails, takes more time
template <typename T>
bool PQTree<T>::SafeReduce(set<T> S) {
  //using a backup copy to enforce safety
  PQTree<T> toCopy(*this);

  if (!Reduce(S)) {
    //reduce failed, so perform a copy
    root_ = new PQNode<T>(*toCopy.root_);
    block_count_ = toCopy.block_count_;
    blocked_nodes_ = toCopy.blocked_nodes_;
    off_the_top_ = toCopy.off_the_top_;
    invalid_ = toCopy.invalid_;
    leaf_address_.clear();
    root_->FindLeaves(leaf_address_);
    return false;
  }
  return true;
}

template <typename T>
bool PQTree<T>::SafeReduceAll(list<set<T> > L) {
  //using a backup copy to enforce safety
  PQTree<T> toCopy(*this);
  if (!ReduceAll(L)) {
    //reduce failed, so perform a copy
    root_ = new PQNode<T>(*toCopy.root_);
    block_count_ = toCopy.block_count_;
    blocked_nodes_ = toCopy.blocked_nodes_;
    off_the_top_ = toCopy.off_the_top_;
    invalid_ = toCopy.invalid_;
    leaf_address_.clear();
    root_->FindLeaves(leaf_address_);
    return false;
  }
  return true;
}


template <typename T>
bool PQTree<T>::Reduce(set<T> reduction_set) {
  if (reduction_set.size() < 2) {
    reductions_.push_back(reduction_set);
    return true;
  }
  if (invalid_)
    return false;
  if (!Bubble(reduction_set)) {
    invalid_ = true;
    return false;
  }
  if (!ReduceStep(reduction_set)) {
    invalid_ = true;
    return false;
  }

  if (pseudonode_) {
    pseudonode_->ForgetChildren();
    delete pseudonode_;
  }

  // Reset all the temporary variables for the next round.
  root_->Reset();

  // Store the reduction set for later lookup.
  reductions_.push_back(reduction_set);
  return true;
}

template <typename T>
bool PQTree<T>::ReduceAll(list<set<T> > L) {
  for (typename list<set<T> >::iterator S = L.begin(); S != L.end(); S++) {
    if (!Reduce(*S))
      return false;
  }
  return true;
}

template <typename T>
list<T> PQTree<T>::Frontier() {
  list<T> out;
  root_->FindFrontier(out);
  return out;
}

template <typename T>
list<T> PQTree<T>::ReducedFrontier() {
  list<T> out, inter;
  root_->FindFrontier(inter);
  set<T> allContained;
  for (typename list<set<T> >::iterator j = reductions_.begin();
      j != reductions_.end(); j++) {
      allContained = SetMethods::SetUnion(allContained, *j);
  }
  for (typename list<T>::iterator j = inter.begin(); j != inter.end(); j++) {
    if (SetMethods::SetFind(allContained, *j))
      out.push_back(*j);
  }
  return out;
}

template <typename T>
list<set<T> > PQTree<T>::GetReductions() {
  return reductions_;
}

template <typename T>
set<T> PQTree<T>::GetContained() {
  set<T> out;
  for (typename list<set<T> >::iterator i = reductions_.begin();
      i!=reductions_.end(); i++) {
    out = SetMethods::SetUnion(out,*i);
  }
  return out;
}

// Default destructor, Needs to delete the root.
template <typename T>
PQTree<T>::~PQTree<T>() {
  delete root_;
}

#endif
