/*
    Copyright (C) 2022-2023 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <t.flouris@ucl.ac.uk>,
    Department of Genetics, Evolution and Environment,
    University College London, Gower Street, London WC1E 6BT, England
*/

#include "bpp-tools.h"

static char * export_newick_recursive(const rnode_t * root,
                                      char * (*cb_serialize)(const rnode_t *))
{
  char * newick;
  long size_alloced;
  assert(root != NULL);

  if (!(root->left) || !(root->right))
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
      size_alloced = (long)strlen(newick);
    }
    else
    {
      size_alloced = xasprintf(&newick, "%s:%.3f", root->label, root->length);
    }
  }
  else
  {
    char * subtree1 = export_newick_recursive(root->left,cb_serialize);
    if (subtree1 == NULL)
    {
      return NULL;
    }
    char * subtree2 = export_newick_recursive(root->right,cb_serialize);
    if (subtree2 == NULL)
    {
      free(subtree1);
      return NULL;
    }

    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      size_alloced = xasprintf(&newick,
                              "(%s,%s)%s",
                              subtree1,
                              subtree2,
                              temp);
      free(temp);
    }
    else
    {
      size_alloced = xasprintf(&newick,
                              "(%s,%s)%s:%.3f",
                              subtree1,
                              subtree2,
                              root->label ? root->label : "",
                              root->length);
    }
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
    fatal("Memory allocation during newick export failed.");

  return newick;
}

char * rtree_export_newick(const rnode_t * root,
                           char * (*cb_serialize)(const rnode_t *))
{
  char * newick;
  long size_alloced;
  if (!root) return NULL;

  if (!(root->left) || !(root->right))
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
      size_alloced = (long)strlen(newick);
    }
    else
    {
      size_alloced = xasprintf(&newick, "%s:%.3f", root->label, root->length);
    }
  }
  else
  {
    char * subtree1 = export_newick_recursive(root->left,cb_serialize);
    if (subtree1 == NULL)
      fatal("Unable to allocate enough memory.");

    char * subtree2 = export_newick_recursive(root->right,cb_serialize);
    if (subtree2 == NULL)
      fatal("Unable to allocate enough memory.");

    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      size_alloced = xasprintf(&newick,
                              "(%s,%s)%s",
                              subtree1,
                              subtree2,
                              temp);
      free(temp);
    }
    else
    {
      size_alloced = xasprintf(&newick,
                              "(%s,%s)%s:%.3f;",
                              subtree1,
                              subtree2,
                              root->label ? root->label : "",
                              root->length);
    }
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0)
    fatal("Memory allocation during newick export failed");

  return newick;
}

static void rtree_traverse_postorder(rnode_t * node,
                                     int (*cbtrav)(rnode_t *),
                                     unsigned int * index,
                                     rnode_t ** outbuffer)
{
  if (!node->left)
  {
    if (cbtrav(node))
    {
      outbuffer[*index] = node;
      *index = *index + 1;
    }
    return;
  }
  if (!cbtrav(node))
    return;

  rtree_traverse_postorder(node->left, cbtrav, index, outbuffer);
  rtree_traverse_postorder(node->right, cbtrav, index, outbuffer);

  outbuffer[*index] = node;
  *index = *index + 1;
}

static void rtree_traverse_preorder(rnode_t * node,
                                    int (*cbtrav)(rnode_t *),
                                    unsigned int * index,
                                    rnode_t ** outbuffer)
{
  if (!node->left)
  {
    if (cbtrav(node))
    {
      outbuffer[*index] = node;
      *index = *index + 1;
    }
    return;
  }
  if (!cbtrav(node))
    return;

  outbuffer[*index] = node;
  *index = *index + 1;

  rtree_traverse_preorder(node->left, cbtrav, index, outbuffer);
  rtree_traverse_preorder(node->right, cbtrav, index, outbuffer);

}

int rtree_traverse(rnode_t * root,
                   int traversal,
                   int (*cbtrav)(rnode_t *),
                   rnode_t ** outbuffer,
                   unsigned int * trav_size)
{
  *trav_size = 0;
  if (!root->left) return BPP_FAILURE;

  /* we will traverse an unrooted tree in the following way

           root
            /\
           /  \
        left   right

     at each node the callback function is called to decide whether we
     are going to traversing the subtree rooted at the specific node */

  if (traversal == TREE_TRAVERSE_POSTORDER)
    rtree_traverse_postorder(root, cbtrav, trav_size, outbuffer);
  else if (traversal == TREE_TRAVERSE_PREORDER)
    rtree_traverse_preorder(root, cbtrav, trav_size, outbuffer);
  else
    fatal("Invalid traversal value.");

  return BPP_SUCCESS;
}

static void rnode_clone(rnode_t * rnode, rnode_t * clone, rtree_t * clone_rtree)
{
  if (clone->label)
    free(clone->label);

  clone->length = rnode->length;
  clone->theta = rnode->theta;
  clone->tau = rnode->tau;
  clone->leaves = rnode->leaves;
  clone->node_index = rnode->node_index;

  if (rnode->attrib)
    clone->attrib = xstrdup(rnode->attrib);

  /* points to relatives */
  if (rnode->parent)
    clone->parent = clone_rtree->nodes[rnode->parent->node_index];
  else
    clone->parent = NULL;

  if (rnode->left)
    clone->left = clone_rtree->nodes[rnode->left->node_index];
  else
    clone->left = NULL;

  if (rnode->right)
    clone->right = clone_rtree->nodes[rnode->right->node_index];
  else
    rnode->right = NULL;

  /* label */
  if (rnode->label)
    clone->label = xstrdup(rnode->label);
  else
    clone->label = NULL;

  /* data  - unused */
  clone->data = NULL;
}

static rtree_t * rtree_clone_init(rtree_t * rtree)
{

  unsigned int i;
  unsigned int nodes_count = rtree->tip_count + rtree->inner_count;
  rtree_t * clone;
  
  clone = (rtree_t *)xcalloc(1, sizeof(rtree_t));
  memcpy(clone, rtree, sizeof(rtree_t));

  /* create cloned species tree nodes */
  clone->nodes = (rnode_t **)xmalloc(nodes_count * sizeof(rnode_t *));
  clone->td = (rnode_t **)xmalloc(nodes_count * sizeof(rnode_t *));
  for (i = 0; i < nodes_count; ++i)
    clone->nodes[i] = (rnode_t *)xcalloc(1, sizeof(rnode_t));

  for (i = 0; i < nodes_count; ++i)
    rnode_clone(rtree->nodes[i], clone->nodes[i], clone);

  clone->root = clone->nodes[rtree->root->node_index];

  return clone;
}

rtree_t * rtree_clone(rtree_t * rtree)
{
  unsigned int i;
  unsigned nodes_count = rtree->tip_count + rtree->inner_count;

  rtree_t * clone = rtree_clone_init(rtree);
  /* clone node contents */
  for (i = 0; i < nodes_count; ++i)
    rnode_clone(rtree->nodes[i], clone->nodes[i], clone);

  clone->root = clone->nodes[rtree->root->node_index];

  return clone;
}


static void fill_leaves_and_indices(rnode_t * node,
                                    unsigned int * tip_index,
                                    unsigned int * inner_index)
{
  if (!node->left)
  {
    /* tip node */
    assert(!node->right);

    node->leaves = 1;
    node->node_index = *tip_index;
    (*tip_index)++;
    return;
  }

  /* inner node */
  assert(node->left && node->right);
  fill_leaves_and_indices(node->left,  tip_index, inner_index);
  fill_leaves_and_indices(node->right, tip_index, inner_index);

  node->leaves = node->left->leaves + node->right->leaves;
  node->node_index = *inner_index;
  (*inner_index)++;
}

static void fill_nodes_update_indices(rnode_t * node, rnode_t ** nodes, unsigned int tip_count)
{
  if (!node->left)
  {
    nodes[node->node_index] = node;
    return;
  }

  /* inner node */
  assert(node->left && node->right);
  fill_nodes_update_indices(node->left,  nodes, tip_count);
  fill_nodes_update_indices(node->right, nodes, tip_count);

  node->node_index += tip_count;
  nodes[node->node_index] = node;
  return;
}

rtree_t * rtree_wraptree(rnode_t * root)
{
  rtree_t * tree = (rtree_t *)xcalloc(1,sizeof(rtree_t));

  fill_leaves_and_indices(root, &(tree->tip_count), &(tree->inner_count));
  tree->edge_count = tree->tip_count+tree->inner_count-1;
  assert(tree->tip_count == tree->inner_count+1);

  size_t total_nodes = tree->tip_count + tree->inner_count;
  tree->nodes = (rnode_t **)xmalloc(total_nodes * sizeof(rnode_t *));

  fill_nodes_update_indices(root, tree->nodes, tree->tip_count);

  tree->root = root;

  return tree;
}
