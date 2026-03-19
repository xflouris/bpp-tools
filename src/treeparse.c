/*
    Copyright (C) 2016-2022 Tomas Flouri, Bruce Rannala and Ziheng Yang

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

/* tokens returned by lexical analyzer */
#define TOKEN_NONE              0
#define TOKEN_ATTR              1
#define TOKEN_OPAR              2
#define TOKEN_CPAR              3
#define TOKEN_COMMA             4
#define TOKEN_COLON             5
#define TOKEN_SEMICOLON         6
#define TOKEN_HASH              7
#define TOKEN_STRING            8

const char * token_names[9] = {
  "Illegal", "TOKEN_ATTR", "TOKEN_OPAR", "TOKEN_CPAR", "TOKEN_COMMA",
  "TOKEN_COLON", "TOKEN_SEMICOLON", "TOKEN_HASH", "TOKEN_STRING"
};

static void node_destroy(node_t * root, void (*cb_data_destroy)(void *));

const unsigned int attrib_map[256] = 
 {
/* 0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F  */

   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 0 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 1 */    
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 2 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 3 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 4 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  0,  1,  1,  /* 5 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 6 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  /* 7 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 8 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 9 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* A */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* B */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* C */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* D */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* E */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* F */
 };

const unsigned int string_map[256] = 
 {
/* 0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F  */

   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 0 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 1 */    
   0,  1,  0,  1,  1,  1,  1,  0,  0,  0,  1,  1,  0,  1,  1,  1,  /* 2 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  0,  1,  1,  1,  1,  /* 3 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 4 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  1,  0,  1,  1,  /* 5 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  /* 6 */
   1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  0,  /* 7 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 8 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* 9 */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* A */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* B */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* C */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* D */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* E */
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  /* F */
 };

/* Syntax parsing table. Entry (i,j) indicates whether token j can follow
   token i (1) or not (0) */
const unsigned int syntax_table[9][9] =
 {
   /* NONE, ATTR, OPAR, CPAR, COMMA, COLON, SEMICOLON, HASH, STRING */
   {  0, 0, 1, 0, 0, 0, 0, 0, 1 },   /* NONE */
   {  0, 0, 0, 1, 1, 1, 1, 1, 0 },   /* ATTR */
   {  0, 0, 1, 0, 0, 0, 0, 0, 1 },   /* OPAR */
   {  0, 1, 0, 1, 1, 1, 1, 1, 1 },   /* CPAR */
   {  0, 0, 1, 0, 0, 0, 0, 0, 1 },   /* COMMA */
   {  0, 0, 0, 0, 0, 0, 0, 0, 1 },   /* COLON */
   {  0, 0, 0, 0, 0, 0, 0, 0, 0 },   /* SEMICOLON */
   {  0, 0, 0, 0, 0, 0, 0, 0, 1 },   /* HASH */
   {  0, 1, 0, 1, 1, 1, 1, 1, 0 },   /* STRING */
 };

typedef struct ltoken_s
{
  char * data;
  int type;
} ltoken_t;


static long parse_attr(char * s, list_t * token_list)
{
  long i = 0;
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));

  assert(*s == '[');
  ++s;

  /* match string */
  while (attrib_map[(int)(s[i])]) ++i;

  if (s[i] != ']')
    fatal("[ERROR] Cannot match closing bracket in attribute: %s", s-1);

  token->data = xstrndup(s-1,i+2);
  token->type = TOKEN_ATTR;

  list_append(token_list,(void *)token);

  /* return length of token */
  return i+2;
}

static void token_clear(void * tokenptr)
{
  ltoken_t * token = (ltoken_t *)tokenptr;
  if (token)
  {
    if (token->data)
      free(token->data);
    free(token);
  }
}

long parse_opar(list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup("(");
  token->type = TOKEN_OPAR;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_cpar(list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup(")");
  token->type = TOKEN_CPAR;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_colon(list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup(":");
  token->type = TOKEN_COLON;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_semicolon(list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup(";");
  token->type = TOKEN_SEMICOLON;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_comma(list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup(",");
  token->type = TOKEN_COMMA;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_theta(list_t * token_list)
{
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  token->data = xstrdup("#");
  token->type = TOKEN_HASH;

  list_append(token_list,(void *)token);

  /* return length of token */
  return 1;
}

static long parse_string(char * s, list_t * token_list)
{
  long i = 0;
  ltoken_t * token = (ltoken_t *)xmalloc(sizeof(ltoken_t));
  
  /* match string */
  while (string_map[(int)(s[i])]) ++i;

  token->data = xstrndup(s,i);
  token->type = TOKEN_STRING;

  list_append(token_list,(void *)token);

  return i;
}

static long get_double(const char * line, double * value)
{
  int ret,len=0;
  size_t ws;
  char * s = xstrdup(line);
  char * p = s;

  /* skip all white-space */
  ws = strspn(p, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!p[ws] || p[ws] == '*' || p[ws] == '#')
  {
    free(s);
    return 0;
  }

  /* store address of value's beginning */
  char * start = p+ws;

  /* skip all characters except star, hash and whitespace */
  char * end = start + strcspn(start," \t\r\n*#");

  *end = 0;

  ret = sscanf(start, "%lf%n", value, &len);
  if ((ret == 0) || (((unsigned int)(len)) < strlen(start)))
  {
    free(s);
    return 0;
  }

  ret = ws + end - start;
  free(s);
  return ret;
}

static void dealloc_data(rnode_t * node,
                         void (*cb_destroy)(void *))
{
  if (node->data)
  {
    if (cb_destroy)
      cb_destroy(node->data);
  }
}

void rtree_destroy(rtree_t * tree,
                   void (*cb_destroy)(void *))
{
  unsigned int i;
  rnode_t * node;

  /* deallocate all nodes */
  for (i = 0; i < tree->tip_count + tree->inner_count; ++i)
  {
    node = tree->nodes[i];
    dealloc_data(node,cb_destroy);

    if (node->label)
      free(node->label);

    free(node);
  }

  if (tree->td)
    free(tree->td);

  /* deallocate tree structure */
  free(tree->nodes);
  free(tree);
}

int opt_precision = 6;
static char * ntree_export_newick_recursive(node_t * root, long print_bl)
{
  int i;
  char * newick;
  char * x;

  if (!root) return NULL;

  if (!root->children_count)
  {
    if (print_bl)
      xasprintf(&newick,
                "%s%s:%.*f",
                root->label,
                root->attr ? root->attr : "",
                opt_precision,
                root->length);
    else
      xasprintf(&newick,
                "%s%s",
                root->label,
                root->attr ? root->attr : "");
  }
  else
  {
    char * subtree = ntree_export_newick_recursive(root->children[0],print_bl);
    xasprintf(&newick, "(%s", subtree);
    free(subtree);
    for (i = 1; i < root->children_count; ++i)
    {
      subtree = ntree_export_newick_recursive(root->children[i],print_bl);
      xasprintf(&x, "%s,%s", newick, subtree);
      free(newick);
      free(subtree);
      newick = x;
    }
    if (print_bl)
      xasprintf(&x,
                "%s)%s%s:%.*f",
                newick,
                root->label ? root->label : "",
                root->attr ? root->attr : "", 
                opt_precision,
                root->length);
    else
      xasprintf(&x,
                "%s)%s%s",
                newick,
                root->label ? root->label : "",
                root->attr ? root->attr : "");
    free(newick);
    newick = x;
  }

  return newick;
}

char * ntree_export_newick(ntree_t * tree, long print_bl)
{
  int i;
  char * newick;
  char * x;

  node_t * root = tree->root;

  if (!root) return NULL;

  if (!root->children_count)
  {
    if (print_bl)
      xasprintf(&newick,
                "%s%s:%.*f",
                root->label,
                root->attr ? root->attr : "",
                opt_precision,
                root->length);
    else
      xasprintf(&newick,
                "%s%s",
                root->label,
                root->attr ? root->attr : "");
  }
  else
  {
    char * subtree = ntree_export_newick_recursive(root->children[0],print_bl);
    xasprintf(&newick, "(%s", subtree);
    free(subtree);
    for (i = 1; i < root->children_count; ++i)
    {
      subtree = ntree_export_newick_recursive(root->children[i],print_bl);
      xasprintf(&x, "%s,%s", newick, subtree);
      free(newick);
      free(subtree);
      newick = x;
    }
    if (print_bl)
      xasprintf(&x,
                "%s)%s%s:%.*f;",
                newick,
                root->label ? root->label : "",
                root->attr ? root->attr : "",
                opt_precision,
                root->length);
    else
      xasprintf(&x,
                "%s)%s%s",
                newick,
                root->label ? root->label : "",
                root->attr ? root->attr : "");
    free(newick);
    newick = x;
  }

  return newick;
}

list_t * parse_tree(char * s)
{
  int rc = 1;
  long count;
  long opar_count = 0;
  long cpar_count = 0;
  list_t * token_list = NULL;

  /* skip all white-space */
  size_t ws = strspn(s, " \t\r\n");

  /* is it a blank line or comment ? */
  if (!s[ws] || s[ws] == '*' || s[ws] == '#')
  {
    rc = 0;
    goto l_unwind;
  }

  char * start = s+ws;

  token_list = (list_t *)xcalloc(1,sizeof(list_t));

  /* TODO: Find end by going backwards and set a zero char */

  while (*start)
  {
    size_t ws = strspn(start, " \t\r\n");
    start += ws;

    if (*start == '(')
    {
      if (opt_debug_parser)
        printf("Parsing OPAR...\n");
      count = parse_opar(token_list);
      start += count;
    }
    else if (*start == ')')
    {
      if (opt_debug_parser)
        printf("Parsing CPAR...\n");
      count = parse_cpar(token_list);
      start += count;
    }
    else if (*start == '[')
    {
      if (opt_debug_parser)
        printf("Parsing ATTR...\n");
      count = parse_attr(start,token_list);
      start += count;
    }
    else if (*start == ':')
    {
      if (opt_debug_parser)
        printf("Parsing COLON...\n");
      count = parse_colon(token_list);
      start += count;
    }
    else if (*start == ';')
    {
      if (opt_debug_parser)
        printf("Parsing SEMICOLON...\n");
      count = parse_semicolon(token_list);
      start += count;
    }
    else if (*start == ',')
    {
      if (opt_debug_parser)
        printf("Parsing COMMA...\n");
      count = parse_comma(token_list);
      start += count;
    }
    else if (*start == '#')
    {
      if (opt_debug_parser)
        printf("Parsing HASH...\n");
      count = parse_theta(token_list);
      start += count;
    }
    else if (string_map[(int)(*start)])
    {
      if (opt_debug_parser)
        printf("Parsing STRING...\n");
      /* TODO : eliminate starting with theta  */
      count = parse_string(start,token_list);
      start += count;
    }
    else
    {
      fatal("Illegal character");
    }

  }
  if (opt_debug_parser)
  {
    printf("Finished\nListing tokens:\n");

    list_item_t * li = token_list->head;
    while (li)
    {
      ltoken_t * token = (ltoken_t *)(li->data);

      printf("%s : %s\n", token_names[token->type], token->data);

      li = li->next;
    }
  }

  
  /* trivial checks */

  /* 1. No tokens means error */
  if (token_list->count == 0)
    fatal("No valid constructs found");


  /* 2. First token must be either OPAR or STRING */
  list_item_t * li = token_list->head;
  ltoken_t * token = (ltoken_t *)(li->data);
  if (token->type != TOKEN_OPAR && token->type != TOKEN_STRING)
    fatal("Newick string must start with '(' or a label");
    
  /* 3. If the first token is not an OPAR then no OPAR and CPAR must exist in
     the list */
  if (token->type != TOKEN_OPAR)
  {
    li = li->next;
    while (li)
    {
      token = (ltoken_t *)(li->data);
      if (token->type == TOKEN_OPAR || token->type == TOKEN_CPAR)
        break;
      li = li->next;
    }
    if (li)
      fatal("Illegal newick format specification - found opening/closing "
            "parenthesis but newick does not start with opening parenthesis");
  }

  /* 4. If first token is OPAR, then check that there is no OPAR when the
     number of open OPARs is 0 */
  li = token_list->head;
  token = (ltoken_t *)(li->data);
  if (token->type == TOKEN_OPAR)
  {
    opar_count++;
    long openpars = 1;
    li = li->next;
    while (li)
    {
      token = (ltoken_t *)(li->data);
      if (token->type == TOKEN_OPAR)
      {
        if (openpars == 0)
        {
          snprintf(bpp_errmsg, 200, "Invalid token type while parsing n-ary tree");
          rc = 0;
          goto l_unwind;
          #if 0
          fatal("Cannot open new parenthesis");
          #endif
        }
        ++openpars;
        ++opar_count;
      }
      else if (token->type == TOKEN_CPAR)
      {
        --openpars;
        ++cpar_count;
      }
      if (openpars < 0)
        fatal("Invalid newick format");

      li = li->next;
    }
  }

  if (opar_count != cpar_count)
  {
    snprintf(bpp_errmsg,
             200,
             "Mismatching number of opening (%ld) and closing (%ld) parentheses)",
             opar_count, cpar_count);
    rc = 0;
    goto l_unwind;
    #if 0
    fatal("Mismatching number of opening (%ld) and closing (%ld) parentheses)",
          opar_count, cpar_count);
    #endif
  }

l_unwind:
  if (!rc)
  {
    if (token_list)
    {
      list_clear(token_list,token_clear);
      free(token_list);
      token_list = NULL;
    }
  }
  return token_list;
}

static void fill_node_lists_recursive(node_t * node,
                                      node_t ** tiplist,
                                      node_t ** innerlist,
                                      int * tipcount,
                                      int * innercount)
{
  int i;

  if (!node) return;

  /* tip node */
  if (!node->children_count)
  {
    tiplist[(*tipcount)++] = node;

    return;
  }

  /* inner node */
  for (i = 0; i < node->children_count; ++i)
  {
    fill_node_lists_recursive(node->children[i],
                              tiplist,
                              innerlist,
                              tipcount,
                              innercount);
  }
  innerlist[(*innercount)++] = node;
}

/* fill the the ntree_t->leaves and ntree_t->inner list of nodes by recursively
   traversing the tree graph starting with the root node */
static void fill_node_lists(node_t * node,
                            node_t ** tiplist,
                            node_t ** innerlist)
{
  int i;
  int tipcount = 0;
  int innercount = 0;

  if (!node) return;

  /* tip node */
  if (!node->children_count)
  {
    tiplist[tipcount++] = node;
    return;
  }

  /* inner nodes */
  for (i = 0; i < node->children_count; ++i)
  {
    fill_node_lists_recursive(node->children[i],
                              tiplist,
                              innerlist,
                              &tipcount,
                              &innercount);
  }
  innerlist[innercount++] = node;
}

ntree_t * ntree_wraptree(node_t * root, int tip_count, int inner_count)
{
  long i;

  ntree_t * tree = (ntree_t *)xcalloc(1,sizeof(ntree_t));

  tree->root = root;
  tree->tip_count = tip_count;
  tree->inner_count = inner_count;

  tree->leaves = (node_t **)xmalloc(tree->tip_count * sizeof(node_t *));
  tree->inner  = (node_t **)xmalloc(tree->inner_count * sizeof(node_t *));

  /* fill node lists in postorder traversal */
  fill_node_lists(tree->root, tree->leaves, tree->inner);

  for (i = 0; i < tree->tip_count; ++i)
    tree->leaves[i]->node_index = i;

  for (i = 0; i < tree->inner_count; ++i)
    tree->inner[i]->node_index = i;

  return tree;
}

static ntree_t * syntax_parse(list_t * token_list)
{
  ltoken_t * token;
  list_item_t * li;
  long prev_token_type = TOKEN_NONE;

  node_t * node = NULL;
  node_t * root = NULL;
  node_t * tmp  = NULL;
  node_t ** children;
  double value;

  int tip_count = 0;
  int inner_count = 0;

  li = token_list->head;

  while (li)
  {
    token = (ltoken_t *)(li->data);

    if (!syntax_table[prev_token_type][token->type])
    {
      #if 0
      printf("%ld %d\n", prev_token_type, token->type);
      #endif
      #if 0
      fatal("Invalid token type while parsing n-ary tree");
      #endif

      snprintf(bpp_errmsg, 200, "Invalid token type while parsing n-ary tree");
      if (root)
        node_destroy(root,NULL);
      return NULL;
    }

    switch (token->type)
    {
      case TOKEN_NONE:
        fatal("Unexpected TOKEN_NONE");
        break;

      case TOKEN_ATTR:
        node->attr = xstrdup(token->data);
        break;

      case TOKEN_OPAR:
        tmp = (node_t *)xcalloc(1,sizeof(node_t));
        tmp->parent = node;
        node = tmp;
        if (!root)
          root = node;
        inner_count++;
        break;

      case TOKEN_CPAR:
        assert(node);
        node = node->parent;

        assert(node);
        if (node->parent)
        {
          node->parent->children_count++;

          if (node->parent->children)
          {
            children = (node_t**)xmalloc((size_t)(node->parent->children_count)*
                                          sizeof(node_t *));
            memcpy(children,node->parent->children,
                   (node->parent->children_count-1)*sizeof(node_t *));
            children[node->parent->children_count-1] = node;
            free(node->parent->children);
            node->parent->children = children;
          }
          else
          {
            assert(node->parent->children_count == 1);
            node->parent->children = (node_t **)xmalloc(sizeof(node_t *));
            node->parent->children[0] = node;
          }
        }
        break;

      case TOKEN_COMMA:
        node = node->parent;
        break;

      case TOKEN_COLON:
        break;

      case TOKEN_SEMICOLON:
        break;

      case TOKEN_HASH:
        break;

      case TOKEN_STRING:
        switch (prev_token_type)
        {
          case TOKEN_NONE:
            assert(node == NULL);
            node = (node_t *)xcalloc(1,sizeof(node_t));
            node->parent = NULL;
            node->label = xstrdup(token->data);
            root = node;
            tip_count++;
            break;

          case TOKEN_OPAR:
            tmp = (node_t *)xcalloc(1,sizeof(node_t));
            tmp->parent = node;
            tmp->label = xstrdup(token->data);
            node = tmp;

            assert(node->parent);
            assert(node->parent->children_count == 0 && !node->parent->children);
            node->parent->children_count = 1;
            node->parent->children = (node_t **)xmalloc(sizeof(node_t *));
            node->parent->children[0] = node;
            tip_count++;
            break;

          case TOKEN_CPAR:
            node->label = xstrdup(token->data);
            break;

          case TOKEN_COMMA:
            tmp = (node_t *)xcalloc(1,sizeof(node_t));
            tmp->parent = node;
            tmp->label = xstrdup(token->data);
            node = tmp;

            node_t * parent = node->parent;
            assert(parent);
            assert(parent->children_count > 0 && parent->children);
            parent->children_count++;

            children = (node_t **)xmalloc((size_t)(parent->children_count) *
                                          sizeof(node_t *));
            memcpy(children,
                   parent->children,
                   (parent->children_count-1)*sizeof(node_t *));
            children[parent->children_count-1] = node;
            free(parent->children);
            parent->children = children;
            tip_count++;
            break;
            
          case TOKEN_COLON:
            if (!get_double(token->data,&value))
              fatal("ERROR: Expected floating point number (branch length) but "
                    "got: %s", token->data);
            if (value < 0)
              fatal("ERROR: Found negative branch length (%f)", value);

            node->length = value;
            break;

          case TOKEN_HASH:
            if (!get_double(token->data,&value))
              fatal("ERROR: Expected floating point number (theta) but got: %s",
                    token->data);
            if (value < 0)
              fatal("ERROR: Found negative theta (%f)", value);

            node->theta = value;
            break;

          default:
            fatal("Internal error when parsing syntax table");

        }
        break;

      default:
        fatal("Unknown token type");
    }

    prev_token_type = token->type;
    li = li->next;
  }

  return ntree_wraptree(root,tip_count,inner_count);
}

static void node_destroy(node_t * root, void (*cb_data_destroy)(void *))
{
  int i;

  if (!root) return;

  if (root->children)
  {
    for (i = 0; i < root->children_count; ++i)
      node_destroy(root->children[i],cb_data_destroy);

    free(root->children);
  }
  
  if (root->data && cb_data_destroy)
    cb_data_destroy(root->data);

  if (root->attr)
    free(root->attr);

  free(root->label);
  free(root);
}


void ntree_destroy(ntree_t * tree, void (*cb_data_destroy)(void *))
{
  if (!tree) return;

  node_destroy(tree->root,cb_data_destroy);

  if (tree->leaves)
    free(tree->leaves);

  if (tree->inner)
    free(tree->inner);

  free(tree);
}

int ntree_check_rbinary(ntree_t * tree)
{
  int i;
  /* checks whether tree is binary rooted */

  assert(!tree->root->parent);

  /* check all inner nodes that they have out-degree 2 */
  for (i = 0; i < tree->inner_count; ++i)
  {
    if (tree->inner[i]->children_count != 2)
      return 0;
  }

  return 1;
}

int ntree_check_ubinary(ntree_t * tree)
{
  int i;
  /* checks whether tree is unrooted binary */

  assert(!tree->root->parent);

  for (i = 0; i < tree->inner_count; ++i)
  {
    if (!tree->inner[i]->parent)
    {
      /* root case */
      if (tree->inner[i]->children_count != 3) return 0;
    }
    else
    {
      /* other inner case */

      if (tree->inner[i]->children_count != 2) return 0;
    }
  }

  return 1;
}

int ntree_check_nary(ntree_t * tree)
{
  int i;
  /* checks whether tree is n-ary */

  assert(!tree->root->parent);

  for (i = 0; i < tree->inner_count; ++i)
  {
    if (tree->inner[i]->children_count > 2) return 1;
  }

  return 0;
}

static rnode_t * rnode_from_node_recursive(node_t * node, rnode_t * parent)
{
  rnode_t * rnode = (rnode_t *)xcalloc(1,sizeof(rnode_t));

  assert(node->children_count <= 2);

  if (node->children_count)
  {
    rnode->left  = rnode_from_node_recursive(node->children[0],rnode);
    if (node->children_count > 1)
      rnode->right = rnode_from_node_recursive(node->children[1],rnode);
    else
      rnode->right = NULL;

  }
  else
  {
    rnode->left = rnode->right = NULL;
  }

  rnode->parent = parent;
  rnode->length = node->length;
  rnode->theta  = node->theta;
  rnode->tau    = node->tau;
  if (node->label)
    rnode->label = xstrdup(node->label);

  if (node->attr)
    rnode->attrib = xstrdup(node->attr);

  return rnode;
}

static void fill_nodes_recursive(rnode_t * rnode,
                                 rnode_t ** array,
                                 unsigned int * tip_index,
                                 unsigned int * inner_index)
{
  if (!rnode->left && !rnode->right)
  {
    array[*tip_index] = rnode;
    *tip_index = *tip_index + 1;
    return;
  }

  array[*inner_index] = rnode;
  *inner_index = *inner_index + 1;

  if (rnode->left)
    fill_nodes_recursive(rnode->left,  array, tip_index, inner_index);
  if (rnode->right)
    fill_nodes_recursive(rnode->right, array, tip_index, inner_index);
}

rtree_t * rtree_from_ntree(ntree_t * ntree)
{
  /* assumes rooted binary ntree */
  unsigned int i;
  unsigned int tip_index = 0;
  unsigned int inner_index = ntree->tip_count;

  rtree_t * rtree = (rtree_t *)xcalloc(1,sizeof(rtree_t));

  /* trivial validation of tree */
  if (!ntree_check_rbinary(ntree))
  {
    if (ntree_check_ubinary(ntree))
    {
      fatal("BPP requires a rooted binary tree as input."
             "Given tree is unrooted binary");
    }
    else
    {
      /* check if tree is n-ary (n>2) */
      if (ntree_check_nary(ntree))
        fatal("BPP requires a rooted binary tree as input."
              "Given tree is multifurcating");
    }
  }

  rtree->root =  rnode_from_node_recursive(ntree->root,NULL);

  rtree->tip_count = ntree->tip_count;
  rtree->inner_count = ntree->inner_count;

  rtree->edge_count = rtree->tip_count+rtree->inner_count-1;

  rtree->nodes = (rnode_t **)xmalloc((size_t)(rtree->tip_count+rtree->inner_count) *
                                     sizeof(rnode_t *));

  fill_nodes_recursive(rtree->root, rtree->nodes, &tip_index, &inner_index);

  for (i = 0; i < rtree->tip_count + rtree->inner_count; ++i)
    rtree->nodes[i]->node_index = i;

  for (i = 0; i < rtree->tip_count + rtree->inner_count; ++i)
    rtree->nodes[i]->node_index = i;

  return rtree;
}

static void rtree_reset_leaves_recursive(rnode_t * rnode)
{
  if (!rnode->left)
  {
    rnode->leaves = 1;
    return;
  }

  if (rnode->left)
    rtree_reset_leaves_recursive(rnode->left);
  if (rnode->right)
    rtree_reset_leaves_recursive(rnode->right);

  rnode->leaves = 0;
  if (rnode->left)
    rnode->leaves = rnode->left->leaves;
  if (rnode->right)
    rnode->leaves += rnode->right->leaves;
}

rtree_t * bpp_parse_newick_string(const char * line)
{
  ntree_t * tree = NULL;
  rtree_t * rtree = NULL;
  list_t * token_list = NULL;

  #if 1
  /* old parser */
  char * s = xstrdup(line);

  if (!(token_list = parse_tree(s)))
    goto l_unwind;

  if (!(tree = syntax_parse(token_list)))
    goto l_unwind;

  rtree = rtree_from_ntree(tree);

  rtree_reset_leaves_recursive(rtree->root);

  rtree->td = (rnode_t **)xmalloc((size_t)(rtree->tip_count+rtree->inner_count) *
                                  sizeof(rnode_t *));

  #else
  /* old parser */

  rtree_t * rtree = rtree_parse_newick_string(line);
  #endif

l_unwind:
  if (s)
    free(s);
  if (token_list)
  {
    list_clear(token_list,token_clear);
    free(token_list);
  }

  if (tree)
    ntree_destroy(tree,NULL);
  
  return rtree;
}

ntree_t * bpp_parse_newick_string_ntree(const char * line)
{
  ntree_t * tree = NULL;
  list_t * token_list = NULL;

  char * s = xstrdup(line);

  if (!(token_list = parse_tree(s)))
    goto l_unwind;

  tree = syntax_parse(token_list);

l_unwind:
  if (token_list)
  {
    list_clear(token_list,token_clear);
    free(token_list);
  }
  free(s);

  return tree;
}
