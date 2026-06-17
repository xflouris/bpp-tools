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
#include "pdfgen.h"

//static double doc_height = PDF_MM_TO_POINT(150.0f);
static double doc_height = PDF_A4_WIDTH;
//static double doc_width = PDF_MM_TO_POINT(198.0f);
static double doc_width = PDF_A4_HEIGHT;
static double margin_top = 40;
static double margin_bottom = 60;
static double margin_left = 50;
static double margin_right = 90;
static double padding_top = 20;
static double padding_bottom = 20;
static double padding_left = 0;
static double padding_right = 0;


static coord_t * coord = NULL;

static int cb_trav_full(rnode_t * x)
{
  if (!x->left && !x->right)
    return 0;

  return 1;
}

static void plot_tree_with_grid(struct pdf_doc * pdf,
                                const rtree_t * rtree,
                                double linewidth)
{
  unsigned int i,node_count = 0;
  unsigned int li,ri,xi;
  unsigned int tip_count = 0;
  double PI = 3.14159265358979323846;

  long xtics = 6;
  double xtic_xshift = 3;
  double xtic_yshift = -3;
  double xtic_fontsize = 8;
  double xtic_angle = -5*PI/12;
  double ci_xshift = 0;
  double ci_yshift = 5;
  double ci_lw = 1;
  double fontsize_label = 8;
  double arrow_offset = 0.05;
  double arrow_width = 8;

  double canvas_width = doc_width - margin_left - margin_right;
  double canvas_height = doc_height - margin_bottom - margin_top;
  double treecanvas_height = canvas_height - (padding_top+padding_bottom);
  double treecanvas_width = canvas_width - (padding_left+padding_right);

  nodepinfo_t * rootni = (nodepinfo_t *)(rtree->root->data);
  double maxage_ext = 0.02;
  double maxage = MAX(rootni->age,rootni->hi)*(1+maxage_ext);
  double xscaler = treecanvas_width / maxage;
  double tipspace = treecanvas_height / (rtree->tip_count-1);
  double cap_size = 6;
  double timeline_lw = linewidth/2.;

  char * label = NULL;
  rnode_t * left;
  rnode_t * right;
  rnode_t * x;

  /* TODO: Rescale timeline if left or right padding */
  assert(padding_left == 0);
  assert(padding_right == 0);

  /* draw background */
  pdf_add_filled_rectangle(pdf,
                           NULL,
                           margin_left,
                           margin_bottom,
                           canvas_width,
                           canvas_height,
                           0,
                           PDF_RGB(0xeb,0xeb,0xeb),
                           PDF_BLACK);

  #if 0
  /* draw canvas border */
  pdf_add_rectangle(pdf,
                    NULL,
                    margin_left,
                    margin_bottom,
                    canvas_width,
                    canvas_height,
                    linewidth,
                    PDF_BLACK);
  #endif

  /* draw timeline */
  pdf_add_line(pdf,
               NULL,
               margin_left,
               margin_bottom,
               margin_left+canvas_width,
               margin_bottom,
               timeline_lw,
               PDF_BLACK);

  /* draw left and right cap as independent line segments until I learn how
     to draw a line with a cap in postscript (I expect to die first) */

  /* draw left cap */
  pdf_add_line(pdf,
               NULL,
               margin_left,
               margin_bottom-cap_size/2.,
               margin_left,
               margin_bottom+cap_size/2.,
               timeline_lw,
               PDF_BLACK);

  /* draw right cap */
  pdf_add_line(pdf,
               NULL,
               margin_left+canvas_width,
               margin_bottom-cap_size/2.,
               margin_left+canvas_width,
               margin_bottom+cap_size/2.,
               timeline_lw,
               PDF_BLACK);

  /* show number of xtics */
  double xtic_space = maxage / (xtics+1);
  for (i = 0; i < xtics; ++i)
  {
    /* convert x position canvas coordinate */
    double x = xtic_space*(i+1)*xscaler;

    /* draw xtic marker */
    pdf_add_line(pdf,
                 NULL,
                 margin_left+canvas_width-x,
                 margin_bottom-cap_size/2.,
                 margin_left+canvas_width-x,
                 margin_bottom+cap_size/2.,
                 timeline_lw,
                 PDF_BLACK);

    /* draw xtic value */
    xasprintf(&label,"%.6f",xtic_space*(i+1));
    pdf_add_text_rotate(pdf,
                        NULL,
                        label,
                        xtic_fontsize,
                        margin_left+canvas_width-x+xtic_xshift,
                        margin_bottom+xtic_yshift,
                        xtic_angle,
                        PDF_BLACK);
    free(label);
  }

  /* draw major xtic gridlines */
  for (i = 0; i < xtics; ++i)
  {
    /* convert x position canvas coordinate */
    double x = xtic_space*(i+1)*xscaler;

    pdf_add_line(pdf,
                 NULL,
                 margin_left+canvas_width-x,
                 margin_bottom+cap_size/2.,
                 margin_left+canvas_width-x,
                 margin_bottom+canvas_height,
                 timeline_lw,
                 PDF_RGB(0xf6,0xf6,0xf6));
  }

  /* get nodes in postorder */
  rtree_traverse(rtree->root,
                 TREE_TRAVERSE_POSTORDER,
                 cb_trav_full,
                 rtree->td,
                 &node_count);

  char trunclabel[256];
  size_t truncsize = 7;
  for (i  = 0; i < node_count; ++i)
  {
    /* get trio of nodes and indices */
    x = rtree->td[i];
    left = x->left;
    right = x->right;

    li = left->node_index;
    ri = right->node_index;
    xi = x->node_index;

    /* if left or right child is a tip, calculate its coordinates */
    if (li < rtree->tip_count)
    {
      if (strlen(rtree->nodes[li]->label) > 7)
      {
        memcpy(trunclabel,rtree->nodes[li]->label,truncsize);
        trunclabel[truncsize+0] = '.';
        trunclabel[truncsize+1] = '.';
        trunclabel[truncsize+2] = '.';
        trunclabel[truncsize+3] = '\0';
      }
      else
      {
        strcpy(trunclabel,rtree->nodes[li]->label);
      }

      coord[li].x = doc_width - margin_right;
      coord[li].y = margin_bottom+padding_bottom + tipspace * tip_count++;

      char * labelindex;
      xasprintf(&labelindex,  "%s (%d)", trunclabel, li);
      pdf_add_text_rotate(pdf, NULL, labelindex, fontsize_label,
                          coord[li].x+5, coord[li].y-5, 0, PDF_BLACK);
      free(labelindex);
    }
    if (ri < rtree->tip_count)
    {
      if (strlen(rtree->nodes[ri]->label) > 7)
      {
        memcpy(trunclabel,rtree->nodes[ri]->label,truncsize);
        trunclabel[truncsize+0] = '.';
        trunclabel[truncsize+1] = '.';
        trunclabel[truncsize+2] = '.';
        trunclabel[truncsize+3] = '\0';
      }
      else
      {
        strcpy(trunclabel,rtree->nodes[ri]->label);
      }

      coord[ri].x = doc_width - margin_right;
      coord[ri].y = margin_bottom+padding_bottom + tipspace * tip_count++;

      char * labelindex;
      xasprintf(&labelindex,  "%s (%d)", trunclabel, ri);
      pdf_add_text_rotate(pdf, NULL, labelindex, fontsize_label,
                          coord[ri].x+5, coord[ri].y-5, 0, PDF_BLACK);
      free(labelindex);
    }

    /* calculate xy coordinates for node x */
    coord[xi].x = margin_left + canvas_width - x->tau * xscaler;
    if (coord[ri].y > coord[li].y)
      coord[xi].y = coord[li].y + (coord[ri].y - coord[li].y) / 2.0f;
    else
      coord[xi].y = coord[ri].y + (coord[li].y - coord[ri].y) / 2.0f;

    /* draw the three line segments around an inner node */
    pdf_add_line(pdf, NULL, coord[li].x, coord[li].y,
                 coord[xi].x, coord[li].y, linewidth, PDF_BLACK);        
    pdf_add_line(pdf, NULL, coord[ri].x, coord[ri].y,
                 coord[xi].x, coord[ri].y, linewidth, PDF_BLACK);        
    pdf_add_line(pdf, NULL, coord[xi].x, coord[li].y,
                 coord[xi].x, coord[ri].y, linewidth, PDF_BLACK);        

    /* print inner node label (node index) */
    xasprintf(&label, "%ld", x->node_index);
    if (x->parent)
      pdf_add_text(pdf, NULL, label, fontsize_label,
                   coord[xi].x+5, coord[xi].y-5, PDF_BLACK);
    else
      pdf_add_text(pdf, NULL, label, fontsize_label,
                   coord[xi].x-13, coord[xi].y-5, PDF_BLACK);
    free(label);
  }

  #if 0
  /* draw 95% HPD interval */
  for (i  = 0; i < node_count; ++i)
  {
    /* get trio of nodes and indices */
    x = rtree->td[i];
    xi = x->node_index;

    /* draw 95% HPD interval */
    nodepinfo_t * ni = (nodepinfo_t *)(x->data);
    double lox = margin_left + canvas_width - ni->lo*xscaler;
    double hix = margin_left + canvas_width - ni->hi*xscaler;
    pdf_add_line(pdf,
                 NULL,
                 lox+ci_xshift,
                 coord[xi].y+ci_yshift,
                 hix+ci_xshift,
                 coord[xi].y+ci_yshift,
                 ci_lw,
                 PDF_BLUE);
    pdf_add_line(pdf,
                 NULL,
                 lox+ci_xshift,
                 coord[xi].y+ci_yshift+cap_size/2.,
                 lox+ci_xshift,
                 coord[xi].y+ci_yshift+cap_size/2.,
                 ci_lw,
                 PDF_BLUE);

  }
  #endif

  /* print max value on timeline */
  xasprintf(&label, "%.6f", maxage );
  pdf_add_text_rotate(pdf,
                      NULL,
                      label,
                      xtic_fontsize,
                      margin_left+xtic_xshift,
                      margin_bottom+xtic_yshift,
                      xtic_angle,
                      PDF_BLACK);
  free(label);


  /* print 0 (present) on timeline */
  xasprintf(&label, "0");
  pdf_add_text_rotate(pdf,
                      NULL,
                      label,
                      xtic_fontsize,
                      margin_left+canvas_width+xtic_xshift,
                      margin_bottom+xtic_yshift,
                      0,
                      PDF_BLACK);
  free(label);

}

void rtree_export_pdf(const rtree_t * rtree, const char * outfile)
{
  unsigned int i,total_nodes;
  struct pdf_doc * pdf;

  struct pdf_info info = {
      .creator = "bpp" "4.7.1",
      .producer = "bpp" "4.7.1",
      .title = "Binary species tree",
      .author = "Tomas Flouri",
      .subject = "",
      .date = "Today"
  };

  total_nodes = rtree->tip_count + rtree->inner_count;

  coord = (coord_t *)xmalloc((size_t)total_nodes * sizeof(coord_t));

  pdf = pdf_create(doc_width, doc_height, &info);
  pdf_append_page(pdf);
  //pdf_set_font(pdf, "Times-Roman");
  pdf_set_font(pdf, "Helvetica");

  for (i = 0; i < rtree->tip_count+rtree->inner_count; ++i)
  {
    rtree->nodes[i]->data = (void *)xcalloc(sizeof(nodepinfo_t),1);

    nodepinfo_t * info = (nodepinfo_t *)(rtree->nodes[i]->data);

    info->age = rtree->nodes[i]->tau;
    info->hi = info->age;
  }

  plot_tree_with_grid(pdf,rtree,3);

  for (i = 0; i < rtree->tip_count+rtree->inner_count; ++i)
   free(rtree->nodes[i]->data);
  /* time signature */
  struct tm * lt = NULL;
  char buffer[256];
  time_t t = time(NULL);
  lt = localtime(&t);
  assert(lt);

  /* strftime(buffer, 256, "%a %b %d %T %Y", lt); */
  strftime(buffer, 256, "%c", lt);

  float buffer_size = 0;
  pdf_get_font_text_width(pdf,"Courier-Bold", buffer, 8, &buffer_size);

  pdf_set_font(pdf, "Courier-Bold");
  pdf_add_text(pdf, NULL, buffer, 8,
               5, 5, PDF_BLACK);
  pdf_add_text(pdf, NULL, cmdline, 8,
               5+buffer_size+5, 5, PDF_BLACK);
  pdf_add_text(pdf,
               NULL, 
               "Created with: bpp-tools " PROG_VERSION,
               8,
               5, doc_height-15, PDF_BLACK);
  
  pdf_get_font_text_width(pdf, "Courier-Bold", PVER_SHA1, 8, &buffer_size);
  pdf_add_text(pdf,
               NULL, 
               PVER_SHA1,
               8,
               doc_width-buffer_size-10, doc_height-15, PDF_BLACK);

  pdf_save(pdf,outfile);
  pdf_destroy(pdf);

  free(coord);
}

