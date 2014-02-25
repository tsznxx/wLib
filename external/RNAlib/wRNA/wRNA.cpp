/*****************************************************************************
  wRNA.cpp
  Last-modified: 21 Oct 2013 04:31:32 PM

  (c) 2012 - Yunfei Wang
  Center for Systems Biology
  Department of Molecular & Cell Biology
  University of Texas at Dallas
  tszn1984@gmail.com

  Licensed under the GNU General Public License 2.0 license.
******************************************************************************/

#include  "Python.h"
#include <iostream>
#include <string>
#include <cstring>
#include <cmath>

#include "wRNA.h"

// RNAStructure headers
#include "ParseCommandLine.h"
#include "StructureImageHandler.h"


// Vienna RNA Package
#include  "utils.h"
#include  "fold_vars.h"
#include  "fold.h"
#include  "cofold.h"
#include  "part_func.h"
#include  "inverse.h"
#include  "RNAstruct.h"
#include  "treedist.h"
#include  "stringdist.h"
#include  "profiledist.h"
#include  "params.h"
#include "PS_dot.h"

static PyObject *
wfold (PyObject * self, PyObject * args)
{
  char *seq, *structure, *cstruct;
  float e;
  float myTemperature = 37.0;
  int myFold_constrained = 0;
  PyObject *ret;
  if (!PyArg_ParseTuple (args, "ssf", &seq, &cstruct, &myTemperature))
    {
      return Py_BuildValue("");;
    }
  temperature = myTemperature;
  initialize_fold (strlen (seq));
  structure = (char *) space (sizeof (char) * (strlen (seq) + 1));
  if (strlen (cstruct) > 1)
    {
      strncpy (structure, cstruct, strlen (seq));
      myFold_constrained = 1;
    }
  fold_constrained = myFold_constrained;
  e = fold (seq, structure);
  ret = Py_BuildValue ("(f,s)", e, structure);
  free(structure);
  free_arrays ();
  return ret;
}

static PyObject *
wcofold (PyObject * self, PyObject * args)
{
  char *seq, *seqA, *seqB, *f;
  float e;
  float myTemperature = 37.0;
  PyObject *ret;
  if (!PyArg_ParseTuple (args, "ss|f", &seqA, &seqB, &myTemperature))
    {
      return Py_BuildValue("");;
    }
  temperature = myTemperature;
  cut_point = strlen (seqA) + 1;
  seq = (char *) space (sizeof (char) * (strlen (seqA) + strlen (seqB) + 1));
  strcpy (seq, seqA);
  strcat (seq, ".");
  strcat (seq, seqB);
  initialize_cofold (strlen (seq));
  f = (char *) space (sizeof (char) * (strlen (seq) + 1));
  e = cofold (seq, f);
  ret = Py_BuildValue ("(f,s)", e / 100, f);
  free_arrays ();
  free (seq);
  free (f);
  return ret;
}

static PyObject *
wdistance (PyObject * self, PyObject * args)
{
  char *structure1, *structure2, *xstruc1, *xstruc2;
  Tree *T1, *T2;
  float tree_dist;
  PyObject *ret;
  if (!PyArg_ParseTuple (args, "ss", &structure1, &structure2))
    {
      return Py_BuildValue("");;
    }
  xstruc1 = expand_Full (structure1);
  xstruc2 = expand_Full (structure2);
  T1 = make_tree (xstruc1);
  T2 = make_tree (xstruc2);
  //S1 = Make_swString (xstruc1);
  //S2 = Make_swString (xstruc2);
  free(xstruc1);
  free(xstruc2);
  edit_backtrack = 1;
  tree_dist = tree_edit_distance (T1, T2);
  free(T1);
  free(T2);
  ret = Py_BuildValue("f",tree_dist);
  return ret;
}

/*
 * Input:
 *   1. Seq
 *   2. structure in dot-bracket format
 *   3. outfile name
 *   4. PS or SVG. 0: PS, 1: SVG
 * Usage:
 *	 wRNA.dotplot("UUAACAUGAUGAAAAAAUAUAUUAACACAGACCUGUACUGAACUUUUCGAAGUUUUGCAGAUAACAAUAUUGCUUUUUUUCUCUGACU",".......((((........))))...(((....)))...((((((....))))))..((((.((.((........)).)).))))...", "test.ps", 0)
 * Output:
 *   test.ps or test.svg
 */

static PyObject *
wdotplot (PyObject * self, PyObject * args)
{
	char *seq, *structure, *outfile;
	int   isnSVG;
	if (!PyArg_ParseTuple (args, "sssi", &seq, &structure, &outfile, &isnSVG)) return Py_BuildValue("s","1");
	if (isnSVG == 0) 
		PS_rna_plot_a(seq, structure, outfile, NULL, NULL);
	else 
		svg_rna_plot(seq, structure, outfile);
	return Py_BuildValue("s","0");
}

static PyObject *
wplot (PyObject * self, PyObject * args)
{
    char *infile, *outfile;
    bool isSVG = true;    
	int  num, isnSVG;
    string rst;
    if (!PyArg_ParseTuple (args, "ssii", &infile, &outfile, &isnSVG, &num)) return Py_BuildValue("s","1");
    if (isnSVG == 0) isSVG = false;
	StructureImageHandler* handler = new StructureImageHandler();
    // Read file
    handler->readRadial( infile, num);
    // plot
    if(strcmp(outfile, "") == 0)
    {
      if (isSVG)
          rst = handler->writeSVG();
      else
          rst = handler->writePostscript();
      delete handler;
      return Py_BuildValue("s", rst.c_str());
    }

    if (isSVG)
      handler->writeSVG(outfile);
    else
      handler->writePostscript(outfile);
    delete handler;
    return Py_BuildValue("s","0");
}


static struct PyMethodDef wRNAMethods[] =
  { {"fold", wfold, 1}, {"cofold", wcofold},{"distance",wdistance}, {"dotplot",wdotplot}, {"plot", wplot}, {NULL, NULL} };

extern "C" 
{
	void initwRNA ()
	{
		PyObject *m;
		m = Py_InitModule ("wRNA", wRNAMethods);
	}
}
