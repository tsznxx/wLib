#include  "Python.h"
#include  "stdio.h"
#include  "math.h"
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
#include  "string.h"
//#include  "malloc.h"

static PyObject *
wRNAfold (PyObject * self, PyObject * args)
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
wRNAcofold (PyObject * self, PyObject * args)
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
wRNAdistance (PyObject * self, PyObject * args)
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

static struct PyMethodDef wRNAMethods[] =
  { {"fold", wRNAfold, 1}, {"cofold", wRNAcofold},{"distance",wRNAdistance}, {NULL, NULL} };
void
initwRNA ()
{
  PyObject *m;
  m = Py_InitModule ("wRNA", wRNAMethods);
}
