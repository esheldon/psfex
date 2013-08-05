/*
   No error checking on the eigens array is done, do that in python
   Make sure it is native byte order and C order.
*/
#include <Python.h>
#include "psfex.h"
#include <numpy/arrayobject.h> 

struct PyPSFExObject {
  PyObject_HEAD
  struct psfex* psfex;
};
static void
PyPSFExObject_dealloc(struct PyPSFExObject* self)
{
    self->psfex = psfex_free(self->psfex);
    self->ob_type->tp_free((PyObject*)self);
}

// do your type checking in python!
static void copy_eigens(struct psfex *self, PyObject *eigens_obj)
{
    long size=0;
    double *data = (double* )PyArray_DATA(eigens_obj);
    int i;

    size=PSFEX_SIZE_TOT(self);
    memcpy(self->eigens->rows[0], data, size);

    //fprintf(stdout,"Copying eigens... [%d]\n",size);
    for (i=0;i<size;i++) {
      self->maskcomp[i] = (float) data[i];
    }
}

static int
PyPSFExObject_init(struct PyPSFExObject* self, PyObject *args)
{
    long neigen=0;
    long nrow=0;   // per eigen
    long ncol=0;   // per eigen
    long poldeg=0;
    double polzero_row=0;
    double polzero_col=0;
    double polscale_row=0;
    double polscale_col=0;
    double psf_samp=0;
    PyObject *eigens_obj=NULL;

    if (!PyArg_ParseTuple(args, 
                          (char*)"lllldddddO", 
                          &neigen,
                          &nrow,
                          &ncol,
                          &poldeg,
                          &polzero_row,
                          &polzero_col,
                          &polscale_row,
                          &polscale_col,
                          &psf_samp,
                          &eigens_obj)) {
        printf("failed to Parse init");
        return -1;
    }

    self->psfex=psfex_new(neigen,
                          nrow,   // per eigen
                          ncol,   // per eigen
                          poldeg,
                          polzero_row,
                          polzero_col,
                          polscale_row,
                          polscale_col,
                          psf_samp);
    if (!self->psfex) {
        return -1;
    }

    copy_eigens(self->psfex, eigens_obj);
    return 0;
}

static PyObject *
PyPSFExObject_repr(struct PyPSFExObject* self) {
    /*
    char repr[255];
    if (self->cosmo != NULL) {
        sprintf(repr, "flat:    %d\n"
                      "DH:      %f\n"
                      "omega_m: %f\n" 
                      "omega_l: %f\n" 
                      "omega_k: %f", 
                      self->cosmo->flat, 
                      self->cosmo->DH, 
                      self->cosmo->omega_m, 
                      self->cosmo->omega_l, 
                      self->cosmo->omega_k);
        return PyString_FromString(repr);
    }  else {
        return PyString_FromString("");
    }
    */
    return PyString_FromString("");
}


static PyObject *make_psf_image(const struct psfex *self)
{
    PyObject *image=NULL;
    int ndims=2;
    npy_intp dims[2];
    dims[0] = PSFEX_NROW(self);
    dims[1] = PSFEX_NCOL(self);
    image = PyArray_ZEROS(ndims, dims, NPY_FLOAT64, 0);
    return image;
}

PyObject *PyPSFExObject_rec(struct PyPSFExObject* self, PyObject *args)
{
    double row=0, col=0;
    if (!PyArg_ParseTuple(args, (char*)"dd", &row, &col)) {
        return NULL;
    }

    PyObject *image=make_psf_image(self->psfex);
    double *data=(double*) PyArray_DATA(image);

    _psfex_rec_fill(self->psfex, row, col, data);

    return image;
}



static PyMethodDef PyPSFExObject_methods[] = {
    {"rec",   (PyCFunction)PyPSFExObject_rec,  METH_VARARGS, "rec\n\nrec(row,col)"},
    {NULL}  /* Sentinel */
};


static PyTypeObject PyPSFExType = {
    PyObject_HEAD_INIT(NULL)
    0,                         /*ob_size*/
    "_psfex_pywrap.PSFEx",             /*tp_name*/
    sizeof(struct PyPSFExObject), /*tp_basicsize*/
    0,                         /*tp_itemsize*/
    (destructor)PyPSFExObject_dealloc, /*tp_dealloc*/
    0,                         /*tp_print*/
    0,                         /*tp_getattr*/
    0,                         /*tp_setattr*/
    0,                         /*tp_compare*/
    //0,                         /*tp_repr*/
    (reprfunc)PyPSFExObject_repr,                         /*tp_repr*/
    0,                         /*tp_as_number*/
    0,                         /*tp_as_sequence*/
    0,                         /*tp_as_mapping*/
    0,                         /*tp_hash */
    0,                         /*tp_call*/
    0,                         /*tp_str*/
    0,                         /*tp_getattro*/
    0,                         /*tp_setattro*/
    0,                         /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    "PSFEx Class",           /* tp_doc */
    0,                     /* tp_traverse */
    0,                     /* tp_clear */
    0,                     /* tp_richcompare */
    0,                     /* tp_weaklistoffset */
    0,                     /* tp_iter */
    0,                     /* tp_iternext */
    PyPSFExObject_methods,             /* tp_methods */
    0,             /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    //0,     /* tp_init */
    (initproc)PyPSFExObject_init,      /* tp_init */
    0,                         /* tp_alloc */
    //PyPSFExObject_new,                 /* tp_new */
    PyType_GenericNew,                 /* tp_new */
};


static PyMethodDef PSFEx_type_methods[] = {
    {NULL}  /* Sentinel */
};


#ifndef PyMODINIT_FUNC  /* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
init_psfex_pywrap(void) 
{
    PyObject* m;

    PyPSFExType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&PyPSFExType) < 0)
        return;

    m = Py_InitModule3("_psfex_pywrap", PSFEx_type_methods, "Define PSFEx type and methods.");

    Py_INCREF(&PyPSFExType);
    PyModule_AddObject(m, "PSFEx", (PyObject *)&PyPSFExType);

    import_array();
}

