/*
   No error checking on the psf_mask array is done, do that in python
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

// do type checking in python!
static void copy_psf_mask(struct psfex *self, PyObject *psf_mask_obj)
{
    long size=0;
    double *data = (double *) PyArray_DATA(psf_mask_obj);

    size=PSFEX_SIZE(self)*PSFEX_NCOMP(self)*sizeof(double);

    memcpy(self->maskcomp, data, size);
}


static int
PyPSFExObject_init(struct PyPSFExObject* self, PyObject *args)
{
    PyObject *masksize_obj = NULL;
    long poldeg;
    PyObject *contextoffset_obj = NULL;
    PyObject *contextscale_obj = NULL;
    double psf_samp;
    PyObject *psf_mask_obj = NULL;

    long *masksize;
    double *contextoffset, *contextscale;

    // can totally redo this object-wise.

    if (!PyArg_ParseTuple(args,
			  (char*)"OlOOdO",
			  &masksize_obj,
			  &poldeg,
			  &contextoffset_obj,
			  &contextscale_obj,
			  &psf_samp,
			  &psf_mask_obj)) {
        printf("failed to Parse init");
        return -1;
    }

    masksize = (long *) PyArray_DATA(masksize_obj);
    contextoffset = (double *) PyArray_DATA(contextoffset_obj);
    contextscale = (double *) PyArray_DATA(contextscale_obj);

    self->psfex=psfex_new(masksize,
			  poldeg,
			  contextoffset,
			  contextscale,
			  psf_samp);

    if (!self->psfex) {
        return -1;
    }

    copy_psf_mask(self->psfex, psf_mask_obj);

    return 0;
}

static PyObject *
PyPSFExObject_repr(struct PyPSFExObject* self) {
    return PyString_FromString("");
}


static PyObject *make_psf_image(const struct psfex *self)
{
    PyObject *image=NULL;
    int ndims=2;
    npy_intp dims[2];
    //dims[0] = PSFEX_NROW(self);
    //dims[1] = PSFEX_NCOL(self);
    dims[0] = RECON_NROW(self);
    dims[1] = RECON_NCOL(self);
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

PyObject *PyPSFExObject_center(struct PyPSFExObject *self, PyObject *args)
{
    PyObject *retval = NULL;
    double row=0, col=0;
    double rowcen, colcen;
    
    int ndims=1;
    npy_intp dims[1];
    double *data;

    
    if (!PyArg_ParseTuple(args, (char*)"dd", &row, &col)) {
	return NULL;
    }

    get_center(RECON_NROW(self->psfex), RECON_NCOL(self->psfex),
	       row, col,
	       self->psfex->pixstep,
	       &rowcen, &colcen);

    dims[0] = 2;
    retval = PyArray_SimpleNew(ndims, dims, NPY_FLOAT64);

    data = (double *) PyArray_DATA(retval);
    data[0] = rowcen;
    data[1] = colcen;

    return retval;
}

static PyMethodDef PyPSFExObject_methods[] = {
    {"rec",   (PyCFunction)PyPSFExObject_rec,  METH_VARARGS, "rec\n\nrec(row,col)"},
    {"center", (PyCFunction)PyPSFExObject_center, METH_VARARGS, "center\n\ncenter(row,col)"},
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

