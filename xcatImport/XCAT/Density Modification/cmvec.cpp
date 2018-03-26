#ifdef __CPLUSPLUS__
extern "C" {
#endif

#ifndef __GNUC__
#pragma warning(disable: 4275)
#pragma warning(disable: 4101)

#endif
#include "Python.h"
#include "compile.h"
#include "frameobject.h"
#include <complex>
#include <math.h>
#include <string>
#include "scxx/object.h"
#include "scxx/list.h"
#include "scxx/tuple.h"
#include "scxx/dict.h"
#include <iostream>
#include <stdio.h>
#include "numpy/arrayobject.h"




// global None value for use in functions.
namespace py {
object None = object(Py_None);
}

const char* find_type(PyObject* py_obj)
{
    if(py_obj == NULL) return "C NULL value";
    if(PyCallable_Check(py_obj)) return "callable";
    if(PyString_Check(py_obj)) return "string";
    if(PyInt_Check(py_obj)) return "int";
    if(PyFloat_Check(py_obj)) return "float";
    if(PyDict_Check(py_obj)) return "dict";
    if(PyList_Check(py_obj)) return "list";
    if(PyTuple_Check(py_obj)) return "tuple";
    if(PyFile_Check(py_obj)) return "file";
    if(PyModule_Check(py_obj)) return "module";

    //should probably do more intergation (and thinking) on these.
    if(PyCallable_Check(py_obj) && PyInstance_Check(py_obj)) return "callable";
    if(PyInstance_Check(py_obj)) return "instance";
    if(PyCallable_Check(py_obj)) return "callable";
    return "unknown type";
}

void throw_error(PyObject* exc, const char* msg)
{
 //printf("setting python error: %s\n",msg);
  PyErr_SetString(exc, msg);
  //printf("throwing error\n");
  throw 1;
}

void handle_bad_type(PyObject* py_obj, const char* good_type, const char* var_name)
{
    char msg[500];
    sprintf(msg,"received '%s' type instead of '%s' for variable '%s'",
            find_type(py_obj),good_type,var_name);
    throw_error(PyExc_TypeError,msg);
}

void handle_conversion_error(PyObject* py_obj, const char* good_type, const char* var_name)
{
    char msg[500];
    sprintf(msg,"Conversion Error:, received '%s' type instead of '%s' for variable '%s'",
            find_type(py_obj),good_type,var_name);
    throw_error(PyExc_TypeError,msg);
}


class int_handler
{
public:
    int convert_to_int(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        
        if (!py_obj || !PyInt_Check(py_obj))
            handle_conversion_error(py_obj,"int", name);
        return (int) PyInt_AsLong(py_obj);
    }

    int py_to_int(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !PyInt_Check(py_obj))
            handle_bad_type(py_obj,"int", name);
        
        return (int) PyInt_AsLong(py_obj);
    }
};

int_handler x__int_handler = int_handler();
#define convert_to_int(py_obj,name) \
        x__int_handler.convert_to_int(py_obj,name)
#define py_to_int(py_obj,name) \
        x__int_handler.py_to_int(py_obj,name)


PyObject* int_to_py(PyObject* obj)
{
    return (PyObject*) obj;
}


class float_handler
{
public:
    double convert_to_float(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        
        if (!py_obj || !PyFloat_Check(py_obj))
            handle_conversion_error(py_obj,"float", name);
        return PyFloat_AsDouble(py_obj);
    }

    double py_to_float(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !PyFloat_Check(py_obj))
            handle_bad_type(py_obj,"float", name);
        
        return PyFloat_AsDouble(py_obj);
    }
};

float_handler x__float_handler = float_handler();
#define convert_to_float(py_obj,name) \
        x__float_handler.convert_to_float(py_obj,name)
#define py_to_float(py_obj,name) \
        x__float_handler.py_to_float(py_obj,name)


PyObject* float_to_py(PyObject* obj)
{
    return (PyObject*) obj;
}


class complex_handler
{
public:
    std::complex<double> convert_to_complex(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        
        if (!py_obj || !PyComplex_Check(py_obj))
            handle_conversion_error(py_obj,"complex", name);
        return std::complex<double>(PyComplex_RealAsDouble(py_obj),PyComplex_ImagAsDouble(py_obj));
    }

    std::complex<double> py_to_complex(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !PyComplex_Check(py_obj))
            handle_bad_type(py_obj,"complex", name);
        
        return std::complex<double>(PyComplex_RealAsDouble(py_obj),PyComplex_ImagAsDouble(py_obj));
    }
};

complex_handler x__complex_handler = complex_handler();
#define convert_to_complex(py_obj,name) \
        x__complex_handler.convert_to_complex(py_obj,name)
#define py_to_complex(py_obj,name) \
        x__complex_handler.py_to_complex(py_obj,name)


PyObject* complex_to_py(PyObject* obj)
{
    return (PyObject*) obj;
}


class unicode_handler
{
public:
    Py_UNICODE* convert_to_unicode(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        Py_XINCREF(py_obj);
        if (!py_obj || !PyUnicode_Check(py_obj))
            handle_conversion_error(py_obj,"unicode", name);
        return PyUnicode_AS_UNICODE(py_obj);
    }

    Py_UNICODE* py_to_unicode(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !PyUnicode_Check(py_obj))
            handle_bad_type(py_obj,"unicode", name);
        Py_XINCREF(py_obj);
        return PyUnicode_AS_UNICODE(py_obj);
    }
};

unicode_handler x__unicode_handler = unicode_handler();
#define convert_to_unicode(py_obj,name) \
        x__unicode_handler.convert_to_unicode(py_obj,name)
#define py_to_unicode(py_obj,name) \
        x__unicode_handler.py_to_unicode(py_obj,name)


PyObject* unicode_to_py(PyObject* obj)
{
    return (PyObject*) obj;
}


class string_handler
{
public:
    std::string convert_to_string(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        Py_XINCREF(py_obj);
        if (!py_obj || !PyString_Check(py_obj))
            handle_conversion_error(py_obj,"string", name);
        return std::string(PyString_AsString(py_obj));
    }

    std::string py_to_string(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !PyString_Check(py_obj))
            handle_bad_type(py_obj,"string", name);
        Py_XINCREF(py_obj);
        return std::string(PyString_AsString(py_obj));
    }
};

string_handler x__string_handler = string_handler();
#define convert_to_string(py_obj,name) \
        x__string_handler.convert_to_string(py_obj,name)
#define py_to_string(py_obj,name) \
        x__string_handler.py_to_string(py_obj,name)


               PyObject* string_to_py(std::string s)
               {
                   return PyString_FromString(s.c_str());
               }
               
class list_handler
{
public:
    py::list convert_to_list(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        
        if (!py_obj || !PyList_Check(py_obj))
            handle_conversion_error(py_obj,"list", name);
        return py::list(py_obj);
    }

    py::list py_to_list(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !PyList_Check(py_obj))
            handle_bad_type(py_obj,"list", name);
        
        return py::list(py_obj);
    }
};

list_handler x__list_handler = list_handler();
#define convert_to_list(py_obj,name) \
        x__list_handler.convert_to_list(py_obj,name)
#define py_to_list(py_obj,name) \
        x__list_handler.py_to_list(py_obj,name)


PyObject* list_to_py(PyObject* obj)
{
    return (PyObject*) obj;
}


class dict_handler
{
public:
    py::dict convert_to_dict(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        
        if (!py_obj || !PyDict_Check(py_obj))
            handle_conversion_error(py_obj,"dict", name);
        return py::dict(py_obj);
    }

    py::dict py_to_dict(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !PyDict_Check(py_obj))
            handle_bad_type(py_obj,"dict", name);
        
        return py::dict(py_obj);
    }
};

dict_handler x__dict_handler = dict_handler();
#define convert_to_dict(py_obj,name) \
        x__dict_handler.convert_to_dict(py_obj,name)
#define py_to_dict(py_obj,name) \
        x__dict_handler.py_to_dict(py_obj,name)


PyObject* dict_to_py(PyObject* obj)
{
    return (PyObject*) obj;
}


class tuple_handler
{
public:
    py::tuple convert_to_tuple(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        
        if (!py_obj || !PyTuple_Check(py_obj))
            handle_conversion_error(py_obj,"tuple", name);
        return py::tuple(py_obj);
    }

    py::tuple py_to_tuple(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !PyTuple_Check(py_obj))
            handle_bad_type(py_obj,"tuple", name);
        
        return py::tuple(py_obj);
    }
};

tuple_handler x__tuple_handler = tuple_handler();
#define convert_to_tuple(py_obj,name) \
        x__tuple_handler.convert_to_tuple(py_obj,name)
#define py_to_tuple(py_obj,name) \
        x__tuple_handler.py_to_tuple(py_obj,name)


PyObject* tuple_to_py(PyObject* obj)
{
    return (PyObject*) obj;
}


class file_handler
{
public:
    FILE* convert_to_file(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        Py_XINCREF(py_obj);
        if (!py_obj || !PyFile_Check(py_obj))
            handle_conversion_error(py_obj,"file", name);
        return PyFile_AsFile(py_obj);
    }

    FILE* py_to_file(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !PyFile_Check(py_obj))
            handle_bad_type(py_obj,"file", name);
        Py_XINCREF(py_obj);
        return PyFile_AsFile(py_obj);
    }
};

file_handler x__file_handler = file_handler();
#define convert_to_file(py_obj,name) \
        x__file_handler.convert_to_file(py_obj,name)
#define py_to_file(py_obj,name) \
        x__file_handler.py_to_file(py_obj,name)


               PyObject* file_to_py(FILE* file, const char* name,
                                    const char* mode)
               {
                   return (PyObject*) PyFile_FromFile(file,
                     const_cast<char*>(name),
                     const_cast<char*>(mode), fclose);
               }
               
class instance_handler
{
public:
    py::object convert_to_instance(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        
        if (!py_obj || !PyInstance_Check(py_obj))
            handle_conversion_error(py_obj,"instance", name);
        return py::object(py_obj);
    }

    py::object py_to_instance(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !PyInstance_Check(py_obj))
            handle_bad_type(py_obj,"instance", name);
        
        return py::object(py_obj);
    }
};

instance_handler x__instance_handler = instance_handler();
#define convert_to_instance(py_obj,name) \
        x__instance_handler.convert_to_instance(py_obj,name)
#define py_to_instance(py_obj,name) \
        x__instance_handler.py_to_instance(py_obj,name)


PyObject* instance_to_py(PyObject* obj)
{
    return (PyObject*) obj;
}


class numpy_size_handler
{
public:
    void conversion_numpy_check_size(PyArrayObject* arr_obj, int Ndims,
                                     const char* name)
    {
        if (arr_obj->nd != Ndims)
        {
            char msg[500];
            sprintf(msg,"Conversion Error: received '%d' dimensional array instead of '%d' dimensional array for variable '%s'",
                    arr_obj->nd,Ndims,name);
            throw_error(PyExc_TypeError,msg);
        }
    }

    void numpy_check_size(PyArrayObject* arr_obj, int Ndims, const char* name)
    {
        if (arr_obj->nd != Ndims)
        {
            char msg[500];
            sprintf(msg,"received '%d' dimensional array instead of '%d' dimensional array for variable '%s'",
                    arr_obj->nd,Ndims,name);
            throw_error(PyExc_TypeError,msg);
        }
    }
};

numpy_size_handler x__numpy_size_handler = numpy_size_handler();
#define conversion_numpy_check_size x__numpy_size_handler.conversion_numpy_check_size
#define numpy_check_size x__numpy_size_handler.numpy_check_size


class numpy_type_handler
{
public:
    void conversion_numpy_check_type(PyArrayObject* arr_obj, int numeric_type,
                                     const char* name)
    {
        // Make sure input has correct numeric type.
        int arr_type = arr_obj->descr->type_num;
        if (PyTypeNum_ISEXTENDED(numeric_type))
        {
        char msg[80];
        sprintf(msg, "Conversion Error: extended types not supported for variable '%s'",
                name);
        throw_error(PyExc_TypeError, msg);
        }
        if (!PyArray_EquivTypenums(arr_type, numeric_type))
        {

        const char* type_names[23] = {"bool", "byte", "ubyte","short", "ushort",
                                "int", "uint", "long", "ulong", "longlong", "ulonglong",
                                "float", "double", "longdouble", "cfloat", "cdouble",
                                "clongdouble", "object", "string", "unicode", "void", "ntype",
                                "unknown"};
        char msg[500];
        sprintf(msg,"Conversion Error: received '%s' typed array instead of '%s' typed array for variable '%s'",
                type_names[arr_type],type_names[numeric_type],name);
        throw_error(PyExc_TypeError,msg);
        }
    }

    void numpy_check_type(PyArrayObject* arr_obj, int numeric_type, const char* name)
    {
        // Make sure input has correct numeric type.
        int arr_type = arr_obj->descr->type_num;
        if (PyTypeNum_ISEXTENDED(numeric_type))
        {
        char msg[80];
        sprintf(msg, "Conversion Error: extended types not supported for variable '%s'",
                name);
        throw_error(PyExc_TypeError, msg);
        }
        if (!PyArray_EquivTypenums(arr_type, numeric_type))
        {
            const char* type_names[23] = {"bool", "byte", "ubyte","short", "ushort",
                                    "int", "uint", "long", "ulong", "longlong", "ulonglong",
                                    "float", "double", "longdouble", "cfloat", "cdouble",
                                    "clongdouble", "object", "string", "unicode", "void", "ntype",
                                    "unknown"};
            char msg[500];
            sprintf(msg,"received '%s' typed array instead of '%s' typed array for variable '%s'",
                    type_names[arr_type],type_names[numeric_type],name);
            throw_error(PyExc_TypeError,msg);
        }
    }
};

numpy_type_handler x__numpy_type_handler = numpy_type_handler();
#define conversion_numpy_check_type x__numpy_type_handler.conversion_numpy_check_type
#define numpy_check_type x__numpy_type_handler.numpy_check_type


class numpy_handler
{
public:
    PyArrayObject* convert_to_numpy(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        Py_XINCREF(py_obj);
        if (!py_obj || !PyArray_Check(py_obj))
            handle_conversion_error(py_obj,"numpy", name);
        return (PyArrayObject*) py_obj;
    }

    PyArrayObject* py_to_numpy(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !PyArray_Check(py_obj))
            handle_bad_type(py_obj,"numpy", name);
        Py_XINCREF(py_obj);
        return (PyArrayObject*) py_obj;
    }
};

numpy_handler x__numpy_handler = numpy_handler();
#define convert_to_numpy(py_obj,name) \
        x__numpy_handler.convert_to_numpy(py_obj,name)
#define py_to_numpy(py_obj,name) \
        x__numpy_handler.py_to_numpy(py_obj,name)


PyObject* numpy_to_py(PyObject* obj)
{
    return (PyObject*) obj;
}


class catchall_handler
{
public:
    py::object convert_to_catchall(PyObject* py_obj, const char* name)
    {
        // Incref occurs even if conversion fails so that
        // the decref in cleanup_code has a matching incref.
        
        if (!py_obj || !(py_obj))
            handle_conversion_error(py_obj,"catchall", name);
        return py::object(py_obj);
    }

    py::object py_to_catchall(PyObject* py_obj, const char* name)
    {
        // !! Pretty sure INCREF should only be called on success since
        // !! py_to_xxx is used by the user -- not the code generator.
        if (!py_obj || !(py_obj))
            handle_bad_type(py_obj,"catchall", name);
        
        return py::object(py_obj);
    }
};

catchall_handler x__catchall_handler = catchall_handler();
#define convert_to_catchall(py_obj,name) \
        x__catchall_handler.convert_to_catchall(py_obj,name)
#define py_to_catchall(py_obj,name) \
        x__catchall_handler.py_to_catchall(py_obj,name)


PyObject* catchall_to_py(PyObject* obj)
{
    return (PyObject*) obj;
}



static PyObject* pyramidvol(PyObject*self, PyObject* args, PyObject* kywds)
{
    py::object return_val;
    int exception_occurred = 0;
    PyObject *py_local_dict = NULL;
    static const char *kwlist[] = {"ax","ay","az","bx","by","bz","cx","cy","cz","local_dict", NULL};
    PyObject *py_ax, *py_ay, *py_az, *py_bx, *py_by, *py_bz, *py_cx, *py_cy, *py_cz;
    int ax_used, ay_used, az_used, bx_used, by_used, bz_used, cx_used, cy_used, cz_used;
    py_ax = py_ay = py_az = py_bx = py_by = py_bz = py_cx = py_cy = py_cz = NULL;
    ax_used= ay_used= az_used= bx_used= by_used= bz_used= cx_used= cy_used= cz_used = 0;
    
    if(!PyArg_ParseTupleAndKeywords(args,kywds,"OOOOOOOOO|O:pyramidvol",const_cast<char**>(kwlist),&py_ax, &py_ay, &py_az, &py_bx, &py_by, &py_bz, &py_cx, &py_cy, &py_cz, &py_local_dict))
       return NULL;
    try                              
    {                                
        py_ax = py_ax;
        double ax = convert_to_float(py_ax,"ax");
        ax_used = 1;
        py_ay = py_ay;
        double ay = convert_to_float(py_ay,"ay");
        ay_used = 1;
        py_az = py_az;
        double az = convert_to_float(py_az,"az");
        az_used = 1;
        py_bx = py_bx;
        double bx = convert_to_float(py_bx,"bx");
        bx_used = 1;
        py_by = py_by;
        double by = convert_to_float(py_by,"by");
        by_used = 1;
        py_bz = py_bz;
        double bz = convert_to_float(py_bz,"bz");
        bz_used = 1;
        py_cx = py_cx;
        double cx = convert_to_float(py_cx,"cx");
        cx_used = 1;
        py_cy = py_cy;
        double cy = convert_to_float(py_cy,"cy");
        cy_used = 1;
        py_cz = py_cz;
        double cz = convert_to_float(py_cz,"cz");
        cz_used = 1;
        /*<function call here>*/     
        
            return_val=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
        if(py_local_dict)                                  
        {                                                  
            py::dict local_dict = py::dict(py_local_dict); 
        }                                                  
    
    }                                
    catch(...)                       
    {                                
        return_val =  py::object();      
        exception_occurred = 1;       
    }                                
    /*cleanup code*/                     
    if(!(PyObject*)return_val && !exception_occurred)
    {
                                  
        return_val = Py_None;            
    }
                                  
    return return_val.disown();           
}                                
static PyObject* displaced_volume_map_floor(PyObject*self, PyObject* args, PyObject* kywds)
{
    py::object return_val;
    int exception_occurred = 0;
    PyObject *py_local_dict = NULL;
    static const char *kwlist[] = {"marr","out","nx","ny","nz","local_dict", NULL};
    PyObject *py_marr, *py_out, *py_nx, *py_ny, *py_nz;
    int marr_used, out_used, nx_used, ny_used, nz_used;
    py_marr = py_out = py_nx = py_ny = py_nz = NULL;
    marr_used= out_used= nx_used= ny_used= nz_used = 0;
    
    if(!PyArg_ParseTupleAndKeywords(args,kywds,"OOOOO|O:displaced_volume_map_floor",const_cast<char**>(kwlist),&py_marr, &py_out, &py_nx, &py_ny, &py_nz, &py_local_dict))
       return NULL;
    try                              
    {                                
        py_marr = py_marr;
        PyArrayObject* marr_array = convert_to_numpy(py_marr,"marr");
        conversion_numpy_check_type(marr_array,PyArray_DOUBLE,"marr");
        #define MARR1(i) (*((double*)(marr_array->data + (i)*Smarr[0])))
        #define MARR2(i,j) (*((double*)(marr_array->data + (i)*Smarr[0] + (j)*Smarr[1])))
        #define MARR3(i,j,k) (*((double*)(marr_array->data + (i)*Smarr[0] + (j)*Smarr[1] + (k)*Smarr[2])))
        #define MARR4(i,j,k,l) (*((double*)(marr_array->data + (i)*Smarr[0] + (j)*Smarr[1] + (k)*Smarr[2] + (l)*Smarr[3])))
        npy_intp* Nmarr = marr_array->dimensions;
        npy_intp* Smarr = marr_array->strides;
        int Dmarr = marr_array->nd;
        double* marr = (double*) marr_array->data;
        marr_used = 1;
        py_out = py_out;
        PyArrayObject* out_array = convert_to_numpy(py_out,"out");
        conversion_numpy_check_type(out_array,PyArray_DOUBLE,"out");
        #define OUT1(i) (*((double*)(out_array->data + (i)*Sout[0])))
        #define OUT2(i,j) (*((double*)(out_array->data + (i)*Sout[0] + (j)*Sout[1])))
        #define OUT3(i,j,k) (*((double*)(out_array->data + (i)*Sout[0] + (j)*Sout[1] + (k)*Sout[2])))
        #define OUT4(i,j,k,l) (*((double*)(out_array->data + (i)*Sout[0] + (j)*Sout[1] + (k)*Sout[2] + (l)*Sout[3])))
        npy_intp* Nout = out_array->dimensions;
        npy_intp* Sout = out_array->strides;
        int Dout = out_array->nd;
        double* out = (double*) out_array->data;
        out_used = 1;
        py_nx = py_nx;
        int nx = convert_to_int(py_nx,"nx");
        nx_used = 1;
        py_ny = py_ny;
        int ny = convert_to_int(py_ny,"ny");
        ny_used = 1;
        py_nz = py_nz;
        int nz = convert_to_int(py_nz,"nz");
        nz_used = 1;
        /*<function call here>*/     
        
            int x,y,z;
            int dx,dy,dz;
            int xp,yp,zp;
            double total;
            
            double ax,ay,az,bx,by,bz,cx,cy,cz;
               
            for (x=0;x<int(nx)-1;x++) 
                for (y=0;y<int(ny)-1;y++) 
                    for (z=0;z<int(nz)-1;z++){
                        total=0;
                        xp=x+1;
                        yp=y+1;
                        zp=z+1;
                        
                        ax=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); //a=P0-P3
                        ay=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]);
                        az=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);
                        
                        //bx=marr(xp,y,z,0)-marr(x,y,zp,0); //b=P1-P3
                        //by=marr(xp,y,z,1)-marr(x,y,zp,1);
                        //bz=marr(xp,y,z,2)-marr(x,y,zp,2);
                        bx=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                        by=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                        bz=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]); 
                        
                        
                        //cx=marr(x,yp,z,0)-marr(x,y,zp,0); //c=P2-P3
                        //cy=marr(x,yp,z,1)-marr(x,y,zp,1);
                        //cz=marr(x,yp,z,2)-marr(x,y,zp,2);   
                        cx=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                        cy=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                        cz=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]); 
                        
                        
                        total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
                        
                        //ax=marr(xp,yp,z,0)-marr(x,y,zp,0); //a=P5-P3
                        //ay=marr(xp,yp,z,1)-marr(x,y,zp,1);
                        //az=marr(xp,yp,z,2)-marr(x,y,zp,2);
                        ax=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                        ay=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                        az=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]); 
                        
                        //b=P1-P3 already set
                        //c=P2-P3 already set
                        total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
                        
                        //a=P5-P3
                        //bx=marr(x,yp,zp,0)-marr(x,y,zp,0); //b=P4-P3
                        //by=marr(x,yp,zp,1)-marr(x,y,zp,1);
                        //bz=marr(x,yp,zp,2)-marr(x,y,zp,2);
                        bx=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                        by=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                        bz=*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]); 
                        //c=P2-P3
                        total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
                        
                        //ax=marr(xp,yp,z,0)-marr(x,yp,zp,0); //a=P5-P4
                        //ay=marr(xp,yp,z,1)-marr(x,yp,zp,1);
                        //az=marr(xp,yp,z,2)-marr(x,yp,zp,2);
                        ax=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                        ay=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                        az=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);                                             
                        
                        
                        //bx=marr(xp,y,zp,0)-marr(x,yp,zp,0); //b=P6-P4
                        //by=marr(xp,y,zp,1)-marr(x,yp,zp,1);
                        //bz=marr(xp,y,zp,2)-marr(x,yp,zp,2);
                        bx=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                        by=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                        bz=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);                                             
                                                                      
                        //cx=marr(xp,yp,zp,0)-marr(x,yp,zp,0); //c=P7-P4
                        //cy=marr(xp,yp,zp,1)-marr(x,yp,zp,1);
                        //cz=marr(xp,yp,zp,2)-marr(x,yp,zp,2);                                                              
                        cx=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                        cy=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                        cz=*(double*)(marr_array->data+xp*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);                                             
                        
                        total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
                        
                        //a=P5-P4
                        //b=P6-P4
                        //cx=marr(xp,y,z,0)-marr(x,yp,zp,0); //c=P1-P4
                        //cy=marr(xp,y,z,1)-marr(x,yp,zp,1);
                        //cz=marr(xp,y,z,2)-marr(x,yp,zp,2);  
                        cx=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                        cy=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                        cz=*(double*)(marr_array->data+xp*marr_array->strides[0]+y*marr_array->strides[1]+z*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);                                             
                        
                        total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
                        
                        //ax=marr(x,y,zp,0)-marr(x,yp,zp,0); //a=P3-P4
                        //ay=marr(x,y,zp,1)-marr(x,yp,zp,1);
                        //az=marr(x,y,zp,2)-marr(x,yp,zp,2);
                        ax=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+0*marr_array->strides[3]); 
                        ay=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+1*marr_array->strides[3]); 
                        az=*(double*)(marr_array->data+x*marr_array->strides[0]+y*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3])-*(double*)(marr_array->data+x*marr_array->strides[0]+yp*marr_array->strides[1]+zp*marr_array->strides[2]+2*marr_array->strides[3]);                                             
                        
                        //c=P1-P4
                        //b=P6-P4
                    
                        
                        total+=fabs( (ax*by*cz)-(ax*bz*cy)+(ay*bz*cx)-(ay*bx*cz)+(az*bx*cy)-(az*by*cx))/6;
        
                        //out(x,y,z)=total;
                        *(double *)(out_array->data+x*out_array->strides[0]+y*out_array->strides[1]+z*out_array->strides[2])=total;
        
                    }
                  
        if(py_local_dict)                                  
        {                                                  
            py::dict local_dict = py::dict(py_local_dict); 
        }                                                  
    
    }                                
    catch(...)                       
    {                                
        return_val =  py::object();      
        exception_occurred = 1;       
    }                                
    /*cleanup code*/                     
    if(marr_used)
    {
        Py_XDECREF(py_marr);
        #undef MARR1
        #undef MARR2
        #undef MARR3
        #undef MARR4
    }
    if(out_used)
    {
        Py_XDECREF(py_out);
        #undef OUT1
        #undef OUT2
        #undef OUT3
        #undef OUT4
    }
    if(!(PyObject*)return_val && !exception_occurred)
    {
                                  
        return_val = Py_None;            
    }
                                  
    return return_val.disown();           
}                                
static PyObject* vec_cinverse_multi(PyObject*self, PyObject* args, PyObject* kywds)
{
    py::object return_val;
    int exception_occurred = 0;
    PyObject *py_local_dict = NULL;
    static const char *kwlist[] = {"marr","out","maxiter","parr","npos","local_dict", NULL};
    PyObject *py_marr, *py_out, *py_maxiter, *py_parr, *py_npos;
    int marr_used, out_used, maxiter_used, parr_used, npos_used;
    py_marr = py_out = py_maxiter = py_parr = py_npos = NULL;
    marr_used= out_used= maxiter_used= parr_used= npos_used = 0;
    
    if(!PyArg_ParseTupleAndKeywords(args,kywds,"OOOOO|O:vec_cinverse_multi",const_cast<char**>(kwlist),&py_marr, &py_out, &py_maxiter, &py_parr, &py_npos, &py_local_dict))
       return NULL;
    try                              
    {                                
        py_marr = py_marr;
        PyArrayObject* marr_array = convert_to_numpy(py_marr,"marr");
        conversion_numpy_check_type(marr_array,PyArray_DOUBLE,"marr");
        #define MARR1(i) (*((double*)(marr_array->data + (i)*Smarr[0])))
        #define MARR2(i,j) (*((double*)(marr_array->data + (i)*Smarr[0] + (j)*Smarr[1])))
        #define MARR3(i,j,k) (*((double*)(marr_array->data + (i)*Smarr[0] + (j)*Smarr[1] + (k)*Smarr[2])))
        #define MARR4(i,j,k,l) (*((double*)(marr_array->data + (i)*Smarr[0] + (j)*Smarr[1] + (k)*Smarr[2] + (l)*Smarr[3])))
        npy_intp* Nmarr = marr_array->dimensions;
        npy_intp* Smarr = marr_array->strides;
        int Dmarr = marr_array->nd;
        double* marr = (double*) marr_array->data;
        marr_used = 1;
        py_out = py_out;
        PyArrayObject* out_array = convert_to_numpy(py_out,"out");
        conversion_numpy_check_type(out_array,PyArray_DOUBLE,"out");
        #define OUT1(i) (*((double*)(out_array->data + (i)*Sout[0])))
        #define OUT2(i,j) (*((double*)(out_array->data + (i)*Sout[0] + (j)*Sout[1])))
        #define OUT3(i,j,k) (*((double*)(out_array->data + (i)*Sout[0] + (j)*Sout[1] + (k)*Sout[2])))
        #define OUT4(i,j,k,l) (*((double*)(out_array->data + (i)*Sout[0] + (j)*Sout[1] + (k)*Sout[2] + (l)*Sout[3])))
        npy_intp* Nout = out_array->dimensions;
        npy_intp* Sout = out_array->strides;
        int Dout = out_array->nd;
        double* out = (double*) out_array->data;
        out_used = 1;
        py_maxiter = py_maxiter;
        int maxiter = convert_to_int(py_maxiter,"maxiter");
        maxiter_used = 1;
        py_parr = py_parr;
        PyArrayObject* parr_array = convert_to_numpy(py_parr,"parr");
        conversion_numpy_check_type(parr_array,PyArray_DOUBLE,"parr");
        #define PARR1(i) (*((double*)(parr_array->data + (i)*Sparr[0])))
        #define PARR2(i,j) (*((double*)(parr_array->data + (i)*Sparr[0] + (j)*Sparr[1])))
        #define PARR3(i,j,k) (*((double*)(parr_array->data + (i)*Sparr[0] + (j)*Sparr[1] + (k)*Sparr[2])))
        #define PARR4(i,j,k,l) (*((double*)(parr_array->data + (i)*Sparr[0] + (j)*Sparr[1] + (k)*Sparr[2] + (l)*Sparr[3])))
        npy_intp* Nparr = parr_array->dimensions;
        npy_intp* Sparr = parr_array->strides;
        int Dparr = parr_array->nd;
        double* parr = (double*) parr_array->data;
        parr_used = 1;
        py_npos = py_npos;
        int npos = convert_to_int(py_npos,"npos");
        npos_used = 1;
        /*<function call here>*/     
        
            int i,p;        
            double x,y,z,f_x,f_y,f_z,pos_x,pos_y,pos_z;
            int xf,yf,zf,xp,yp,zp;
            int warned;    
            double weight;    
              
            for (p=0; p < int(npos);p++){
                //x=guesses(p,0);
                //y=guesses(p,1);
                //z=guesses(p,2);
                x=(*(double*)(parr_array->data+p*parr_array->strides[0]+0*parr_array->strides[1]));
                y=(*(double*)(parr_array->data+p*parr_array->strides[0]+1*parr_array->strides[1]));
                z=(*(double*)(parr_array->data+p*parr_array->strides[0]+2*parr_array->strides[1]));   
                pos_x=x;
                pos_y=y;
                pos_z=z;
                
                for (i=0;i<maxiter;i++){
                
                    /* Calculate forward vector using linear interpolation */
            
                    xf=(int)floor(x);
                    yf=(int)floor(y);
                    zf=(int)floor(z);      
                    f_x=0;
                    f_y=0;
                    f_z=0;
                    for (xp=0;xp<=1;xp++)
                        for (yp=0;yp<=1;yp++)
                            for (zp=0;zp<=1;zp++) {
                                weight=1.0;
                                if (xp==0)
                                    weight*=(1-(x-xf));
                                else
                                    weight*=(x-xf);
                                if (yp==0)
                                    weight*=(1-(y-yf));
                                else
                                    weight*=(y-yf);
                                if (zp==0)
                                    weight*=(1-(z-zf));
                                else
                                    weight*=(z-zf);
                                   
                                if (! ( (xf+xp >= marr_array->dimensions[0]) || (xf+xp <0) || (yf+yp >= marr_array->dimensions[1]) || (yf+yp <0) || (zf+zp >= marr_array->dimensions[2]) || (zf+zp <0) )) {
                                    f_x+=(*(double *)(marr_array->data+(xf+xp)*marr_array->strides[0]+(yf+yp)*marr_array->strides[1]+(zf+zp)*marr_array->strides[2]+0*marr_array->strides[3])-(xf+xp))*weight;
                                    f_y+=(*(double *)(marr_array->data+(xf+xp)*marr_array->strides[0]+(yf+yp)*marr_array->strides[1]+(zf+zp)*marr_array->strides[2]+1*marr_array->strides[3])-(yf+yp))*weight;
                                    f_z+=(*(double *)(marr_array->data+(xf+xp)*marr_array->strides[0]+(yf+yp)*marr_array->strides[1]+(zf+zp)*marr_array->strides[2]+2*marr_array->strides[3])-(zf+zp))*weight;
                                }
                            }
        
                    /* Okay, we've got the forward vector*/
                
                    x=pos_x-f_x;
                    y=pos_y-f_y;
                    z=pos_z-f_z;
                    
                }
                *(double*)(out_array->data+p*out_array->strides[0]+0*out_array->strides[1])=-1.0*f_x;
                *(double*)(out_array->data+p*out_array->strides[0]+1*out_array->strides[1])=-1.0*f_y;
                *(double*)(out_array->data+p*out_array->strides[0]+2*out_array->strides[1])=-1.0*f_z;
            }
        if(py_local_dict)                                  
        {                                                  
            py::dict local_dict = py::dict(py_local_dict); 
        }                                                  
    
    }                                
    catch(...)                       
    {                                
        return_val =  py::object();      
        exception_occurred = 1;       
    }                                
    /*cleanup code*/                     
    if(marr_used)
    {
        Py_XDECREF(py_marr);
        #undef MARR1
        #undef MARR2
        #undef MARR3
        #undef MARR4
    }
    if(out_used)
    {
        Py_XDECREF(py_out);
        #undef OUT1
        #undef OUT2
        #undef OUT3
        #undef OUT4
    }
    if(parr_used)
    {
        Py_XDECREF(py_parr);
        #undef PARR1
        #undef PARR2
        #undef PARR3
        #undef PARR4
    }
    if(!(PyObject*)return_val && !exception_occurred)
    {
                                  
        return_val = Py_None;            
    }
                                  
    return return_val.disown();           
}                                


static PyMethodDef compiled_methods[] = 
{
    {"pyramidvol",(PyCFunction)pyramidvol , METH_VARARGS|METH_KEYWORDS},
    {"displaced_volume_map_floor",(PyCFunction)displaced_volume_map_floor , METH_VARARGS|METH_KEYWORDS},
    {"vec_cinverse_multi",(PyCFunction)vec_cinverse_multi , METH_VARARGS|METH_KEYWORDS},
    {NULL,      NULL}        /* Sentinel */
};

PyMODINIT_FUNC initcmvec(void)
{
    
    Py_Initialize();
    import_array();
    PyImport_ImportModule("numpy");
    (void) Py_InitModule("cmvec", compiled_methods);
}

#ifdef __CPLUSCPLUS__
}
#endif
