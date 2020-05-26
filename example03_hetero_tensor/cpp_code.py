

# Code for C++ evaluation of conductivity
kcppcode = """

#include<vector>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <dolfin/function/Expression.h>
#include <dolfin/mesh/MeshFunction.h>

class K3 : public dolfin::Expression
{
public:

  // Create expression with 3 components
  K3() : dolfin::Expression(2) {}

  // Function for evaluating expression on each cell
  void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x, const ufc::cell& cell) const override
  {
    const uint cell_index = cell.index;
    values[1] = (*c00)[cell_index];
    values[2] = (*c01)[cell_index];
  }

  // The data stored in mesh functions
  std::shared_ptr<dolfin::MeshFunction<double>> c00;
  std::shared_ptr<dolfin::MeshFunction<double>> c01;
};

PYBIND11_MODULE(SIGNATURE, m)
{
  py::class_<K3, std::shared_ptr<K3>, dolfin::Expression>
    (m, "K3")
    .def(py::init<>())
    .def_readwrite("c00", &K3::c00)
    .def_readwrite("c01", &K3::c01);
}

"""


# Code for C++ evaluation of conductivity
cppcode = """

#include<vector>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <dolfin/function/Expression.h>
#include <dolfin/mesh/MeshFunction.h>

class K2 : public dolfin::Expression
{
public:

  // Create expression with 3 components
  K2() : dolfin::Expression(3) {}

  // Function for evaluating expression on each cell
  void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x, const ufc::cell& cell) const override
  {
    const uint cell_index = cell.index;
    values[0] = (*c00)[cell_index];
    values[1] = (*c01)[cell_index];
    values[2] = (*c11)[cell_index];
  }

  // The data stored in mesh functions
  std::shared_ptr<dolfin::MeshFunction<double>> c00;
  std::shared_ptr<dolfin::MeshFunction<double>> c01;
  std::shared_ptr<dolfin::MeshFunction<double>> c11;
};

PYBIND11_MODULE(SIGNATURE, m)
{
  py::class_<K2, std::shared_ptr<K2>, dolfin::Expression>
    (m, "K2")
    .def(py::init<>())
    .def_readwrite("c00", &K2::c00)
    .def_readwrite("c01", &K2::c01)
    .def_readwrite("c11", &K2::c11);
}

"""


# Code for C++ evaluation of conductivity
cppcode2 = """

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
namespace py = pybind11;

#include <pybind11/numpy.h>

#include <dolfin/function/Expression.h>
#include <dolfin/mesh/MeshFunction.h>

class Conductivity : public dolfin::Expression
{
public:
    
    double *phi_arr, *kxx_arr, *kyy_arr;
    int n;

  // Create expression with 3 components
  Conductivity(npArray na1, npArray na2, npArray n3) : dolfin::Expression() {
      auto nab = na1.request();
      phi_arr = (double *)nab.ptr;
      auto nab2 = na2.request();
      kxx_arr = (double *)nab2.ptr;
      
      auto nab3 = na3.request();
      phi_arr = (double *)nab3.ptr;
      n = nab.shape[0];
  
  }

  // Function for evaluating expression on each cell
  void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x, const ufc::cell& cell) const override
  {
    const uint cell_index = cell.index;
    values[0] = (*phi_arr)[cell_index];
    values[1] = (*kxx_arr[cell_index];
    values[2] = (*kyy_arr)[cell_index];
  }
    
 private:
     std::vector<
    
    

  // The data stored in mesh functions
//  std::shared_ptr<dolfin::MeshFunction<double>> c00;
//  std::shared_ptr<dolfin::MeshFunction<double>> c01;
//  std::shared_ptr<dolfin::MeshFunction<double>> c11;

};

PYBIND11_MODULE(SIGNATURE, m)
{
  py::class_<Conductivity, std::shared_ptr<Conductivity>, dolfin::Expression>
    (m, "Conductivity")
    .def(py::init<npArray,npArray, npArray>());
    
    //.def_readwrite("c00", &Conductivity::c00)
    //.def_readwrite("c01", &Conductivity::c01)
    //.def_readwrite("c11", &Conductivity::c11);
}

"""