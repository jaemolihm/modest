
// C.f. https://numpy.org/doc/1.21/reference/c-api/array.html#importing-the-api
#define PY_ARRAY_UNIQUE_SYMBOL _cpp2py_ARRAY_API
#ifndef CLAIR_C2PY_WRAP_GEN
#ifdef __clang__
// #pragma clang diagnostic ignored "-W#warnings"
#endif
#ifdef __GNUC__
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wcast-function-type"
#pragma GCC diagnostic ignored "-Wcpp"
#endif

#define C2PY_VERSION_MAJOR 0
#define C2PY_VERSION_MINOR 1

#include <c2py/c2py.hpp>

using c2py::operator""_a;

// ==================== Wrapped classes =====================

// ==================== enums =====================

// ==================== module classes =====================

// ==================== module functions ====================

// delta_G0_C_k_sigma
static auto const fun_0 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](const triqs::modest::one_body_elements_on_grid &obe, double mu, long k_idx, long sigma,
                 const std::vector<triqs::dcomplex> &omegas) { return triqs::modest::delta_G0_C_k_sigma(obe, mu, k_idx, sigma, omegas); },
              "obe", "mu", "k_idx", "sigma", "omegas")};

// density
static auto const fun_1 = c2py::dispatcher_f_kw_t{
   c2py::cfun(
      [](const triqs::modest::one_body_elements_on_grid &obe, double mu,
         const triqs::gfs::block2_gf<triqs::mesh::imfreq, triqs::gfs::matrix_valued> &Sigma_dynamic,
         const nda::array<nda::matrix<triqs::dcomplex>, 2> &Sigma_static) { return triqs::modest::density(obe, mu, Sigma_dynamic, Sigma_static); },
      "obe", "mu", "Sigma_dynamic", "Sigma_static"),
   c2py::cfun(
      [](const triqs::modest::one_body_elements_on_grid &obe, double mu,
         const triqs::gfs::block2_gf<triqs::mesh::dlr_imfreq, triqs::gfs::matrix_valued> &Sigma_dynamic,
         const nda::array<nda::matrix<triqs::dcomplex>, 2> &Sigma_static) { return triqs::modest::density(obe, mu, Sigma_dynamic, Sigma_static); },
      "obe", "mu", "Sigma_dynamic", "Sigma_static"),
   c2py::cfun([](const triqs::modest::one_body_elements_tb &obe, double mu,
                 const triqs::gfs::block2_gf<triqs::mesh::imfreq, triqs::gfs::matrix_valued> &Sigma_dynamic,
                 const nda::array<nda::matrix<triqs::dcomplex>, 2> &Sigma_static,
                 const triqs::lattice::bz_int_options &opt) { return triqs::modest::density(obe, mu, Sigma_dynamic, Sigma_static, opt); },
              "obe", "mu", "Sigma_dynamic", "Sigma_static", "opt")};

// density_nk
static auto const fun_2 = c2py::dispatcher_f_kw_t{
   c2py::cfun([](const triqs::modest::one_body_elements_on_grid &obe, double mu, double beta) { return triqs::modest::density_nk(obe, mu, beta); },
              "obe", "mu", "beta")};

// find_chemical_potential
static auto const fun_3 =
   c2py::dispatcher_f_kw_t{
      c2py::cfun([](double target_density, const triqs::modest::one_body_elements_on_grid &obe, double beta, std::string method, double precision,
                    bool verbosity) { return triqs::modest::find_chemical_potential(target_density, obe, beta, method, precision, verbosity); },
                 "target_density", "obe", "beta", "method"_a = "dichotomy", "precision"_a = 1.e-5, "verbosity"_a = true),
      c2py::cfun(
         [](double target_density, const triqs::modest::one_body_elements_on_grid &obe,
            const triqs::gfs::block2_gf<triqs::mesh::imfreq, triqs::gfs::matrix_valued> &Sigma_dynamic,
            const nda::array<nda::matrix<triqs::dcomplex>, 2> &Sigma_static, std::string method, double precision, bool verbosity) {
           return triqs::modest::find_chemical_potential(target_density, obe, Sigma_dynamic, Sigma_static, method, precision, verbosity);
         },
         "target_density", "obe", "Sigma_dynamic", "Sigma_static", "method"_a = "dichotomy", "precision"_a = 1.e-5, "verbosity"_a = true),
      c2py::cfun(
         [](double target_density, const triqs::modest::one_body_elements_on_grid &obe,
            const triqs::gfs::block2_gf<triqs::mesh::dlr_imfreq, triqs::gfs::matrix_valued> &Sigma_dynamic,
            const nda::array<nda::matrix<triqs::dcomplex>, 2> &Sigma_static, std::string method, double precision, bool verbosity) {
           return triqs::modest::find_chemical_potential(target_density, obe, Sigma_dynamic, Sigma_static, method, precision, verbosity);
         },
         "target_density", "obe", "Sigma_dynamic", "Sigma_static", "method"_a = "dichotomy", "precision"_a = 1.e-5, "verbosity"_a = true),
      c2py::cfun(
         [](double target_density, const triqs::modest::one_body_elements_tb &obe,
            const triqs::gfs::block2_gf<triqs::mesh::imfreq, triqs::gfs::matrix_valued> &Sigma_dynamic,
            const nda::array<nda::matrix<triqs::dcomplex>, 2> &Sigma_static, const triqs::lattice::bz_int_options &opt, std::string method,
            double precision, bool verbosity) {
           return triqs::modest::find_chemical_potential(target_density, obe, Sigma_dynamic, Sigma_static, opt, method, precision, verbosity);
         },
         "target_density", "obe", "Sigma_dynamic", "Sigma_static", "opt", "method"_a = "dichotomy", "precision"_a = 1.e-5, "verbosity"_a = true),
      c2py::cfun([](double target_density, const triqs::modest::one_body_elements_tb &obe, const triqs::mesh::imfreq &mesh,
                    const triqs::lattice::bz_int_options &opt, std::string method, double precision,
                    bool verbosity) { return triqs::modest::find_chemical_potential(target_density, obe, mesh, opt, method, precision, verbosity); },
                 "target_density", "obe", "mesh", "opt", "method"_a = "dichotomy", "precision"_a = 1.e-5, "verbosity"_a = true)};

// forces
static auto const fun_4 = c2py::dispatcher_f_kw_t{
   c2py::cfun(
      [](const triqs::modest::one_body_elements_on_grid &obe, double mu,
         const triqs::gfs::block2_gf<triqs::mesh::imfreq, triqs::gfs::matrix_valued> &Sigma_dynamic,
         const nda::array<nda::matrix<triqs::dcomplex>, 2> &Sigma_static) { return triqs::modest::forces(obe, mu, Sigma_dynamic, Sigma_static); },
      "obe", "mu", "Sigma_dynamic", "Sigma_static"),
   c2py::cfun(
      [](const triqs::modest::one_body_elements_on_grid &obe, double mu,
         const triqs::gfs::block2_gf<triqs::mesh::dlr_imfreq, triqs::gfs::matrix_valued> &Sigma_dynamic,
         const nda::array<nda::matrix<triqs::dcomplex>, 2> &Sigma_static) { return triqs::modest::forces(obe, mu, Sigma_dynamic, Sigma_static); },
      "obe", "mu", "Sigma_dynamic", "Sigma_static")};

// trace_forces
static auto const fun_5 =
   c2py::dispatcher_f_kw_t{c2py::cfun(
                              [](const triqs::modest::one_body_elements_on_grid &obe, double mu, long k_idx, long sigma,
                                 const triqs::gfs::block2_gf<triqs::mesh::imfreq, triqs::gfs::matrix_valued> &Sigma_dynamic,
                                 const nda::array<nda::matrix<triqs::dcomplex>, 2> &Sigma_static) {
                                return triqs::modest::trace_forces(obe, mu, k_idx, sigma, Sigma_dynamic, Sigma_static);
                              },
                              "obe", "mu", "k_idx", "sigma", "Sigma_dynamic", "Sigma_static"),
                           c2py::cfun(
                              [](const triqs::modest::one_body_elements_on_grid &obe, double mu, long k_idx, long sigma,
                                 const triqs::gfs::block2_gf<triqs::mesh::dlr_imfreq, triqs::gfs::matrix_valued> &Sigma_dynamic,
                                 const nda::array<nda::matrix<triqs::dcomplex>, 2> &Sigma_static) {
                                return triqs::modest::trace_forces(obe, mu, k_idx, sigma, Sigma_dynamic, Sigma_static);
                              },
                              "obe", "mu", "k_idx", "sigma", "Sigma_dynamic", "Sigma_static")};

static const auto doc_d_0 = fun_0.doc(
   R"DOC(
Compute delta_G Pulay term using P and delta_P for force calculations.

Parameters
----------
obe : {par_0}
   One-body elements (must contain delta_P projector derivatives).
mu : {par_1}
   Chemical potential.
k_idx : {par_2}
   K-point index.
sigma : {par_3}
   Spin index.
omegas : {par_4}
   Vector of complex frequencies.

Returns
-------
{ret_0}
   Delta_G array with shape [n_delta, n_omega, M, M].
)DOC",
   std::vector<std::string>{c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::modest::one_body_elements_on_grid &>()}, ", "),
                            c2py::join(std::vector<std::string>{c2py::python_typename<double>()}, ", "),
                            c2py::join(std::vector<std::string>{c2py::python_typename<long>()}, ", "),
                            c2py::join(std::vector<std::string>{c2py::python_typename<long>()}, ", "),
                            c2py::join(std::vector<std::string>{c2py::python_typename<const std::vector<triqs::dcomplex> &>()}, ", ")},
   std::vector<std::string>{std::vector<std::string>{c2py::python_typename<nda::array<triqs::dcomplex, 4>>()}});
static const auto doc_d_1 = fun_1.doc(
   R"DOC(
[1, 2] Compute the density of the lattice Green's function with a self-energy using Woodbury.

------

[3] Compute the density of the lattice Green's function with a self-energy.

------

Parameters
----------
obe : {par_0}
   One-body elements.
mu : {par_1}
   Chemical potential.
Sigma_dynamic : {par_2}
   Dynamic part of the embedded self-energy.
Sigma_static : {par_3}
   Static part of the embedded self-energy.
opt : {par_4}
   Container for options related to the integration of the BZ.

Returns
-------
{ret_0}
   Electron density of the lattice Green's function.
)DOC",
   std::vector<std::string>{
      c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::modest::one_body_elements_on_grid &>(),
                                          c2py::python_typename<const triqs::modest::one_body_elements_tb &>()},
                 ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<double>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::gfs::block2_gf<triqs::mesh::imfreq, triqs::gfs::matrix_valued> &>(),
                                          c2py::python_typename<const triqs::gfs::block2_gf<triqs::mesh::dlr_imfreq, triqs::gfs::matrix_valued> &>()},
                 ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const nda::array<nda::matrix<triqs::dcomplex>, 2> &>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::lattice::bz_int_options &>()}, ", ")},
   std::vector<std::string>{std::vector<std::string>{c2py::python_typename<double>()}});
static const auto doc_d_2 = fun_2.doc(R"DOC(
Compute number of particles :math:`n = \sum f(\beta(\varepsilon(k) - μ))`.
)DOC",
                                      std::vector<std::string>{}, std::vector<std::string>{});
static const auto doc_d_3 = fun_3.doc(
   R"DOC(
[1] Find the chemical potenital from the local Green's function given a target density.

------

[2, 3, 4, 5] Find the chemical potenital from the local Green's function and self-energy given a target density.

------

Parameters
----------
target_density : {par_0}
   Total electron density.
obe : {par_1}
   One-body elements.
beta : {par_2}
   Inverse temperature (units 1/eV).
method : {par_3}
   Root finding method to use (default = `dichotomy`).
precision : {par_4}
   Precision to end search (default = 1e-5).
verbosity : {par_5}
   Printing of the root finder's progress (default = true).
Sigma_dynamic : {par_6}
   Dynamic part of the embedded self-energy.
Sigma_static : {par_7}
   Static part of the embedded self-energy.
opt : {par_8}
   Container for options related to integration of the BZ.
mesh : {par_9}
   Mesh on which local GF will be computed.

Returns
-------
{ret_0}
   Chemical potential corresponding to target density.
)DOC",
   std::vector<std::string>{
      c2py::join(std::vector<std::string>{c2py::python_typename<double>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::modest::one_body_elements_on_grid &>(),
                                          c2py::python_typename<const triqs::modest::one_body_elements_tb &>()},
                 ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<double>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<std::string>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<double>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<bool>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::gfs::block2_gf<triqs::mesh::imfreq, triqs::gfs::matrix_valued> &>(),
                                          c2py::python_typename<const triqs::gfs::block2_gf<triqs::mesh::dlr_imfreq, triqs::gfs::matrix_valued> &>()},
                 ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const nda::array<nda::matrix<triqs::dcomplex>, 2> &>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::lattice::bz_int_options &>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::mesh::imfreq &>()}, ", ")},
   std::vector<std::string>{std::vector<std::string>{c2py::python_typename<double>()}});
static const auto doc_d_4 = fun_4.doc(
   R"DOC(
Compute force contributions from the lattice Green's function with a self-energy using Woodbury.

Computes the correction term to forces (no Kohn-Sham term included).

Parameters
----------
obe : {par_0}
   One-body elements (must contain delta_P projector derivatives).
mu : {par_1}
   Chemical potential.
Sigma_dynamic : {par_2}
   Dynamic part of the embedded self-energy.
Sigma_static : {par_3}
   Static part of the embedded self-energy.

Returns
-------
{ret_0}
   Force contributions as array of size [n_delta] where n_delta = n_atoms * 3.
)DOC",
   std::vector<std::string>{
      c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::modest::one_body_elements_on_grid &>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<double>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::gfs::block2_gf<triqs::mesh::imfreq, triqs::gfs::matrix_valued> &>(),
                                          c2py::python_typename<const triqs::gfs::block2_gf<triqs::mesh::dlr_imfreq, triqs::gfs::matrix_valued> &>()},
                 ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const nda::array<nda::matrix<triqs::dcomplex>, 2> &>()}, ", ")},
   std::vector<std::string>{std::vector<std::string>{c2py::python_typename<nda::array<double, 1>>()}});
static const auto doc_d_5 = fun_5.doc(
   R"DOC(
Compute force contributions for a given k-point and spin.

Parameters
----------
obe : {par_0}
   One-body elements.
mu : {par_1}
   Chemical potential.
k_idx : {par_2}
   K-point index.
sigma : {par_3}
   Spin index.
Sigma_dynamic : {par_4}
   Dynamic part of the embedded self-energy.
Sigma_static : {par_5}
   Static part of the embedded self-energy.

Returns
-------
{ret_0}
   Force contributions as array[n_delta, n_omega].
)DOC",
   std::vector<std::string>{
      c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::modest::one_body_elements_on_grid &>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<double>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<long>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<long>()}, ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const triqs::gfs::block2_gf<triqs::mesh::imfreq, triqs::gfs::matrix_valued> &>(),
                                          c2py::python_typename<const triqs::gfs::block2_gf<triqs::mesh::dlr_imfreq, triqs::gfs::matrix_valued> &>()},
                 ", "),
      c2py::join(std::vector<std::string>{c2py::python_typename<const nda::array<nda::matrix<triqs::dcomplex>, 2> &>()}, ", ")},
   std::vector<std::string>{std::vector<std::string>{c2py::python_typename<nda::array<triqs::dcomplex, 2>>()}});
//--------------------- module function table  -----------------------------

static PyMethodDef module_methods[] = {
   {"delta_G0_C_k_sigma", (PyCFunction)c2py::pyfkw<fun_0>, METH_VARARGS | METH_KEYWORDS, doc_d_0.c_str()},
   {"density", (PyCFunction)c2py::pyfkw<fun_1>, METH_VARARGS | METH_KEYWORDS, doc_d_1.c_str()},
   {"density_nk", (PyCFunction)c2py::pyfkw<fun_2>, METH_VARARGS | METH_KEYWORDS, doc_d_2.c_str()},
   {"find_chemical_potential", (PyCFunction)c2py::pyfkw<fun_3>, METH_VARARGS | METH_KEYWORDS, doc_d_3.c_str()},
   {"forces", (PyCFunction)c2py::pyfkw<fun_4>, METH_VARARGS | METH_KEYWORDS, doc_d_4.c_str()},
   {"trace_forces", (PyCFunction)c2py::pyfkw<fun_5>, METH_VARARGS | METH_KEYWORDS, doc_d_5.c_str()},
   {nullptr, nullptr, 0, nullptr} // Sentinel
};

//--------------------- module struct & init error definition ------------

//// module doc directly in the code or "" if not present...
/// Or mandatory ?
static struct PyModuleDef module_def = {
   PyModuleDef_HEAD_INIT,
   "rho_and_mu", /* name of module */
   R"RAWDOC(Functions to calculate the total electron density of lattice Green's functions using different methods,
methods for chemical potenital searching given a target density, and force calculations.

Calculate the total electronic density from the lattice Green's function and search for the chemical potential that gives a
target electron count. Compute force contributions from projector derivatives. We provide several algorithms optimized for
speed and robustness. Due to the modularity you can also leverage external root finding routines to write a custom chemical
potential search.
)RAWDOC",        /* module documentation, may be NULL */
   -1,           /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
   module_methods,
   NULL,
   NULL,
   NULL,
   NULL};

//--------------------- module init function -----------------------------

extern "C" __attribute__((visibility("default"))) PyObject *PyInit_rho_and_mu() {

  if (not c2py::check_python_version("rho_and_mu")) return NULL;

  // import numpy iff 'numpy/arrayobject.h' included
#ifdef Py_ARRAYOBJECT_H
  import_array();
#endif

  PyObject *m;

  if (PyType_Ready(&c2py::wrap_pytype<c2py::py_range>) < 0) return NULL;

  m = PyModule_Create(&module_def);
  if (m == NULL) return NULL;

  auto &conv_table = *c2py::conv_table_sptr.get();

  conv_table[std::type_index(typeid(c2py::py_range)).name()] = &c2py::wrap_pytype<c2py::py_range>;

  return m;
}
#endif
// CLAIR_WRAP_GEN
