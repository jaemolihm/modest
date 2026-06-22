// Copyright (c) 2025--present, The Simons Foundation
// This file is part of TRIQS/modest and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <nda/nda.hpp>
#include <triqs/gfs.hpp>
#include <itertools/omp_chunk.hpp>
#include "./downfolding.hpp"
#include "utils/defs.hpp"
#include "utils/enumerate_slice.hpp"
#include "utils/gf_supp.hpp"
#include "utils/nda_supp.hpp"

namespace triqs::modest {

  // JML TODO: The following functions are temporarily exposed to Python API for debugging purposes.
  // JML TODO: These should be moved to detail namespace and made private once debugging is complete.

  //-------------------------------------------------------------------------------------------
  /**
   * @ingroup forces
   * @brief Compute delta_G Pulay term using P and delta_P for force calculations.
   *
   * @param obe One-body elements (must contain delta_P projector derivatives).
   * @param mu Chemical potential.
   * @param k_idx K-point index.
   * @param sigma Spin index.
   * @param omegas Vector of complex frequencies.
   * @return Delta_G array with shape [n_delta, n_omega, M, M].
   */
  nda::array<dcomplex, 4> delta_G0_C_k_sigma(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                              std::vector<dcomplex> const &omegas);

  //-------------------------------------------------------------------------------------------
  /**
   * @ingroup forces
   * @brief Compute force contributions for a given k-point and spin.
   *
   * @tparam Mesh The mesh type.
   * @param obe One-body elements.
   * @param mu Chemical potential.
   * @param k_idx K-point index.
   * @param sigma Spin index.
   * @param Sigma_dynamic Dynamic part of the embedded self-energy.
   * @param Sigma_static Static part of the embedded self-energy.
   * @return Force contributions as array[n_delta, n_omega].
   */
  template <typename Mesh>
  nda::array<dcomplex, 2> trace_forces(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                        block2_gf<Mesh, matrix_valued> const &Sigma_dynamic,
                                        nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);

  // Explicit template instantiations for trace_forces
  extern template nda::array<dcomplex, 2> trace_forces(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                                        block2_gf<imfreq, matrix_valued> const &Sigma_dynamic,
                                                        nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);
  extern template nda::array<dcomplex, 2> trace_forces(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                                        block2_gf<dlr_imfreq, matrix_valued> const &Sigma_dynamic,
                                                        nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);

  //-------------------------------------------------------------------------------------------
  /**
   * @ingroup forces
   * @brief Compute force contributions from the lattice Green's function with a self-energy using Woodbury.
   *
   * Computes the correction term to forces (no Kohn-Sham term included).
   *
   * @tparam Mesh The mesh type.
   * @param obe One-body elements (must contain delta_P projector derivatives).
   * @param mu Chemical potential.
   * @param Sigma_dynamic Dynamic part of the embedded self-energy.
   * @param Sigma_static Static part of the embedded self-energy.
   * @return Force contributions as array of size [n_delta] where n_delta = n_atoms * 3.
   */
  template <typename Mesh>
  nda::array<double, 1> forces(one_body_elements_on_grid const &obe, double mu,
                                block2_gf<Mesh, matrix_valued> const &Sigma_dynamic,
                                nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);

  // Explicit template instantiations
  extern template nda::array<double, 1> forces(one_body_elements_on_grid const &obe, double mu,
                                                block2_gf<imfreq, matrix_valued> const &Sigma_dynamic,
                                                nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);
  extern template nda::array<double, 1> forces(one_body_elements_on_grid const &obe, double mu,
                                                block2_gf<dlr_imfreq, matrix_valued> const &Sigma_dynamic,
                                                nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);

} // namespace triqs::modest
