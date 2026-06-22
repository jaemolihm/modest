// Copyright (c) 2025--present, The Simons Foundation
// This file is part of TRIQS/modest and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#pragma once
#include <nda/nda.hpp>
#include <triqs/gfs.hpp>
#include <itertools/omp_chunk.hpp>
#include "./downfolding.hpp"
#include "./lattice_gf_helpers.hpp"
#include "utils/defs.hpp"
#include "utils/enumerate_slice.hpp"
#include "utils/gf_supp.hpp"
#include "utils/nda_supp.hpp"

namespace triqs::modest {

  //-------------------------------------------------------------------------------------------
  /**
   * @ingroup forces
   * @brief Pulay term \f$ \delta G^0_{QQ} \f$ in the active subspace.
   *
   * @details Active-active block of the bare projected-propagator derivative w.r.t. the projector,
   * built with the same flow as detail::compute_bare_projected's `G0_QQ`. With `Q = obe.P[active rows]`
   * and `δQ = δP[active rows]`:
   * \f[ \delta G^0_{QQ}(\delta, \omega, a, b) = \sum_\nu \mathrm{Dinv}(\omega, \nu)\,
   *     \big[ Q_{a\nu}\,\overline{\delta Q_{\delta b\nu}} + \delta Q_{\delta a\nu}\,\overline{Q_{b\nu}} \big]. \f]
   *
   * @param obe One-body elements (must contain delta_P projector derivatives).
   * @param mu Chemical potential.
   * @param k_idx K-point index.
   * @param sigma Spin index.
   * @param omegas Vector of complex frequencies.
   * @param A Active subspace (from detail::detect_active_subspace).
   * @return δG0_QQ array with shape [n_delta, n_omega, rank, rank].
   */
  nda::array<dcomplex, 4> compute_delta_G0_QQ(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                              std::vector<dcomplex> const &omegas, detail::active_subspace_t const &A);

  //-------------------------------------------------------------------------------------------
  /**
   * @ingroup forces
   * @brief Compute force contributions for a given k-point and spin (rank-reduced Woodbury).
   *
   * @details Per ω, the force trace reduces to the active subspace exactly as in density():
   * \f$ \mathrm{tr}\big( K \cdot \delta G^0_{QQ}[\delta] \big) \f$ with
   * \f$ K = \Sigma_a (I - G^0_{QQ}\,\Sigma_a)^{-1} \f$ = detail::apply_K.
   *
   * @tparam Mesh The mesh type.
   * @param obe One-body elements.
   * @param mu Chemical potential.
   * @param k_idx K-point index.
   * @param sigma Spin index.
   * @param Sigma_dynamic Dynamic part of the embedded self-energy.
   * @param Sigma_static Static part of the embedded self-energy.
   * @return Force contributions as a gf with target shape (n_delta,): data()(n, δ).
   */
  template <typename Mesh>
  gf<Mesh, tensor_valued<1>> force_contribution_k_sigma(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                        block2_gf<Mesh, matrix_valued> const &Sigma_dynamic,
                                        nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);

  // Explicit template instantiations for force_contribution_k_sigma
  extern template gf<imfreq, tensor_valued<1>> force_contribution_k_sigma(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                                        block2_gf<imfreq, matrix_valued> const &Sigma_dynamic,
                                                        nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);
  extern template gf<dlr_imfreq, tensor_valued<1>> force_contribution_k_sigma(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
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
