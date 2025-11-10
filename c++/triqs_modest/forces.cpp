// Copyright (c) 2025--present, The Simons Foundation
// This file is part of TRIQS/modest and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./forces.hpp"
#include "./root_finder.hpp"

namespace triqs::modest {

  // JML TODO: The following functions are temporarily exposed to Python API for debugging purposes.
  // JML TODO: These should be moved to detail namespace and made private once debugging is complete.

  // ============================================
  // Compute delta_G Pulay term using P and delta_P for force calculations
  // Formula: delta_G[delta_i, omega, m, m'] = P * Dinv * delta_P[delta_i]^dagger + delta_P[delta_i] * Dinv * P^dagger
  // Using the same efficient BLAS pattern as G0_C_k_sigma
  nda::array<dcomplex, 4> delta_G0_C_k_sigma(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                              std::vector<dcomplex> const &omegas) {
    auto n_nu = obe.H.N_nu(sigma, k_idx);
    auto R_nu = nda::range(n_nu);
    auto M    = obe.C_space.dim();

    // Get P and delta_P
    auto P       = obe.P.P(sigma, k_idx);       // [M, n_nu]
    auto delta_P = obe.P.delta_P(sigma, k_idx); // optional<[n_delta, M, n_nu]>

    if (!delta_P.has_value()) { throw std::runtime_error("delta_P is not present in the projector. Cannot compute forces."); }

    auto n_delta = delta_P->shape()[0];

    // Compute Dinv(omega, nu) = 1/(om + mu - Hk(nu))
    auto n_omega = omegas.size();
    auto Dinv    = nda::matrix<dcomplex>(n_omega, n_nu);
    auto Hk      = obe.H.H(sigma, k_idx); // Hk[nu,nu']

    for (auto [n, om] : enumerate(omegas))
      for (auto nu : R_nu) Dinv(n, nu) = 1.0 / (om + mu - Hk(nu, nu));

    // Precompute combined product matrix PdP[delta_i, nu, m, m']
    // PdP = P[m, nu] * conj(delta_P[delta_i, m', nu]) + delta_P[delta_i, m, nu] * conj(P[m', nu])
    auto PdP = nda::zeros<dcomplex>(n_delta, n_nu, M, M);
    for (auto delta_i : range(n_delta))
      for (auto nu : R_nu)
        for (auto m : range(M))
          for (auto mp : range(M))
            PdP(delta_i, nu, m, mp) = P(m, nu) * conj((*delta_P)(delta_i, mp, nu)) + (*delta_P)(delta_i, m, nu) * conj(P(mp, nu));

    // Compute dG0[delta_i, omega, m, m'] = Dinv[omega, nu] * PdP[delta_i, nu, m, m']
    auto dG0 = nda::zeros<dcomplex>(n_delta, n_omega, M, M);
    for (auto delta_i : range(n_delta)) {
      // Extract PdP for this delta_i: [nu, m, m']
      auto PdP_i = PdP(delta_i, r_all, r_all, r_all);

      // Temporary storage for this delta_i's result
      auto Y = nda::zeros<dcomplex>(n_omega, M, M);

      // Reshape to [nu, M*M] and [omega, M*M] for BLAS gemm
      auto PdP_asmat = cmat_vt{nda::group_indices_view(PdP_i, nda::idx_group<0>, nda::idx_group<1, 2>)};
      auto Y_asmat   = cmat_vt{nda::group_indices_view(Y, nda::idx_group<0>, nda::idx_group<1, 2>)};

      // Y[omega, m, m'] = Dinv[omega, nu] * PdP_i[nu, m, m']
      nda::blas::gemm(1, Dinv, PdP_asmat, 0, Y_asmat);

      // Store result
      dG0(delta_i, r_all, r_all, r_all) = Y;
    }

    return dG0;
  }

  //-------------------------------------------------------------------------------------------
  // Compute force contributions for a given k-point and spin
  template <typename Mesh>
  nda::array<dcomplex, 2> trace_forces(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                        block2_gf<Mesh, matrix_valued> const &Sigma_dynamic,
                                        nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static) {

    auto M                = obe.C_space.dim();
    auto &mesh            = Sigma_dynamic(0, 0).mesh();
    auto omegas           = mesh | tl::to<std::vector<dcomplex>>();
    auto embedding_decomp = get_struct(Sigma_dynamic).dims(r_all, 0) | tl::to<std::vector>();

    // Get Y1 and Y2
    auto Y1 = detail::G0_C_k_sigma(obe, mu, k_idx, sigma, omegas, false); // Y1 = G0_𝓒 [n_omega, M, M]
    auto Y2 = delta_G0_C_k_sigma(obe, mu, k_idx, sigma, omegas);          // delta_G Pulay term [n_delta, n_omega, M, M]

    auto n_delta = Y2.shape()[0];
    auto result  = nda::zeros<dcomplex>(n_delta, omegas.size());

    for (auto &&[n, om] : itertools::enumerate(mesh)) {
      for (auto delta_i : range(n_delta)) {
        // Compute (1- Y1 Sigma)^{-1} Y2[delta_i]
        // Reuse calc_inv_G_G0 logic from density.hpp
        auto YS = nda::matrix<dcomplex, nda::F_layout>::zeros(M, M);

        for (auto &&[alpha, R] : enumerated_sub_slices(embedding_decomp)) {
          auto Sigma = nda::matrix<dcomplex>{Sigma_dynamic(alpha, sigma)[om] + Sigma_static(alpha, sigma)};
          nda::blas::gemm(-1, Y1(n, r_all, R), Sigma, 0, YS(r_all, R));
        }

        // Z = (1 -YS)^{-1} * Y2[delta_i]
        for (auto m : range(M)) YS(m, m) += 1;
        auto B = nda::matrix<dcomplex, nda::F_layout>{Y2(delta_i, n, r_all, r_all)};
        Ainv_B(YS, B);

        // Tr (Sigma * B)
        dcomplex tr_Sigma_B = 0;
        for (auto &&[alpha, R] : enumerated_sub_slices(embedding_decomp)) {
          auto [m_dim, mp_dim] = Sigma_dynamic(alpha, sigma).target_shape();
          auto A               = Sigma_dynamic(alpha, sigma).data()(n, r_all, r_all) + Sigma_static(alpha, sigma);
          auto C               = B(R, R);
          for (auto m1 : range(m_dim))
            for (auto m2 : range(mp_dim)) tr_Sigma_B += A(m1, m2) * C(m2, m1);
        }
        result(delta_i, n) = tr_Sigma_B;
      }
    }
    return result;
  }

  //-------------------------------------------------------------------------------------------
  /**
   * @brief Compute force contributions from the lattice Green's function with a self-energy using Woodbury.
   */
  template <typename Mesh>
  nda::array<double, 1> forces(one_body_elements_on_grid const &obe, double mu,
                                block2_gf<Mesh, matrix_valued> const &Sigma_dynamic,
                                nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static) {

    // Check for matrix-valued case
    if (obe.H.matrix_valued) { throw std::runtime_error("Matrix-valued Hamiltonian not supported for forces calculation."); }

    // Check if delta_P exists
    auto test_delta_P = obe.P.delta_P(0, 0);
    if (!test_delta_P.has_value()) { throw std::runtime_error("delta_P is not present in the projector. Cannot compute forces."); }

    auto n_delta     = test_delta_P->shape()[0];
    auto n_sigma     = Sigma_dynamic.size2();
    auto n_k         = obe.H.n_k();
    auto const &mesh = Sigma_dynamic(0, 0).mesh();

    // Create n_delta gf objects to accumulate force contributions
    std::vector<gf<Mesh, scalar_valued>> force_gfs;
    for (auto i : range(n_delta)) { force_gfs.push_back(gf{mesh}); }

    // ---------
    mpi::communicator comm = {}; // for now using default comm in MPI

#pragma omp parallel default(none) shared(n_k, comm, n_sigma, obe, mu, Sigma_dynamic, Sigma_static, n_delta, force_gfs, mesh)
    {
      // Thread-local gf objects
      std::vector<gf<Mesh, scalar_valued>> thread_force_gfs;
      for (auto i : range(n_delta)) { thread_force_gfs.push_back(gf{mesh}); }

#pragma omp for collapse(2)
      for (auto k_idx : mpi::chunk(range(n_k), comm)) {
        for (auto sigma : range(n_sigma)) {
          //  Correction term
          auto contrib = trace_forces(obe, mu, k_idx, sigma, Sigma_dynamic, Sigma_static); // [n_delta, n_omega]

          // Accumulate into gf objects with k-weights
          auto k_weight = obe.H.k_weights(k_idx);
          for (auto delta_i : range(n_delta)) {
            for (auto n : range(mesh.size())) { thread_force_gfs[delta_i].data()(n) += k_weight * contrib(delta_i, n); }
          }
        }
      }

// Reduce thread results
#pragma omp critical
      {
        for (auto delta_i : range(n_delta)) { force_gfs[delta_i].data() += thread_force_gfs[delta_i].data(); }
      }
    }

    // MPI reduce and compute final result
    auto result = nda::zeros<double>(n_delta);
    for (auto delta_i : range(n_delta)) {
      force_gfs[delta_i] = mpi::all_reduce(force_gfs[delta_i]);
      result(delta_i)    = real(density(force_gfs[delta_i]));
    }

    return result;
  }

  // ------------------------------------------------------------------------------------
  // Explicit template instantiations
  template nda::array<dcomplex, 2> trace_forces(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                                 block2_gf<imfreq, matrix_valued> const &Sigma_dynamic,
                                                 nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);
  template nda::array<dcomplex, 2> trace_forces(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                                 block2_gf<dlr_imfreq, matrix_valued> const &Sigma_dynamic,
                                                 nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);

  template nda::array<double, 1> forces(one_body_elements_on_grid const &obe, double mu,
                                        block2_gf<imfreq, matrix_valued> const &Sigma_dynamic,
                                        nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);
  template nda::array<double, 1> forces(one_body_elements_on_grid const &obe, double mu,
                                        block2_gf<dlr_imfreq, matrix_valued> const &Sigma_dynamic,
                                        nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static);

} // namespace triqs::modest
