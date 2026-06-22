// Copyright (c) 2025--present, The Simons Foundation
// This file is part of TRIQS/modest and is licensed under the terms of GPLv3 or later.
// SPDX-License-Identifier: GPL-3.0-or-later
// See LICENSE in the root of this distribution for details.

#include "./forces.hpp"

namespace triqs::modest {

  // ===================================================================================
  // Pulay term δG0_QQ in the active subspace, built with the same flow as
  // detail::compute_bare_projected's G0_QQ (Dinv · index-tensor, one batched gemm per δ).
  // ===================================================================================
  nda::array<dcomplex, 4> delta_G0_C_k_sigma(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                             std::vector<dcomplex> const &omegas, detail::active_subspace_t const &A) {
    auto n_nu          = obe.H.N_nu(sigma, k_idx);
    auto R_nu          = nda::range(n_nu);
    auto rank          = A.rank;
    auto n_omega       = long(omegas.size());
    auto P_obe         = obe.P.P(sigma, k_idx);       // [M, n_nu]
    auto Hk            = obe.H.H(sigma, k_idx);        // Hk[nu, nu']
    auto const &active = A.c_indices;

    auto delta_P = obe.P.delta_P(sigma, k_idx); // optional<[n_delta, M, n_nu]>
    if (!delta_P.has_value()) { throw std::runtime_error("delta_P is not present in the projector. Cannot compute forces."); }
    auto n_delta = delta_P->shape()[0];

    // Dinv(n, ν) = 1/(ω_n + μ − ε_ν(k)).
    auto Dinv = nda::matrix<dcomplex>(n_omega, n_nu);
    for (auto n : range(n_omega))
      for (auto nu : R_nu) Dinv(n, nu) = 1.0 / (omegas[n] + mu - Hk(nu, nu));

    // Index tensor over band ν, batched over δ.  Q = obe.P[active rows], δQ = delta_P[active rows]:
    //   dQQ(δ, ν, a, b) = Q(a,ν)·conj(δQ(δ,b,ν)) + δQ(δ,a,ν)·conj(Q(b,ν))
    auto dQQ = nda::zeros<dcomplex>(n_delta, n_nu, rank, rank);
    for (auto delta_i : range(n_delta))
      for (auto a : range(rank))
        for (auto b : range(rank))
          for (auto nu : R_nu)
            dQQ(delta_i, nu, a, b) = P_obe(active[a], nu) * conj((*delta_P)(delta_i, active[b], nu))
               + (*delta_P)(delta_i, active[a], nu) * conj(P_obe(active[b], nu));

    // δG0_QQ(δ, n, a, b) = Σ_ν Dinv(n, ν) · dQQ(δ, ν, a, b) — one batched gemm per δ.
    auto dG0_QQ = nda::zeros<dcomplex>(n_delta, n_omega, rank, rank);
    if (rank > 0)
      for (auto delta_i : range(n_delta)) {
        auto dQQ_i = dQQ(delta_i, r_all, r_all, r_all);
        auto out_i = dG0_QQ(delta_i, r_all, r_all, r_all);
        nda::blas::gemm(1, Dinv, detail::as_2d(dQQ_i), 0, detail::as_2d(out_i));
      }
    return dG0_QQ;
  }

  //-------------------------------------------------------------------------------------------
  // Compute force contributions for a given k-point and spin (rank-reduced Woodbury).
  template <typename Mesh>
  nda::array<dcomplex, 2> trace_forces(one_body_elements_on_grid const &obe, double mu, long k_idx, long sigma,
                                       block2_gf<Mesh, matrix_valued> const &Sigma_dynamic,
                                       nda::array<nda::matrix<dcomplex>, 2> const &Sigma_static) {

    auto const &mesh = Sigma_dynamic(0, 0).mesh();
    auto omegas      = mesh | tl::to<std::vector<dcomplex>>();

    // delta_P must be present for the Pulay term.
    auto delta_P = obe.P.delta_P(sigma, k_idx);
    if (!delta_P.has_value()) { throw std::runtime_error("delta_P is not present in the projector. Cannot compute forces."); }
    auto n_delta = delta_P->shape()[0];

    auto result = nda::zeros<dcomplex>(n_delta, omegas.size());

    // Active subspace + block-diagonal Σ_active (Σ_total = Σ_dynamic + Σ_static), as in density().
    auto Sigma_total = detail::make_sigma_total(Sigma_dynamic, Sigma_static);
    auto decomp      = get_struct(Sigma_total).dims(r_all, 0) | tl::to<std::vector>();
    auto active      = detail::detect_active_subspace(Sigma_total, decomp);
    if (active.rank == 0) return result; // Σ ≡ 0 ⇒ no force correction.
    auto Sa_n = detail::compute_sigma_active(Sigma_total, active, sigma); // [n_w, rank, rank]

    // Bare G0_QQ (active block of P M₀⁻¹ P†) and the Pulay term δG0_QQ, both rank × rank.
    auto bare   = detail::compute_bare_projected(obe, obe.P, mu, k_idx, sigma, omegas, active);
    auto dG0_QQ = delta_G0_C_k_sigma(obe, mu, k_idx, sigma, omegas, active);

    // Per ω: K = Σ_a · (I − G0_QQ · Σ_a)⁻¹ (apply_K); force contribution = tr(K · δG0_QQ[δ]).
    for (long n = 0; n < long(omegas.size()); ++n) {
      auto Sa  = nda::matrix<dcomplex>{Sa_n(n, r_all, r_all)};
      auto Yaa = nda::matrix<dcomplex>{bare.G0_QQ(n, r_all, r_all)};
      for (auto delta_i : range(n_delta)) {
        auto K_dQQ         = detail::apply_K(Sa, Yaa, dG0_QQ(delta_i, n, r_all, r_all)); // K · δG0_QQ[δ]
        result(delta_i, n) = nda::trace(K_dQQ);
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
    force_gfs.reserve(n_delta);
    for (auto i : range(n_delta)) { force_gfs.push_back(gf{mesh}); }

    // ---------
    mpi::communicator comm = {}; // for now using default comm in MPI

#pragma omp parallel default(none) shared(n_k, comm, n_sigma, obe, mu, Sigma_dynamic, Sigma_static, n_delta, force_gfs, mesh)
    {
      // Thread-local gf objects
      std::vector<gf<Mesh, scalar_valued>> thread_force_gfs;
      thread_force_gfs.reserve(n_delta);
      for (auto i : range(n_delta)) { thread_force_gfs.push_back(gf{mesh}); }

#pragma omp for collapse(2)
      for (auto k_idx : mpi::chunk(range(n_k), comm)) {
        for (auto sigma : range(n_sigma)) {
          //  Correction term
          auto contrib = trace_forces(obe, mu, k_idx, sigma, Sigma_dynamic, Sigma_static); // [n_delta, n_omega]

          // Accumulate into gf objects with k-weights
          auto w_k = obe.H.k_weights(k_idx);
          for (auto delta_i : range(n_delta)) {
            for (auto n : range(mesh.size())) { thread_force_gfs[delta_i].data()(n) += w_k * contrib(delta_i, n); }
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
      result(delta_i)    = -real(density(force_gfs[delta_i])); // density(force_gfs) = dF/dtau, and forces = -dF/dtau.
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
