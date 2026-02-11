# =============================================================================
# Vignette Comparison Tests
# =============================================================================
# Mirrors ALL R multinma vignette examples and compares data pipeline results
# against the R package. MCMC-dependent results will be testable once CmdStan
# is connected.
# =============================================================================

using multinma
using Test
using DataFrames
using CategoricalArrays
using Distributions
using Graphs
using Statistics

@testset "Vignette Comparison Tests" begin

# =========================================================================
# 1. Smoking Cessation (Hasselblad 1998) - TSD 4
# =========================================================================
@testset "Smoking Cessation" begin
    smoking = load_dataset("smoking")

    @testset "Data loading" begin
        @test nrow(smoking) == 50
        @test ncol(smoking) == 5
        @test "studyn" in names(smoking)
        @test "trtc" in names(smoking)
        @test "r" in names(smoking)
        @test "n" in names(smoking)
    end

    @testset "Network setup (set_agd_arm)" begin
        smknet = set_agd_arm(smoking;
                             study = :studyn,
                             trt = :trtc,
                             r = :r,
                             n = :n,
                             trt_ref = "No intervention")

        # R: 4 treatments, 24 studies, ref = "No intervention"
        @test length(levels(smknet.treatments)) == 4
        @test length(levels(smknet.studies)) == 24
        @test levels(smknet.treatments)[1] == "No intervention"
        @test !isnothing(smknet.agd_arm)
        @test isnothing(smknet.agd_contrast)
        @test isnothing(smknet.ipd)

        trts = levels(smknet.treatments)
        @test "Group counselling" in trts
        @test "Individual counselling" in trts
        @test "Self-help" in trts

        # Binary outcome auto-detected
        @test smknet.outcome[:agd_arm] == multinma.OUTCOME_BINARY
    end

    @testset "Network connectivity" begin
        smknet = set_agd_arm(smoking;
                             study = :studyn, trt = :trtc, r = :r, n = :n,
                             trt_ref = "No intervention")

        @test is_network_connected(smknet)
        g, trt_labels = to_graph(smknet)
        @test Graphs.nv(g) == 4
        @test Graphs.ne(g) >= 3
    end

    @testset "Node-splitting" begin
        smknet = set_agd_arm(smoking;
                             study = :studyn, trt = :trtc, r = :r, n = :n,
                             trt_ref = "No intervention")

        ns = get_nodesplits(smknet)
        # R: 2+ node-splits (we may find more valid splits than R's default)
        @test length(ns) >= 2
    end

    @testset "Prior construction" begin
        p_trt = normal(; scale = 100.0)
        @test p_trt.dist == "Normal"
        @test p_trt.location == 0.0
        @test p_trt.scale == 100.0

        p_het = half_normal(; scale = 5.0)
        @test p_het.dist == "half-Normal"
        @test p_het.scale == 5.0
    end

    @testset "Model configuration (RE consistency)" begin
        smknet = set_agd_arm(smoking;
                             study = :studyn, trt = :trtc, r = :r, n = :n,
                             trt_ref = "No intervention")

        smkfit = nma(smknet;
                     trt_effects = :random,
                     prior_intercept = normal(; scale = 100.0),
                     prior_trt = normal(; scale = 100.0),
                     prior_het = half_normal(; scale = 5.0))

        @test smkfit.likelihood == "binomial_1par"
        @test smkfit.link == "logit"
        @test smkfit.trt_effects == :random
        @test smkfit.consistency == :consistency
        @test isnothing(smkfit.stanfit)
    end

    @testset "Model configuration (UME)" begin
        smknet = set_agd_arm(smoking;
                             study = :studyn, trt = :trtc, r = :r, n = :n,
                             trt_ref = "No intervention")

        smkfit_ume = nma(smknet;
                         consistency = :ume,
                         trt_effects = :random,
                         prior_intercept = normal(; scale = 100.0),
                         prior_trt = normal(; scale = 100.0),
                         prior_het = half_normal(; scale = 5.0))

        @test smkfit_ume.consistency == :ume
        @test smkfit_ume.trt_effects == :random
    end
end

# =========================================================================
# 2. Beta Blockers (Carlin 1992) - TSD 2
# =========================================================================
@testset "Beta Blockers" begin
    blocker = load_dataset("blocker")

    @testset "Data loading" begin
        @test nrow(blocker) == 44
        @test ncol(blocker) == 5
    end

    @testset "Network setup" begin
        blocker_net = set_agd_arm(blocker;
                                  study = :studyn, trt = :trtc,
                                  r = :r, n = :n,
                                  trt_ref = "Control")

        @test length(levels(blocker_net.treatments)) == 2
        @test length(levels(blocker_net.studies)) == 22
        @test levels(blocker_net.treatments)[1] == "Control"
        @test "Beta Blocker" in levels(blocker_net.treatments)
        @test blocker_net.outcome[:agd_arm] == multinma.OUTCOME_BINARY
    end

    @testset "Model configuration (FE)" begin
        blocker_net = set_agd_arm(blocker;
                                  study = :studyn, trt = :trtc,
                                  r = :r, n = :n,
                                  trt_ref = "Control")

        fit_FE = nma(blocker_net;
                     trt_effects = :fixed,
                     prior_intercept = normal(; scale = 100.0),
                     prior_trt = normal(; scale = 100.0))

        @test fit_FE.likelihood == "binomial_1par"
        @test fit_FE.link == "logit"
        @test fit_FE.trt_effects == :fixed
    end

    @testset "Model configuration (RE)" begin
        blocker_net = set_agd_arm(blocker;
                                  study = :studyn, trt = :trtc,
                                  r = :r, n = :n,
                                  trt_ref = "Control")

        fit_RE = nma(blocker_net;
                     trt_effects = :random,
                     prior_intercept = normal(; scale = 100.0),
                     prior_trt = normal(; scale = 100.0),
                     prior_het = half_normal(; scale = 5.0))

        @test fit_RE.trt_effects == :random
        @test fit_RE.priors[:het].dist == "half-Normal"
    end

    @testset "Network is pairwise" begin
        blocker_net = set_agd_arm(blocker;
                                  study = :studyn, trt = :trtc,
                                  r = :r, n = :n,
                                  trt_ref = "Control")

        g, _ = to_graph(blocker_net)
        @test Graphs.nv(g) == 2
        @test Graphs.ne(g) == 1
        ns = get_nodesplits(blocker_net)
        @test length(ns) == 0
    end
end

# =========================================================================
# 3. Thrombolytics (Boland 2003) - TSD 4
# =========================================================================
@testset "Thrombolytics" begin
    thrombo = load_dataset("thrombolytics")

    @testset "Data loading" begin
        @test nrow(thrombo) == 102
        @test ncol(thrombo) == 5
    end

    @testset "Network setup" begin
        thrombo_net = set_agd_arm(thrombo;
                                  study = :studyn, trt = :trtc,
                                  r = :r, n = :n,
                                  trt_ref = "SK")

        @test length(levels(thrombo_net.treatments)) == 9
        @test length(levels(thrombo_net.studies)) == 50
        @test levels(thrombo_net.treatments)[1] == "SK"
        @test thrombo_net.outcome[:agd_arm] == multinma.OUTCOME_BINARY

        trts = levels(thrombo_net.treatments)
        for t in ["SK", "t-PA", "Acc t-PA", "SK + t-PA", "r-PA",
                   "TNK", "PTCA", "UK", "ASPAC"]
            @test t in trts
        end
    end

    @testset "Network connectivity" begin
        thrombo_net = set_agd_arm(thrombo;
                                  study = :studyn, trt = :trtc,
                                  r = :r, n = :n,
                                  trt_ref = "SK")

        @test is_network_connected(thrombo_net)
        g, _ = to_graph(thrombo_net)
        @test Graphs.nv(g) == 9
    end

    @testset "Node-splitting" begin
        thrombo_net = set_agd_arm(thrombo;
                                  study = :studyn, trt = :trtc,
                                  r = :r, n = :n,
                                  trt_ref = "SK")

        ns = get_nodesplits(thrombo_net)
        @test length(ns) >= 2

        @test has_direct(thrombo_net, "SK", "t-PA")
        @test has_direct(thrombo_net, "SK", "Acc t-PA")
        @test has_direct(thrombo_net, "Acc t-PA", "ASPAC")
    end

    @testset "FE model" begin
        thrombo_net = set_agd_arm(thrombo;
                                  study = :studyn, trt = :trtc,
                                  r = :r, n = :n,
                                  trt_ref = "SK")

        thrombo_fit = nma(thrombo_net;
                          trt_effects = :fixed,
                          prior_intercept = normal(; scale = 100.0),
                          prior_trt = normal(; scale = 100.0))

        @test thrombo_fit.likelihood == "binomial_1par"
        @test thrombo_fit.link == "logit"
        @test thrombo_fit.trt_effects == :fixed
    end

    @testset "UME model" begin
        thrombo_net = set_agd_arm(thrombo;
                                  study = :studyn, trt = :trtc,
                                  r = :r, n = :n,
                                  trt_ref = "SK")

        thrombo_ume = nma(thrombo_net;
                          consistency = :ume,
                          trt_effects = :fixed,
                          prior_intercept = normal(; scale = 100.0),
                          prior_trt = normal(; scale = 100.0))

        @test thrombo_ume.consistency == :ume
    end
end

# =========================================================================
# 4. Parkinson's Disease (TSD 2)
# =========================================================================
@testset "Parkinson's Disease" begin
    parkinsons = load_dataset("parkinsons")

    @testset "Data loading" begin
        @test nrow(parkinsons) == 15
        @test ncol(parkinsons) == 7
        @test "y" in names(parkinsons)
        @test "se" in names(parkinsons)
        @test "diff" in names(parkinsons)
        @test "se_diff" in names(parkinsons)
    end

    @testset "Arm-based network" begin
        arm_net = set_agd_arm(parkinsons;
                              study = :studyn, trt = :trtn,
                              y = :y, se = :se, sample_size = :n,
                              trt_ref = "4")

        @test length(levels(arm_net.treatments)) == 5
        @test length(levels(arm_net.studies)) == 7
        @test levels(arm_net.treatments)[1] == "4"
        @test arm_net.outcome[:agd_arm] == multinma.OUTCOME_CONTINUOUS
    end

    @testset "Arm-based FE model" begin
        arm_net = set_agd_arm(parkinsons;
                              study = :studyn, trt = :trtn,
                              y = :y, se = :se, sample_size = :n,
                              trt_ref = "4")

        arm_fit_FE = nma(arm_net;
                         trt_effects = :fixed,
                         prior_intercept = normal(; scale = 100.0),
                         prior_trt = normal(; scale = 10.0))

        @test arm_fit_FE.likelihood == "normal"
        @test arm_fit_FE.link == "identity"
        @test arm_fit_FE.trt_effects == :fixed
    end

    @testset "Arm-based RE model" begin
        arm_net = set_agd_arm(parkinsons;
                              study = :studyn, trt = :trtn,
                              y = :y, se = :se, sample_size = :n,
                              trt_ref = "4")

        arm_fit_RE = nma(arm_net;
                         trt_effects = :random,
                         prior_intercept = normal(; scale = 100.0),
                         prior_trt = normal(; scale = 100.0),
                         prior_het = half_normal(; scale = 5.0))

        @test arm_fit_RE.trt_effects == :random
        @test arm_fit_RE.likelihood == "normal"
    end

    @testset "Contrast-based network" begin
        contr_net = set_agd_contrast(parkinsons;
                                     study = :studyn, trt = :trtn,
                                     y = :diff, se = :se_diff,
                                     sample_size = :n,
                                     trt_ref = "4")

        @test length(levels(contr_net.treatments)) == 5
        @test length(levels(contr_net.studies)) == 7
        @test levels(contr_net.treatments)[1] == "4"
        @test contr_net.outcome[:agd_contrast] == multinma.OUTCOME_CONTINUOUS
        @test isnothing(contr_net.agd_arm)
        @test !isnothing(contr_net.agd_contrast)
    end

    @testset "Contrast-based FE model" begin
        contr_net = set_agd_contrast(parkinsons;
                                     study = :studyn, trt = :trtn,
                                     y = :diff, se = :se_diff,
                                     sample_size = :n,
                                     trt_ref = "4")

        contr_fit_FE = nma(contr_net;
                           trt_effects = :fixed,
                           prior_trt = normal(; scale = 100.0))

        @test contr_fit_FE.likelihood == "normal"
        @test contr_fit_FE.link == "identity"
    end

    @testset "Mixed arm+contrast network" begin
        parkinsons_arm = filter(r -> r.studyn in [1, 2, 3], parkinsons)
        parkinsons_contr = filter(r -> r.studyn in [4, 5, 6, 7], parkinsons)

        mix_arm_net = set_agd_arm(parkinsons_arm;
                                   study = :studyn, trt = :trtn,
                                   y = :y, se = :se, sample_size = :n,
                                   trt_ref = "4")

        mix_contr_net = set_agd_contrast(parkinsons_contr;
                                         study = :studyn, trt = :trtn,
                                         y = :diff, se = :se_diff,
                                         sample_size = :n,
                                         trt_ref = "4")

        mix_net = combine_network(mix_arm_net, mix_contr_net)

        @test length(levels(mix_net.treatments)) == 5
        @test length(levels(mix_net.studies)) == 7
        @test !isnothing(mix_net.agd_arm)
        @test !isnothing(mix_net.agd_contrast)
    end

    @testset "Mixed FE model" begin
        parkinsons_arm = filter(r -> r.studyn in [1, 2, 3], parkinsons)
        parkinsons_contr = filter(r -> r.studyn in [4, 5, 6, 7], parkinsons)

        mix_arm_net = set_agd_arm(parkinsons_arm;
                                   study = :studyn, trt = :trtn,
                                   y = :y, se = :se, sample_size = :n,
                                   trt_ref = "4")

        mix_contr_net = set_agd_contrast(parkinsons_contr;
                                         study = :studyn, trt = :trtn,
                                         y = :diff, se = :se_diff,
                                         sample_size = :n,
                                         trt_ref = "4")

        mix_net = combine_network(mix_arm_net, mix_contr_net)

        mix_fit_FE = nma(mix_net;
                         trt_effects = :fixed,
                         prior_intercept = normal(; scale = 100.0),
                         prior_trt = normal(; scale = 100.0))

        @test mix_fit_FE.likelihood == "normal"
        @test mix_fit_FE.link == "identity"
        @test mix_fit_FE.trt_effects == :fixed
    end
end

# =========================================================================
# 5. Diabetes (Elliott 2007) - TSD 2
# =========================================================================
@testset "Diabetes" begin
    diabetes = load_dataset("diabetes")

    @testset "Data loading" begin
        @test nrow(diabetes) == 48
        @test ncol(diabetes) == 7
        @test "studyc" in names(diabetes)
        @test "trtc" in names(diabetes)
        @test "r" in names(diabetes)
        @test "n" in names(diabetes)
        @test "time" in names(diabetes)
    end

    @testset "Network setup" begin
        db_net = set_agd_arm(diabetes;
                             study = :studyc, trt = :trtc,
                             r = :r, n = :n,
                             trt_ref = "Placebo")

        @test length(levels(db_net.treatments)) == 6
        @test length(levels(db_net.studies)) == 22
        @test levels(db_net.treatments)[1] == "Placebo"
        @test db_net.outcome[:agd_arm] == multinma.OUTCOME_BINARY

        for t in ["ACE Inhibitor", "ARB", "Beta Blocker", "CCB",
                   "Diuretic", "Placebo"]
            @test t in levels(db_net.treatments)
        end
    end

    @testset "Network connectivity" begin
        db_net = set_agd_arm(diabetes;
                             study = :studyc, trt = :trtc,
                             r = :r, n = :n,
                             trt_ref = "Placebo")
        @test is_network_connected(db_net)
        g, _ = to_graph(db_net)
        @test Graphs.nv(g) == 6
    end

    @testset "FE model with cloglog link" begin
        db_net = set_agd_arm(diabetes;
                             study = :studyc, trt = :trtc,
                             r = :r, n = :n,
                             trt_ref = "Placebo")

        db_fit_FE = nma(db_net;
                        trt_effects = :fixed,
                        link = "cloglog",
                        prior_intercept = normal(; scale = 100.0),
                        prior_trt = normal(; scale = 100.0))

        @test db_fit_FE.link == "cloglog"
        @test db_fit_FE.likelihood == "binomial_1par"
        @test db_fit_FE.trt_effects == :fixed
    end

    @testset "RE model with cloglog link" begin
        db_net = set_agd_arm(diabetes;
                             study = :studyc, trt = :trtc,
                             r = :r, n = :n,
                             trt_ref = "Placebo")

        db_fit_RE = nma(db_net;
                        trt_effects = :random,
                        link = "cloglog",
                        prior_intercept = normal(; scale = 10.0),
                        prior_trt = normal(; scale = 10.0),
                        prior_het = half_normal(; scale = 5.0))

        @test db_fit_RE.trt_effects == :random
        @test db_fit_RE.link == "cloglog"
    end
end

# =========================================================================
# 6. Statins
# =========================================================================
@testset "Statins" begin
    statins = load_dataset("statins")

    @testset "Data loading" begin
        @test nrow(statins) == 38
        @test ncol(statins) == 7
    end

    @testset "Network setup" begin
        stat_net = set_agd_arm(statins;
                               study = :studyc, trt = :trtc,
                               r = :r, n = :n)

        @test length(levels(stat_net.treatments)) >= 2
        @test is_network_connected(stat_net)
        @test stat_net.outcome[:agd_arm] == multinma.OUTCOME_BINARY
        g, _ = to_graph(stat_net)
        @test Graphs.nv(g) >= 2
    end
end

# =========================================================================
# 7. Transfusion
# =========================================================================
@testset "Transfusion" begin
    transfusion = load_dataset("transfusion")

    @testset "Data loading" begin
        @test nrow(transfusion) == 12
        @test ncol(transfusion) == 4
    end

    @testset "Network setup" begin
        trans_net = set_agd_arm(transfusion;
                                study = :studyc, trt = :trtc,
                                r = :r, n = :n)

        @test length(levels(trans_net.treatments)) >= 2
        @test is_network_connected(trans_net)
        @test trans_net.outcome[:agd_arm] == multinma.OUTCOME_BINARY
    end
end

# =========================================================================
# 8. Dietary Fat
# =========================================================================
@testset "Dietary Fat" begin
    dietary = load_dataset("dietary_fat")

    @testset "Data loading" begin
        @test nrow(dietary) == 21
        @test ncol(dietary) == 7
    end

    @testset "Network setup" begin
        diet_net = set_agd_arm(dietary;
                               study = :studyc, trt = :trtc,
                               r = :r, n = :n)

        @test length(levels(diet_net.treatments)) >= 2
        @test is_network_connected(diet_net)
        @test diet_net.outcome[:agd_arm] == multinma.OUTCOME_BINARY
    end
end

# =========================================================================
# 9. Atrial Fibrillation
# =========================================================================
@testset "Atrial Fibrillation" begin
    af = load_dataset("atrial_fibrillation")

    @testset "Data loading" begin
        @test nrow(af) == 63
        @test ncol(af) == 11
    end
end

# =========================================================================
# 10. All Prior Types
# =========================================================================
@testset "All Prior Types" begin
    @test normal(; scale = 100.0).dist == "Normal"
    @test half_normal(; scale = 5.0).dist == "half-Normal"
    @test cauchy(; scale = 2.5).dist == "Cauchy"
    @test half_cauchy(; scale = 2.5).dist == "half-Cauchy"
    @test student_t(; df = 3.0, scale = 2.5).dist == "Student t"
    @test student_t(; df = 3.0, scale = 2.5).df == 3.0
    @test half_student_t(; df = 3.0, scale = 2.5).dist == "half-Student t"
    @test log_normal(; location = 0.0, scale = 1.0).dist == "log-Normal"
    @test exponential_prior(; scale = 1.0).dist == "Exponential"
    @test flat().dist == "flat (implicit)"
end

# =========================================================================
# 11. Link Functions
# =========================================================================
@testset "Link Functions" begin
    for link_name in ["logit", "probit", "cloglog", "log", "identity"]
        f = get_link(link_name)
        inv_f = get_inv_link(link_name)
        @test f isa Function
        @test inv_f isa Function
    end

    @test get_link("logit")(0.5) ≈ 0.0 atol = 1e-10
    @test get_inv_link("logit")(0.0) ≈ 0.5 atol = 1e-10
    @test get_link("identity")(42.0) == 42.0
    @test get_inv_link("identity")(42.0) == 42.0
    @test get_link("log")(1.0) ≈ 0.0 atol = 1e-10
    @test get_inv_link("log")(0.0) ≈ 1.0 atol = 1e-10
    @test get_link("probit")(0.5) ≈ 0.0 atol = 1e-10
    @test get_inv_link("probit")(0.0) ≈ 0.5 atol = 1e-10

    cll = get_link("cloglog")
    icll = get_inv_link("cloglog")
    @test icll(cll(0.3)) ≈ 0.3 atol = 1e-10
end

# =========================================================================
# 12. Survival Distributions
# =========================================================================
@testset "Survival Distributions" begin
    exp_dist = ExponentialSurv(0.5)
    @test surv_hazard(exp_dist, 1.0) ≈ 0.5
    @test surv_survival(exp_dist, 0.0) ≈ 1.0

    wb = WeibullSurv(1.0, 0.5)
    @test surv_survival(wb, 0.0) ≈ 1.0
    @test surv_hazard(wb, 1.0) > 0

    gp = GompertzSurv(0.1, 0.5)
    @test surv_survival(gp, 0.0) ≈ 1.0

    ln = LogNormalSurv(0.0, 1.0)
    @test surv_survival(ln, 0.0) ≈ 1.0

    ll = LogLogisticSurv(1.0, 1.0)
    @test surv_survival(ll, 0.0) ≈ 1.0

    ga = GammaSurv(1.0, 1.0)
    @test surv_survival(ga, 0.0) ≈ 1.0

    gg = GenGammaSurv(0.0, 1.0, 1.0)
    @test surv_survival(gg, 0.0) ≈ 1.0
end

# =========================================================================
# 13. M-spline Basis
# =========================================================================
@testset "M-spline Basis" begin
    knots = make_knots([1.0, 2.0, 3.0, 5.0, 8.0, 10.0]; df = 5)
    @test length(knots) > 0

    bs = mspline_basis(5.0, knots; order = 3)
    @test length(bs) > 0
    @test all(bs .>= 0)
end

# =========================================================================
# 14. MCMCArray
# =========================================================================
@testset "MCMCArray" begin
    data = randn(100, 4, 3)
    pnames = ["d[1]", "d[2]", "tau"]
    mc = MCMCArray(data, pnames)

    @test size(mc.data) == (100, 4, 3)
    @test mc.parameters == pnames
end

# =========================================================================
# 15. Summary with MCMCArray
# =========================================================================
@testset "Summary functions" begin
    data = randn(1000, 4, 3)
    pnames = ["d[1]", "d[2]", "tau"]
    mc = MCMCArray(data, pnames)

    summ = summarise_draws(mc)
    @test nrow(summ) == 3
    @test "parameter" in names(summ)
    @test "mean" in names(summ)
    @test "sd" in names(summ)
    @test "2%" in names(summ) || "3%" in names(summ)
end

# =========================================================================
# 16. Custom Distributions (function-based API)
# =========================================================================
@testset "Custom distributions" begin
    # Logit-Normal
    @test 0.0 < qlogitnorm(0.5) < 1.0
    @test qlogitnorm(0.5; location = 0.0, scale = 1.0) ≈ 0.5 atol = 1e-10
    @test dlogitnorm(0.5) > 0
    @test plogitnorm(0.5) ≈ 0.5 atol = 1e-10

    # Generalised Student's t
    @test dgent(0.0, 3.0) > 0
    @test pgent(0.0, 3.0) ≈ 0.5 atol = 1e-10
    @test qgent(0.5, 3.0) ≈ 0.0 atol = 1e-10

    # Log Student's t
    @test dlogt(1.0, 3.0) > 0
    @test plogt(1.0, 3.0) ≈ 0.5 atol = 1e-10
    @test qlogt(0.5, 3.0) ≈ 1.0 atol = 1e-10

    # Bernoulli quantile
    @test qbern(0.3, 0.5) == 0.0
    @test qbern(0.7, 0.5) == 1.0
end

# =========================================================================
# 17. Plotting data generation
# =========================================================================
@testset "Plotting data" begin
    smoking = load_dataset("smoking")
    smknet = set_agd_arm(smoking;
                         study = :studyn, trt = :trtc, r = :r, n = :n,
                         trt_ref = "No intervention")

    plot_data = plot_network(smknet)
    @test haskey(plot_data, :graph)
    @test haskey(plot_data, :labels)
    @test length(plot_data[:labels]) == 4
end

# =========================================================================
# 18. Stan interface
# =========================================================================
@testset "Stan interface" begin
    @test isdir(stan_model_dir())
    models = available_stan_models()
    @test length(models) > 0
end

# =========================================================================
# 19. All vignette model configs
# =========================================================================
@testset "All vignette model configs" begin
    smoking = load_dataset("smoking")
    blocker = load_dataset("blocker")
    thrombo = load_dataset("thrombolytics")
    parkinsons = load_dataset("parkinsons")
    diabetes = load_dataset("diabetes")

    # Smoking: RE, binomial, logit
    smknet = set_agd_arm(smoking; study=:studyn, trt=:trtc, r=:r, n=:n,
                         trt_ref="No intervention")
    smkfit = nma(smknet; trt_effects=:random,
                 prior_intercept=normal(;scale=100.0),
                 prior_trt=normal(;scale=100.0),
                 prior_het=half_normal(;scale=5.0))
    @test smkfit.likelihood == "binomial_1par"
    @test smkfit.link == "logit"

    # Blocker FE: binomial, logit
    bnet = set_agd_arm(blocker; study=:studyn, trt=:trtc, r=:r, n=:n,
                       trt_ref="Control")
    bfit = nma(bnet; trt_effects=:fixed,
               prior_intercept=normal(;scale=100.0),
               prior_trt=normal(;scale=100.0))
    @test bfit.likelihood == "binomial_1par"
    @test bfit.link == "logit"

    # Thrombo FE: binomial, logit
    tnet = set_agd_arm(thrombo; study=:studyn, trt=:trtc, r=:r, n=:n)
    tfit = nma(tnet; trt_effects=:fixed,
               prior_intercept=normal(;scale=100.0),
               prior_trt=normal(;scale=100.0))
    @test tfit.likelihood == "binomial_1par"
    @test tfit.link == "logit"

    # Parkinsons FE: normal, identity
    anet = set_agd_arm(parkinsons; study=:studyn, trt=:trtn, y=:y, se=:se,
                       sample_size=:n)
    afit = nma(anet; trt_effects=:fixed,
               prior_intercept=normal(;scale=100.0),
               prior_trt=normal(;scale=10.0))
    @test afit.likelihood == "normal"
    @test afit.link == "identity"

    # Diabetes FE: binomial, cloglog
    dnet = set_agd_arm(diabetes; study=:studyc, trt=:trtc, r=:r, n=:n)
    dfit = nma(dnet; trt_effects=:fixed, link="cloglog",
               prior_intercept=normal(;scale=100.0),
               prior_trt=normal(;scale=100.0))
    @test dfit.likelihood == "binomial_1par"
    @test dfit.link == "cloglog"
end

end  # outer testset
