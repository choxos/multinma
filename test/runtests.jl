using multinma
using Test
using DataFrames
using CategoricalArrays
using Distributions
using Graphs

@testset "multinma.jl" begin

    # =========================================================================
    @testset "Priors" begin
        @testset "Normal prior" begin
            p = normal(; scale=10.0)
            @test p.dist == "Normal"
            @test p.location == 0.0
            @test p.scale == 10.0
            @test isnan(p.df)
        end

        @testset "Half-Normal prior" begin
            p = half_normal(; scale=5.0)
            @test p.dist == "half-Normal"
            @test p.scale == 5.0
        end

        @testset "Cauchy prior" begin
            p = cauchy(; location=1.0, scale=2.5)
            @test p.dist == "Cauchy"
            @test p.location == 1.0
            @test p.scale == 2.5
        end

        @testset "Half-Cauchy prior" begin
            p = half_cauchy(; scale=2.5)
            @test p.dist == "half-Cauchy"
            @test p.scale == 2.5
        end

        @testset "Student-t prior" begin
            p = student_t(; scale=2.5, df=3.0)
            @test p.dist == "Student t"
            @test p.df == 3.0
        end

        @testset "Half-Student-t prior" begin
            p = half_student_t(; scale=2.5, df=3.0)
            @test p.dist == "half-Student t"
        end

        @testset "Log-Normal prior" begin
            p = log_normal(; location=0.0, scale=1.0)
            @test p.dist == "log-Normal"
        end

        @testset "Log-Student-t prior" begin
            p = log_student_t(; location=0.0, scale=1.0, df=3.0)
            @test p.dist == "log-Student t"
        end

        @testset "Exponential prior" begin
            p = exponential_prior(; scale=1.0)
            @test p.dist == "Exponential"
            @test p.scale == 1.0

            p2 = exponential_prior(; rate=2.0)
            @test p2.scale == 0.5
        end

        @testset "Flat prior" begin
            p = flat()
            @test p.dist == "flat (implicit)"
        end

        @testset "Prior validation" begin
            @test_throws ArgumentError normal(; scale=-1.0)
            @test_throws ArgumentError half_normal(; scale=0.0)
            @test_throws ArgumentError normal(; location=Inf, scale=1.0)
        end

        @testset "Stan data translation" begin
            p = normal(; location=0.0, scale=10.0)
            sd = prior_to_stan_data(p, "prior_trt")
            @test sd["prior_trt_dist"] == 1
            @test sd["prior_trt_location"] == 0.0
            @test sd["prior_trt_scale"] == 10.0
        end

        @testset "Default marking" begin
            p = normal(; scale=10.0)
            @test !is_default(p)
            pd = set_default(p)
            @test is_default(pd)
        end
    end

    # =========================================================================
    @testset "Link functions" begin
        @test link_identity(5.0) == 5.0
        @test inv_identity(5.0) == 5.0

        @test link_log(exp(1.0)) ≈ 1.0
        @test inv_log(1.0) ≈ exp(1.0)

        @test link_logit(0.5) ≈ 0.0
        @test inv_logit(0.0) ≈ 0.5

        @test link_probit(0.5) ≈ 0.0 atol=1e-10
        @test inv_probit(0.0) ≈ 0.5

        @test link_cloglog(1 - exp(-1.0)) ≈ 0.0 atol=1e-10

        # Dispatch by name
        @test get_link("logit")(0.5) ≈ 0.0
        @test get_inv_link("log")(0.0) ≈ 1.0

        @test_throws ArgumentError get_link("invalid_link")
    end

    # =========================================================================
    @testset "Custom distributions" begin
        @testset "Generalised t" begin
            # Standard t with df=100 should be close to normal
            @test dgent(0.0, 100.0) ≈ pdf(Normal(), 0.0) atol=0.01
            @test pgent(0.0, 100.0) ≈ 0.5 atol=0.01
        end

        @testset "Log-t" begin
            @test dlogt(-1.0, 3.0) == 0.0
            @test plogt(0.0, 3.0) == 0.0
            @test qlogt(0.5, 3.0) > 0
        end

        @testset "Logit-Normal" begin
            @test dlogitnorm(0.0) == 0.0
            @test dlogitnorm(1.0) == 0.0
            @test dlogitnorm(0.5) > 0
            @test plogitnorm(0.5) ≈ 0.5 atol=0.01
            @test 0 < qlogitnorm(0.5) < 1
        end

        @testset "Bernoulli quantile" begin
            @test qbern(0.3, 0.5) == 0.0
            @test qbern(0.7, 0.5) == 1.0
        end
    end

    # =========================================================================
    @testset "Data setup" begin
        @testset "set_agd_arm (binary)" begin
            df = DataFrame(
                study = ["S1","S1","S2","S2"],
                trt = ["A","B","A","B"],
                r = [10, 20, 15, 25],
                n = [100, 100, 100, 100]
            )
            net = set_agd_arm(df; study=:study, trt=:trt, r=:r, n=:n, trt_ref="A")
            @test net isa NMAData
            @test has_agd_arm(net)
            @test !has_ipd(net)
            @test !has_agd_contrast(net)
            @test levels(net.treatments)[1] == "A"
            @test net.outcome[:agd_arm] == OUTCOME_BINARY
        end

        @testset "set_agd_arm (continuous)" begin
            df = DataFrame(
                study = ["S1","S1","S2","S2"],
                trt = ["A","B","A","B"],
                y = [1.0, 2.0, 1.5, 2.5],
                se = [0.5, 0.5, 0.5, 0.5]
            )
            net = set_agd_arm(df; study=:study, trt=:trt, y=:y, se=:se)
            @test net.outcome[:agd_arm] == OUTCOME_CONTINUOUS
        end

        @testset "set_agd_contrast" begin
            df = DataFrame(
                study = ["S1","S1","S2","S2"],
                trt = ["A","B","A","B"],
                y = [NaN, -0.5, NaN, -0.3],
                se = [NaN, 0.2, NaN, 0.15]
            )
            net = set_agd_contrast(df; study=:study, trt=:trt, y=:y, se=:se)
            @test has_agd_contrast(net)
            @test net.outcome[:agd_contrast] == OUTCOME_CONTINUOUS
        end

        @testset "set_ipd" begin
            df = DataFrame(
                study = repeat(["S1","S1"], 50),
                trt = repeat(["A","B"], 50),
                r = rand([0,1], 100)
            )
            net = set_ipd(df; study=:study, trt=:trt, r=:r)
            @test has_ipd(net)
            @test net.outcome[:ipd] == OUTCOME_BINARY
        end

        @testset "Outcome detection" begin
            @test detect_outcome_type(y=[1.0], se=[0.5]) == OUTCOME_CONTINUOUS
            @test detect_outcome_type(r=[1], n=[10]) == OUTCOME_BINARY
            @test detect_outcome_type(r=[1], E=[1.0]) == OUTCOME_RATE
            @test detect_outcome_type() == OUTCOME_NONE
        end
    end

    # =========================================================================
    @testset "Network operations" begin
        # Create a simple network
        df = DataFrame(
            study = ["S1","S1","S2","S2","S3","S3"],
            trt = ["A","B","B","C","A","C"],
            r = [10, 20, 15, 25, 12, 22],
            n = [100, 100, 100, 100, 100, 100]
        )
        net = set_agd_arm(df; study=:study, trt=:trt, r=:r, n=:n, trt_ref="A")

        @testset "Network graph" begin
            g, labels = to_graph(net)
            @test length(labels) == 3  # A, B, C
            @test Graphs.nv(g) == 3
            @test Graphs.ne(g) == 3  # A-B, B-C, A-C (triangle)
        end

        @testset "Connectivity" begin
            @test is_network_connected(net)
        end

        @testset "Direct/indirect evidence" begin
            @test has_direct(net, "A", "B")
            @test has_direct(net, "B", "C")
            @test has_direct(net, "A", "C")
            @test has_indirect(net, "A", "B")  # via C
        end

        @testset "Node-splits" begin
            splits = get_nodesplits(net)
            @test length(splits) >= 1
        end

        @testset "combine_network" begin
            df2 = DataFrame(
                study = ["S4","S4"],
                trt = ["A","D"],
                r = [10, 30],
                n = [100, 100]
            )
            net2 = set_agd_arm(df2; study=:study, trt=:trt, r=:r, n=:n, trt_ref="A")
            combined = combine_network(net, net2; trt_ref="A")
            @test length(levels(combined.treatments)) == 4  # A, B, C, D
        end
    end

    # =========================================================================
    @testset "NMA model (no Stan)" begin
        df = DataFrame(
            study = ["S1","S1","S2","S2","S3","S3"],
            trt = ["A","B","B","C","A","C"],
            r = [10, 20, 15, 25, 12, 22],
            n = [100, 100, 100, 100, 100, 100]
        )
        net = set_agd_arm(df; study=:study, trt=:trt, r=:r, n=:n, trt_ref="A")

        # nma() should succeed but with stanfit=nothing
        fit = nma(net; trt_effects=:fixed)
        @test fit isa StanNMA
        @test fit.trt_effects == :fixed
        @test fit.consistency == :consistency
        @test fit.likelihood == "binomial_1par"
        @test fit.link == "logit"
        @test isnothing(fit.stanfit)

        # Random effects
        fit_re = nma(net; trt_effects=:random)
        @test fit_re.trt_effects == :random

        # Invalid args
        @test_throws ArgumentError nma(net; trt_effects=:invalid)
        @test_throws ArgumentError nma(net; consistency=:invalid)
    end

    # =========================================================================
    @testset "Stan interface" begin
        @testset "Model directory" begin
            dir = stan_model_dir()
            @test isdir(dir)
        end

        @testset "Available models" begin
            models = available_stan_models()
            @test "binomial_1par" in models || "normal" in models || length(models) >= 0
        end
    end

    # =========================================================================
    @testset "Survival distributions" begin
        @testset "Exponential" begin
            d = ExponentialSurv(0.5)
            @test surv_survival(d, 0.0) == 1.0
            @test surv_survival(d, Inf) == 0.0
            @test surv_hazard(d, 1.0) == 0.5
            @test surv_cumhaz(d, 2.0) ≈ 1.0
            @test surv_quantile(d, 0.5) ≈ log(2) / 0.5
            @test surv_pdf(d, 1.0) ≈ 0.5 * exp(-0.5)
            @test surv_cdf(d, 1.0) ≈ 1 - exp(-0.5)
        end

        @testset "Weibull" begin
            d = WeibullSurv(2.0, 1.0)
            @test surv_survival(d, 0.0) == 1.0
            @test surv_hazard(d, 1.0) ≈ 2.0
            @test surv_cumhaz(d, 1.0) ≈ 1.0
        end

        @testset "Gompertz" begin
            d = GompertzSurv(0.1, 0.05)
            @test surv_survival(d, 0.0) == 1.0
            @test surv_hazard(d, 0.0) ≈ 0.1
        end

        @testset "Log-Normal" begin
            d = LogNormalSurv(0.0, 1.0)
            @test surv_survival(d, 0.0) == 1.0
            @test surv_quantile(d, 0.5) ≈ 1.0
        end

        @testset "Log-Logistic" begin
            d = LogLogisticSurv(2.0, 1.0)
            @test surv_survival(d, 0.0) == 1.0
            @test surv_survival(d, 1.0) ≈ 0.5
        end

        @testset "Gamma" begin
            d = GammaSurv(2.0, 1.0)
            @test surv_survival(d, 0.0) == 1.0
            @test surv_cdf(d, 0.0) == 0.0
        end

        @testset "RMST" begin
            d = ExponentialSurv(1.0)
            r = rmst(d, 10.0)
            @test r > 0
            @test r ≈ 1 - exp(-10.0) atol=0.01
        end

        @testset "SurvObs" begin
            obs = surv([1.0, 2.0, 3.0], [1, 0, 1])
            @test length(obs) == 3
            @test obs[1].event == true
            @test obs[2].event == false
        end

        @testset "Kaplan-Meier" begin
            obs = surv([1.0, 2.0, 3.0, 1.5, 2.5], [1, 1, 1, 0, 1])
            km = kaplan_meier(obs)
            @test km.survival[1] == 1.0  # time 0
            @test all(diff(km.survival) .<= 0)  # non-increasing
        end
    end

    # =========================================================================
    @testset "M-splines" begin
        @testset "make_knots" begin
            ki = make_knots([1.0, 2.0, 3.0, 4.0, 5.0]; df=5)
            @test ki.df == 5
            @test ki.boundary_knots[1] == 0.0
            @test ki.boundary_knots[2] > 5.0
        end

        @testset "M-spline basis" begin
            ki = make_knots([1.0, 2.0, 3.0, 4.0, 5.0]; df=3)
            basis = mspline_basis(2.5, ki)
            @test length(basis) == 3
            @test all(basis .>= 0)  # M-splines are non-negative
        end

        @testset "MSplineHazard" begin
            ki = make_knots([1.0, 2.0, 3.0, 4.0, 5.0]; df=3)
            model = MSplineHazard([1.0, 1.0, 1.0], ki)
            @test mspline_hazard(model, 2.5) >= 0
            @test mspline_survival(model, 0.0) ≈ 1.0 atol=0.01
        end
    end

    # =========================================================================
    @testset "MCMCArray" begin
        data = randn(100, 4, 3)
        params = ["d[1]", "d[2]", "d[3]"]
        mcmc = MCMCArray(data, params)
        @test n_iterations(mcmc) == 100
        @test n_chains(mcmc) == 4
        @test n_parameters(mcmc) == 3
        @test size(mcmc) == (100, 4, 3)

        @test_throws DimensionMismatch MCMCArray(data, ["a", "b"])
    end

    # =========================================================================
    @testset "Summary functions" begin
        data = randn(100, 2, 3)
        params = ["p1", "p2", "p3"]
        mcmc = MCMCArray(data, params)
        df = summarise_draws(mcmc)
        @test nrow(df) == 3
        @test hasproperty(df, :mean)
        @test hasproperty(df, :sd)
        @test hasproperty(df, Symbol("50%"))
    end

    # =========================================================================
    @testset "Ranking helper" begin
        ranks = multinma._rank([3.0, 1.0, 2.0])
        @test ranks == [3.0, 1.0, 2.0]

        # Ties
        ranks_tied = multinma._rank([1.0, 2.0, 2.0])
        @test ranks_tied[1] == 1.0
        @test ranks_tied[2] == 2.5
        @test ranks_tied[3] == 2.5
    end

    # =========================================================================
    @testset "Datasets" begin
        @testset "example_smoking" begin
            df = example_smoking()
            @test df isa DataFrame
            @test nrow(df) == 49
            @test hasproperty(df, :trt)
            @test hasproperty(df, :r)
            @test hasproperty(df, :n)
        end

        @testset "example_blocker" begin
            df = example_blocker()
            @test df isa DataFrame
            @test nrow(df) == 10
        end

        @testset "Smoking NMA setup" begin
            smoking = example_smoking()
            net = set_agd_arm(smoking; study=:studyn, trt=:trt,
                              r=:r, n=:n, trt_ref="No contact")
            @test is_network_connected(net)
            @test length(levels(net.treatments)) == 4

            fit = nma(net; trt_effects=:fixed)
            @test fit.likelihood == "binomial_1par"
            @test fit.link == "logit"
        end
    end

    # =========================================================================
    @testset "Utilities" begin
        @testset "softmax" begin
            x = [1.0, 2.0, 3.0]
            s = multinma.softmax(x)
            @test sum(s) ≈ 1.0
            @test all(s .> 0)
        end

        @testset "make_treatment_factor" begin
            f = make_treatment_factor(["B", "A", "C"]; trt_ref="A")
            @test levels(f)[1] == "A"
        end

        @testset "make_contrasts" begin
            c = make_contrasts(["A", "B", "C"])
            @test nrow(c) == 3  # AB, AC, BC
        end
    end

    # =========================================================================
    @testset "Display methods" begin
        # Just test they don't error
        p = normal(; scale=10.0)
        io = IOBuffer()
        show(io, p)
        @test length(take!(io)) > 0

        df = DataFrame(
            study = ["S1","S1"],
            trt = ["A","B"],
            r = [10, 20],
            n = [100, 100]
        )
        net = set_agd_arm(df; study=:study, trt=:trt, r=:r, n=:n)
        show(io, net)
        @test length(take!(io)) > 0

        fit = nma(net; trt_effects=:fixed)
        show(io, fit)
        @test length(take!(io)) > 0
    end
end
