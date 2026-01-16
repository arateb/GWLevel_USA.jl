using Test
using GWLevel_USA
using DataFrames
using Statistics

@testset "GWLevel_USA.jl" begin

    @testset "Config" begin
        # Default config
        cfg = GWConfig()
        @test cfg.baseline_start == 1940
        @test cfg.baseline_end == 1955
        @test cfg.idw_power == 2.0

        # Custom config
        cfg2 = GWConfig(baseline_start=1950, mc_samples=5000)
        @test cfg2.baseline_start == 1950
        @test cfg2.mc_samples == 5000

        # Sy priors
        sy = get_sy_prior("C_S_High_Plains")
        @test sy.mean == 0.15
        @test sy.sd == 0.04
        @test sy.confinement == :unconfined

        # Unknown aquifer returns default
        sy_unknown = get_sy_prior("Unknown_Aquifer")
        @test sy_unknown.mean == 0.10
    end

    @testset "Spatial - Projection" begin
        # Test projection at known point (Austin, TX)
        lon, lat = -97.7431, 30.2672
        x, y = project_to_albers(lon, lat)

        # Should be in reasonable range for CONUS
        @test -2e6 < x < 3e6
        @test -2e6 < y < 4e6

        # Round-trip test
        lon2, lat2 = inverse_project(x, y)
        @test abs(lon2 - lon) < 0.01  # Within 0.01 degrees
        @test abs(lat2 - lat) < 0.01

        # Vectorized version
        lons = [-97.0, -96.0, -95.0]
        lats = [30.0, 31.0, 32.0]
        xs, ys = project_to_albers(lons, lats)
        @test length(xs) == 3
        @test length(ys) == 3
    end

    @testset "Spatial - IDW" begin
        # Simple test: 4 corners, equal distance from center
        x_obs = [0.0, 100.0, 0.0, 100.0]
        y_obs = [0.0, 0.0, 100.0, 100.0]
        z_obs = [10.0, 20.0, 30.0, 40.0]

        # Center point should get mean of all (equal distances)
        x_pred = [50.0]
        y_pred = [50.0]

        z_pred = idw_interpolate(x_obs, y_obs, z_obs, x_pred, y_pred; idp=2, nmax=4)
        @test abs(z_pred[1] - 25.0) < 0.01  # Mean of 10,20,30,40

        # Test with unequal distances
        x_pred2 = [10.0]
        y_pred2 = [10.0]
        z_pred2 = idw_interpolate(x_obs, y_obs, z_obs, x_pred2, y_pred2; idp=2, nmax=4)
        @test z_pred2[1] < 25.0  # Should be closer to 10 (first point)
    end

    @testset "Spatial - Nyquist" begin
        # Regular grid should give expected spacing
        x = [0.0, 100.0, 200.0, 0.0, 100.0, 200.0]
        y = [0.0, 0.0, 0.0, 100.0, 100.0, 100.0]

        spacing = compute_nyquist_spacing(x, y)
        @test spacing ≈ 50.0 atol=1.0  # NN distance = 100, Nyquist = 50
    end

    @testset "Bayesian - Anomaly" begin
        # Anomaly = Baseline - Current
        @test compute_anomaly(100.0, 90.0) == -10.0  # Depletion (deeper)
        @test compute_anomaly(80.0, 90.0) == 10.0   # Recovery (shallower)
        @test compute_anomaly(90.0, 90.0) == 0.0    # No change
    end

    @testset "Bayesian - Storage" begin
        # Test storage estimation
        sy = SyPrior(0.15, 0.04, "Test", :unconfined)
        result = estimate_storage(-10.0, 100000.0, sy; n_samples=1000)

        # ΔS = Sy × Area × Δh / 1000 ≈ 0.15 × 100000 × -10 / 1000 = -150 km³
        @test result.storage_mean_km3 ≈ -150.0 atol=20.0
        @test result.storage_sd_km3 > 0
        @test result.ci_95[1] < result.storage_mean_km3
        @test result.ci_95[2] > result.storage_mean_km3
    end

    @testset "Bayesian - Trend" begin
        years = [2000.0, 2005.0, 2010.0, 2015.0, 2020.0]
        values = [0.0, -1.0, -2.0, -3.0, -4.0]  # -1 m per 5 years = -2 m/decade

        trend = compute_trend(years, values)
        @test trend.slope_decade ≈ -2.0 atol=0.01
        @test trend.r2 ≈ 1.0 atol=0.01  # Perfect linear fit
    end

    @testset "QualityControl - Outliers" begin
        # Normal data with one outlier
        values = [1.0, 2.0, 1.5, 2.5, 1.8, 100.0]  # 100 is outlier

        outliers = detect_outliers(values; z_threshold=3.0)
        @test outliers[6] == true   # 100 is outlier
        @test sum(outliers[1:5]) == 0  # Others are not
    end

    @testset "QualityControl - MAD" begin
        # Test robust MAD
        x = [1.0, 2.0, 3.0, 4.0, 5.0]
        mad = robust_mad(x)
        @test mad ≈ 1.4826 atol=0.01  # MAD of symmetric data

        # With outlier - should be robust
        x_outlier = [1.0, 2.0, 3.0, 4.0, 100.0]
        mad_outlier = robust_mad(x_outlier)
        @test mad_outlier < 10.0  # Should not be affected much by outlier
    end

    @testset "Quality Classification" begin
        @test classify_quality(80, 150, 0.6, 0.3) == :HIGH
        @test classify_quality(60, 50, 0.3, 0.8) == :MEDIUM
        @test classify_quality(30, 10, 0.1, 1.5) == :LOW
    end

end

println("All tests passed!")
