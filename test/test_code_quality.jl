using QuantumInputOutput, Test

@testset "best practices" begin
    using Aqua

    Aqua.test_ambiguities([QuantumInputOutput]; broken = false)
    Aqua.test_piracies(QuantumInputOutput; broken = false)
    Aqua.test_all(QuantumInputOutput; ambiguities = false, piracies = false)
end

@testset "ExplicitImports" begin
    using ExplicitImports

    @test check_no_implicit_imports(QuantumInputOutput) == nothing
    @test check_all_explicit_imports_via_owners(QuantumInputOutput) == nothing
    @test check_no_stale_explicit_imports(QuantumInputOutput) == nothing
    @test check_all_qualified_accesses_via_owners(QuantumInputOutput) == nothing
    @test check_no_self_qualified_accesses(QuantumInputOutput) == nothing
end

@static if isempty(VERSION.prerelease)
    @testset "Code linting" begin
        using JET

        rep = report_package(
            QuantumInputOutput;
            target_modules = (QuantumInputOutput,),
            ignore_missing_comparison = true,
        )
        @show rep
        @test isempty(JET.get_reports(rep))
    end
end

@testset "Concretely typed" begin
    import QuantumInputOutput as QIO
    using CheckConcreteStructs

    all_concrete(typeof(QIO.SLH(1, 0, 0))) # TODO
    all_concrete(QIO.Gaussian)
    all_concrete(QIO.PulseCoupling)
end
