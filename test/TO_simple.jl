
@testset "Transfer operator: simple" begin
    A = rand(10, 3)
    A[1, :] = A[5, :]
    A[2, :] = A[7, :]
    n_bins = 12
    TO = TO_from_binning(A, n_bins)

    @test typeof(TO) == SparseMatrixCSC{Float64, Int64}
end
