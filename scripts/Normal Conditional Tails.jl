using Random, Statistics

# Set the seed for reproducibility (optional)
Random.seed!(123)

M = 10000

# Generate M observations from the standard normal distribution
X = randn(M)


u = 1.0

# Compute conditional EDF of the upper tail

# Compute the number of observations in the upper tail
n = sum(X .> u)


function cond_edf(x, u, X)
    n = sum(X .> u)
    return sum(X .> x + u) / n
end


cond_edf(1, u, X)


# Extract values of X which are greater than u
Xu = X[X .> u]

X_tail = Xu .- u

mean(X_tail)