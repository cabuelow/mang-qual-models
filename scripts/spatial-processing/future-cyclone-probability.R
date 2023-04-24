# Define the parameters
freq <- 1/10000      # Annual frequency of cyclones
years <- 2050-2023      # Number of years in the period

# Calculate the probability of a cyclone occurring in any given year
p <- freq * 1

# Calculate the probability of the first cyclone occurring by the year 2100
prob <- 1 - (1 - p) ^ years

# Print the result
cat("The probability of the first cyclone occurring by the year 2100 is", round(prob, 4))
