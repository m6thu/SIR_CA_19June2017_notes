setwd('~/Desktop/Vaccine_Model')
source("~/Desktop/Vaccine_Model/SIR_CA.R")

## SIS model script for comparing with mass model

########### Initialize ##############
# Default settings loaded from SIR_CA.R
#prob_infection = 0.5
#prob_recover = 0.5
#prob_mutate = 0.5
#max_iter_rec = 5
#iterations = 500

kernel <- genKernel(neighbor = 100, c(1,1,1,1,1,1,1,1,1,1,
                                      1,1,1,1,1,1,1,1,1,1,
                                      1,1,1,1,1,1,1,1,1,1,
                                      1,1,1,1,1,1,1,1,1,1,
                                      1,1,1,1,1,1,1,1,1,1,
                                      1,1,1,1,1,1,1,1,1,1,
                                      1,1,1,1,1,1,1,1,1,1,
                                      1,1,1,1,1,1,1,1,1,1,
                                      1,1,1,1,1,1,1,1,1,1,
                                      1,1,1,1,1,1,1,1,1,1))

plots <- list()
s_n <- vector()
i_n <- vector()
r_n <- vector()
plots[1] <- list(initGrid())
track <- cleanGrid(dim(grid_map)[1], dim(grid_map)[2])

####### Select R0 here ########
R0 = 10
n = dim(grid_map)[1]*dim(grid_map)[2]
prob_infection = R0/(n*max_iter_rec)
prob_infection

######## CA Model #############
for(i in 1:iterations){
    # Get locations to be updated as infected this turn mark as 'i1'
    new_infect_grid <- calcInfect(grid_map, prob_infection, kernel, 'i1')
    # Get grid locations that are 'i1' longer than max_iter_rec timesteps
    output <- calcSusceptible(grid_map, i = max_iter_rec, track, 'i1')
    new_cured <- output[[1]]
    track <- output[[2]]
    
    # Update grid map with new infected for this timestep
    grid_map <- addInfect(grid_map, new_infect_grid, 'i1')
    # Update grid map with new susceptible for this timestep
    grid_map <- addCuredtoSusceptible(grid_map, new_cured, 'i1')
    # Save grid plots for movie
    plots[i+1] <- list(plotGrid(grid_map)) # shift to make space for first frame from initGrid
    
    # Collect summary for plotting
    s_n[i] = sum(grid_map == 0)
    i_n[i] = sum(grid_map == 'i1')
    
}
idx <- 1:length(s_n)
df <- data.frame(idx, s_n, i_n)
png("r010_di5.png")
ggplot(df) + geom_line(aes(x=idx, y=s_n), color='green') + geom_line(aes(x=idx, y=i_n), color='red')
dev.off()

# makeMovie(plots) # Run this to printout jpg and make movies

######### Mass model ###########
png("r010_di5_mass.png")
calcSIStheoretical(beta = R0, mu = 1, t_series = seq(0, 10, 0.001), n = n)
dev.off()
