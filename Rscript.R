#######################################################################################################
#                                                                                                     #
#                       REEF LIFE SURVEY DIVERSITY EFFECTS - PRICE EQUATION                           #
#                                                                                                     #
#######################################################################################################

# Author: Jon Lefcheck
# Last updated: 11 October 2021
# Contact: lefcheckJ@si.edu

#######################################################################################################

# Load required libraries
library(compiler)
library(cowplot)
library(ggdist)
library(grid)
library(gridExtra)
library(piecewiseSEM)
library(randomForest)
library(sp)
library(parallel)
library(plotrix)
library(tidyverse)
library(readxl)

source("https://gist.githubusercontent.com/jslefche/76c076c1c7c5d200e5cb87113cdb9fb4/raw/300ed0c310cd366f2e3d1250c010e390f46fc970/price2.R")

source("https://gist.githubusercontent.com/jslefche/66d461cdf6ffbbc81c31/raw/5a9cb1c1f7ab915b43bf4d96897aec9c3740a377/price.R")

# Compile Price function to increase speed
enableJIT(3)

price2 <- cmpfun(price2) 

price <- cmpfun(price)

#######################################################################################################
#                                IMPORTING AND FORMATTING THE DATA                                    #
#######################################################################################################

# Read in fish data
rls <- read_excel("Source Data.xlsx", sheet = 1)

# Read in environmental data
env <-  read_excel("Source Data.xlsx", sheet = 2)

# Summarize by SiteCode
survey <- rls %>% 
  
  group_by(Realm, SiteCode, SiteLat, SiteLong, Ecoregion, SPECIES_NAME) %>%
  
  summarize(Depth = mean(Depth, na.rm = T), Abundance = mean(Abundance, na.rm = T), Biomass = mean(Biomass, na.rm = T)) %>%
  
  ungroup()

# Store colors for plotting
colors <- data.frame(
  variable = c("RICH-L", "RICH-G", "COMP-L", "COMP-G", "CDE", "DIV"),
  col = c("dodgerblue", "cornflowerblue", "firebrick", "firebrick2", "darkorange", "darkorchid1")
)

#######################################################################################################

# Use fixed distance (kms) to account for potentially variable species pools

# Vector of kms to test from the baseline site
kms <- 100 #c(15, 25, 50, 100)

# Matrix of longlats
longlats <- as.matrix(unique(survey[, c("SiteCode", "SiteLong", "SiteLat")]))[, -1, drop = F]

longlats <- apply(longlats, 2, as.numeric)

# Create vector of site codes ordered by total biomass
sitecodes <- survey %>% group_by(SiteCode) %>% summarize(biomass = sum(Biomass, na.rm = T)) %>%
  
  ungroup() %>% arrange(desc(biomass)) %>% collect() %>% .[["SiteCode"]]

# Loop over kms
km.mat.list <- lapply(kms, function(km) {
  
  # Create vector of site codes to ignore
  ignore <<- c()
  
  # Compute Price components for each site within km radius
  mat.list <- lapply(sitecodes, function(i) {
    
    if(i %in% ignore) NULL else {
      
      # Get latlong of site
      longlat <- unique(subset(survey, SiteCode == i)[, c("SiteLong", "SiteLat")])
      
      # Compute distance between baseline and all other sites
      d <- spDistsN1(longlats, as.numeric(longlat), longlat = T)
      
      # Get vector of comparison sites
      n <- unique(survey$SiteCode)[d <= km]
      
      n <- n[!duplicated(n)]
      
      # Get matrix of biomasses
      bmat <- survey %>% 
        
        filter(SiteCode %in% n) %>% 
        
        filter(!SiteCode %in% ignore) %>%
        
        select(SiteCode, SPECIES_NAME, Biomass) %>%
        
        spread(SPECIES_NAME, Biomass, fill = 0) %>%
        
        as.data.frame()
      
      # Set rownames
      rownames(bmat) <- bmat$SiteCode
      
      bmat <- bmat[, -1, drop = F]
      
      # Get matrix of abundances
      amat <- survey %>% 
        
        filter(SiteCode %in% n) %>% 
        
        filter(!SiteCode %in% ignore) %>%
        
        select(SiteCode, SPECIES_NAME, Abundance) %>%
        
        spread(SPECIES_NAME, Abundance, fill = 0) %>%
        
        as.data.frame()
      
      # Set rownames
      rownames(amat) <- amat$SiteCode
      
      amat <- amat[, -1, drop = F]
      
      # Remove any species with biomass/abundance == 0
      bmat <- bmat[, colSums(bmat) != 0, drop = F]
      
      amat <- amat[, colSums(amat) != 0, drop = F]
      
      # Add SiteCode to ignore
      ignore <<- c(ignore, rownames(amat))
      
      if(nrow(amat) == 0 | nrow(bmat) == 0) NULL else
        
        list(
          amat = amat,
          bmat = bmat
        )
      
    }
    
  } )
  
  mat.list[!sapply(mat.list, is.null)]
  
} )

names(km.mat.list) <- kms

# Get summary statistics on site replication for 100 km
summ <- sapply(km.mat.list[["100"]], function(i) ifelse(nrow(i$bmat) >= 5, nrow(i$bmat) - 1, NA))

mean(summ, na.rm = T); range(summ, na.rm = T)

# Apply Price equation to biomass matrices
price.km.list <- lapply(names(km.mat.list), function(i) {
  
  do.call(rbind, lapply(km.mat.list[[i]], function(j) {
    
    j <- j$bmat
    
    # Remove all matrices with fewer than 5 communities
    if(dim(j)[1] < 5) data.frame() else {
      
      # replace with `price` function for original decomposition
      out <- cbind.data.frame(km = i, type = "observed", price2(j, standardize = TRUE)) 
      
      out <- out[!out$sharedS == 0, ]
      
      return(out)
      
    }
    
  } ) )
  
} )

names(price.km.list) <- kms

# # Summarize by baseline site
# price.km.list <- lapply(names(price.km.list), function(i) {
# 
#   price.km.list[[i]] %>% filter(sharedS != 0) %>%
# 
#     group_by(km, type, baseline) %>%
# 
#     summarize_at(vars(baselineS, comparisonS, sharedS, RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV),
#                  funs(mean), na.rm = T) %>%
# 
#     ungroup()
# 
# } )
# 
# names(price.km.list) <- kms

# Get summary data
lapply(names(price.km.list), function(i) {
  
  # Create summary data for means +/- SEs
  price.km.list[[i]] %>%
    
    summarize_at(vars(RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), lst(mean, std.error), na.rm = T) %>%
    
    gather() %>%
    
    mutate(variable = unlist(strsplit(key, "_mean|_std.error")), type = c(rep("mean", 6), rep("se", 6))) %>%
    
    select(-key) %>%
    
    spread(type, value)
  
} )

# function to compute the geometric mean
gm_mean <- function(x, na.rm = TRUE) exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))

# Get ratio of RICH to COMP, and DIV to CDE with bootstrapped CI's
lapply(names(price.km.list), function(i) {
  
  dat <- price.km.list[[i]] 
  
  # Get vector of ratios
  RICH_COMP_vec <- with(dat, (COMP_L + RICH_L)/RICH_L) 
  
  DIV_CDE_vec <- with(dat, TOTAL_DIV/CDE) # abs(TOTAL_DIV/CDE))
  
  # DIV_CDE_vec <- DIV_CDE_vec[!is.infinite(DIV_CDE_vec)] # for abs to run above
  
  # Sample 5000 times (with replacement) to generate distribution of means
  RICH_COMP_means <- sapply(1:5000, function(i) {
    
    new_vec <- sample(RICH_COMP_vec[RICH_COMP_vec > 0], replace = T)
    
    # Take geometric mean
    gm_mean(new_vec, na.rm = T)
    
  } )
  
  RICH_COMP_means <- RICH_COMP_means[order(RICH_COMP_means)]
  
  # Sample 5000 times (with replacement) to generate distribution of means
  DIV_CDE_means <- sapply(1:5000, function(i) {
    
    new_vec <- sample(DIV_CDE_vec[DIV_CDE_vec >0], replace = T)
    
    # Take geometric mean
    gm_mean(new_vec, na.rm = T)
    
  } )
  
  DIV_CDE_means <- DIV_CDE_means[order(DIV_CDE_means)]
  
  data.frame(
    ratio = c("COMP_RICH", "DIV_CDE"),
    mean = c(exp(mean(log(RICH_COMP_vec[!is.infinite(RICH_COMP_vec)]), na.rm = T)),
             exp(mean(log(DIV_CDE_vec[!is.infinite(DIV_CDE_vec)]), na.rm = T))),
    lower = c(RICH_COMP_means[5000*0.025], DIV_CDE_means[5000*0.025]),
    upper = c(RICH_COMP_means[5000*0.975], DIV_CDE_means[5000*0.975])
  )
  
} )

# Plot results
lapply(as.character(kms), function(i) { # names(price.km.list)
  
  # Melt data so it is shortways (raw components)
  price.km.df.melt <- price.km.list[[i]] %>% 
    
    gather(variable, value, RICH_L:TOTAL_DIV, factor_key = T)
  
  # Create summary data for means +/- SEs
  price.km.df.summary <- price.km.list[[i]] %>%
    
    summarize_at(vars(RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), funs(mean, std.error), na.rm = T) %>%
    
    gather() %>%
    
    mutate(variable = unlist(strsplit(key, "_mean|_std.error")), type = c(rep("mean", 6), rep("se", 6))) %>%
    
    select(-key) %>%
    
    spread(type, value)
  
  # Re-level for better plotting
  price.km.df.melt$variable <- factor(price.km.df.melt$variable, levels = c("RICH_L", "COMP_L", "RICH_G", "COMP_G", "CDE", "TOTAL_DIV"))
  
  levels(price.km.df.melt$variable) <- c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV")
  
  price.km.df.summary$variable <- factor(price.km.df.summary$variable, levels = c("RICH_L", "COMP_L", "RICH_G", "COMP_G", "CDE", "TOTAL_DIV"))
  
  levels(price.km.df.summary$variable) <- c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV")
  
  # Generate plots where x-axis is free
  price.km.plot.list <- lapply(levels(price.km.df.melt$variable), function(j) {
    
    ylimit <- max(density(na.omit(subset(price.km.df.melt, variable == j)$value))$y) * 0.2
    
    p <- ggplot(subset(price.km.df.melt, variable == j)) +
      geom_density(aes(value, fill = variable), col = NA) +
      geom_rug(aes(value), col = "grey80") +
      geom_vline(xintercept = 0, lwd = 0.5) +
      geom_point(
        data = subset(price.km.df.summary, variable == j),
        aes(y = ylimit, x = mean), size = 2) +
      geom_errorbarh(
        data = subset(price.km.df.summary, variable == j),
        aes(y = ylimit, x = mean, xmin = mean - 2 * se, xmax = mean + 2 * se),
        height = 0) +
      coord_flip() +
      scale_fill_manual(values = as.character(subset(colors, variable == j)$col)) +
      scale_x_continuous(limits = c(-1, 1)) +
      labs(y = j, x = "Standardized\nloss of biomass") +
      theme_bw(base_size = 11) + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x = unit(0, "lines")
      )
    
    if(j != "RICH-L") {
      
      p <- p + theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
      )
      
    }
    
    if(j %in% levels(price.km.df.melt$variable)[-1]) {
      
      if(j == "COMP-L")
        
        p <- p + theme(plot.margin = unit(c(1, 1, 1, -0.3), "cm")) else
          
          p <- p + theme(plot.margin = unit(c(1, 1, 1, -1.1), "cm"))
        
    }
    
    # Turn off clipping
    p <- ggplot_gtable(ggplot_build(p))
    p$layout$clip[p$layout$name == "panel"] <- "off"
    
    return(p)
    
  } )
  
  pdf(paste0("Figure 3_", i, "km-new.PDF"),
      width = 7.2, height = 2.5)
  bquiet = print(grid.draw(do.call(cbind, c(price.km.plot.list, size = "last"))))
  dev.off() 
  
} )

dev.off()

# Mean raw change in diversity +/- SE
mean(with(price.km.list[["100"]], comparisonS - baselineS))

std.error(with(price.km.list[["100"]], comparisonS - baselineS))

# Mean proportional change in diversity +/- SE
mean(with(price.km.list[["100"]], (comparisonS - baselineS) / baselineS)*100)

std.error(with(price.km.list[["100"]], (comparisonS - baselineS) / baselineS)*100)

# Mean number of shared species +/- SE
mean(price.km.list[["100"]]$sharedS)

std.error(price.km.list[["100"]]$sharedS)

# Mean proportion of of shared species +/- SE
sapply(names(price.km.list), function(x) mean(with(price.km.list[[x]], sharedS/baselineS)))

std.error(with(price.km.list[["100"]], sharedS/baselineS))

# Plot map of 100 km sites

# Get lat-longs of included sites
final_sites <- unique(unlist(unique(price.km.list[["100"]][, c("baseline", "comparison")])))

final_longlats <- survey %>% filter(SiteCode %in% final_sites) %>% group_by(SiteCode) %>%
  
  summarize(SiteLat = mean(SiteLat), SiteLong = mean(SiteLong))

map.world <- map_data("world")

(map <- ggplot() +
    geom_polygon(data = map.world, aes(x = long, y = lat, group = group), fill = "cornsilk") +
    geom_point(data = final_longlats, aes(x = SiteLong, y = SiteLat), 
               shape = 21, col = "black", fill = "cornflowerblue", size = 1) +
    # theme_map() +
    labs(x = "", y = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "lightblue1"),
      panel.border = element_rect(fill = NA)
    )
)

ggsave("Figure 2.PDF", map, width = 7, height = 4)

#######################################################################################################

# Conduct null simulations

library(parallel)

no_cores <- detectCores()

# Initiate cluster
cl <- makeCluster(no_cores)

clusterExport(cl, c("km.mat.list", "randomizeMat", "price2"))

# Shuffle observations and recompute Price equations
null.km.list <- lapply("100", function(i) { # names(km.mat.list)
  
  do.call(rbind, parLapply(cl, 1:length(km.mat.list[[i]]), function(j) { # length(km.mat.list[[i]])
    
    amat <- km.mat.list[[i]][[j]]$amat
    
    bmat <- km.mat.list[[i]][[j]]$bmat
    
    amat <- amat[rownames(amat) %in% rownames(bmat), colnames(amat) %in% colnames(bmat)]
    
    if(dim(amat)[1] < 5) data.frame() else {
      
      # Generate null matrix
      ret <- do.call(rbind, lapply(1:1000, function(k) {
        
        newmat <- randomizeMat(amat, bmat)
        
        price2(newmat)
        
      } ) )
      
      write.csv(ret, paste("Price null 100 km", j, ".csv"))
      
    }
    
  } ) ) 
  
} )

stopCluster(cl)

# Read in .csv's
temp <- list.files(".", pattern = "*.csv")

null.list <- lapply(temp, function(x) {
  
  error <- try(read.csv(paste0(".", x)))
  
  if(class(error) == "try-error") data.frame() else {
    
    error
    
  }
  
} )

# Compute means by iteration
null.split.list <- lapply(null.list, function(i) {
  
  if(nrow(i) == 0) NULL else {
    
    # Identify which is the last comparison in the sequence
    last <- i[which(i$comparison == 1) - 1, "comparison"][1] - 1
    
    # Split into list every certain number of rows
    i.list <- split(i, rep(1:ceiling(nrow(i)/last), each = last, length.out = nrow(i)))
    
  }
  
} )

null.split.list <- null.split.list[!sapply(null.split.list, is.null)]

null.bind.list <- transpose(null.split.list) %>% map(bind_rows)

# Compute means for each iteration
null.df <- do.call(rbind, lapply(null.bind.list, function(x) colMeans(x[, 8:13])))

# Melt null data so it is shortways (raw components)
price.null.df.melt <- null.df %>% 
  
  as.data.frame() %>%
  
  # filter(sharedS != 0) %>%
  
  # select(-comparison, -baseline, -deltaFunc) %>%
  
  gather(variable, value, RICH_L:TOTAL_DIV, factor_key = T) %>%
  
  na.omit() 

# Create summary data for observed means +/- SEs
price.km.df.summary <- price.km.list[["100"]] %>%
  
  filter(sharedS != 0) %>%
  
  summarize_at(vars(RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), funs(mean, std.error)) %>%
  
  gather() %>%
  
  mutate(variable = unlist(strsplit(key, "_mean|_std.error")), type = c(rep("mean", 6), rep("se", 6))) %>%
  
  select(-key) %>%
  
  spread(type, value) %>%
  
  mutate(n = nrow(price.km.list[["100"]]))

# Create summary for null values
price.km.null.df.summary <- null.df %>%
  
  as.data.frame() %>%
  
  # filter(sharedS != 0) %>%
  
  summarize_at(vars(RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), funs(mean, std.error), na.rm = T) %>%
  
  gather() %>%
  
  mutate(variable = unlist(strsplit(key, "_mean|_std.error")), type = c(rep("mean", 6), rep("se", 6))) %>%
  
  select(-key) %>%
  
  spread(type, value) %>%
  
  mutate(n = nrow(null.df))

# Re-level for better plotting
levels(price.null.df.melt$variable) <- c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV")

price.null.df.melt$variable <- factor(price.null.df.melt$variable, levels = c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV"))

price.km.df.summary$variable <- c("CDE", "COMP-G", "COMP-L", "RICH-G", "RICH-L", "DIV")

price.km.df.summary$variable <- factor(price.km.df.summary$variable, levels = c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV"))

price.km.null.df.summary$variable <- c("CDE", "COMP-G", "COMP-L", "RICH-G", "RICH-L", "DIV")

price.km.null.df.summary$variable <- factor(price.km.null.df.summary$variable, levels = c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV"))

# Generate plots where x-axis is free
price.null.plot.list <- lapply(levels(price.null.df.melt$variable), function(j) {
  
  ylimit <- max(density(subset(price.null.df.melt, variable == j)$value)$y) * 0.2
  
  p <- ggplot(subset(price.null.df.melt, variable == j)) +
    geom_density(aes(value, fill = variable), col = NA) +
    geom_rug(aes(value), col = "grey80") +
    geom_vline(xintercept = 0, lwd = 0.5) +
    geom_point(
      data = subset(price.km.null.df.summary, variable == j),
      aes(y = ylimit, x = mean), pch = 21, fill = "white", size = 2) +
    geom_point(
      data = subset(price.km.df.summary, variable == j),
      aes(y = ylimit, x = mean), size = 2) +
    geom_errorbarh(
      data = subset(price.km.df.summary, variable == j),
      aes(y = ylimit, x = mean, xmin = mean - 2 * se, xmax = mean + 2 * se),
      height = 0) +
    coord_flip() +
    scale_fill_manual(values = as.character(subset(colors, variable == j)$col)) +
    scale_x_continuous(limits = c(-0.6, 0.2)) +
    labs(y = j, x = "Standardized\nloss of biomass") +
    theme_bw(base_size = 11) + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.ticks = element_blank(), 
      axis.text.x = element_blank(),
      strip.text = element_blank(),
      panel.spacing.x = unit(0, "lines")
    )
  
  if(j != "RICH-L") {
    
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )
    
  }
  
  if(j %in% levels(price.null.df.melt$variable)[-1]) {
    
    if(j == "COMP-L")
      
      p <- p + theme(plot.margin = unit(c(1, 1, 1, -0.3), "cm")) else
        
        p <- p + theme(plot.margin = unit(c(1, 1, 1, -1.1), "cm"))
      
  }
  
  # Turn off clipping
  p <- ggplot_gtable(ggplot_build(p))
  p$layout$clip[p$layout$name == "panel"] <- "off"
  
  return(p)
  
} )

pdf(paste0("Figure S2_revised.PDF"),
    width = 7.2, height = 2.5)
bquiet = print(grid.draw(do.call(cbind, c(price.null.plot.list, size = "last"))))
dev.off() 

# Compuate one-tailed tests
price.km.df.summary <- price.km.list[["100"]] %>%
  
  filter(sharedS != 0) %>%
  
  summarize_at(vars(RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), funs(mean, sd)) %>%
  
  gather() %>%
  
  mutate(variable = unlist(strsplit(key, "_mean|_sd")), type = c(rep("mean", 6), rep("sd", 6))) %>%
  
  select(-key) %>%
  
  spread(type, value) %>%
  
  mutate(n = nrow(price.km.list[["100"]]))

price.km.null.df.summary <- null.df %>%
  
  as.data.frame() %>%
  
  # filter(sharedS != 0) %>%
  
  summarize_at(vars(RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), funs(mean, sd), na.rm = T) %>%
  
  gather() %>%
  
  mutate(variable = unlist(strsplit(key, "_mean|_sd")), type = c(rep("mean", 6), rep("sd", 6))) %>%
  
  select(-key) %>%
  
  spread(type, value) %>%
  
  mutate(n = nrow(null.df))

lapply(unique(price.km.df.summary$variable), function(i) {
  
  x1 <- subset(price.km.df.summary, variable == i)
  
  x2 <- subset(price.km.null.df.summary, variable == i)
  
  df <- x1$n + x2$n - 2
  
  t <- (x1$mean - x2$mean) / ((sqrt(((x1$n - 1)*x1$sd + (x2$n - 1)*x2$sd)/df)) * sqrt((1/x1$n) + (1/x2$n)))
  
  pt(t, df, lower.tail = TRUE)
  
} )

#######################################################################################################

# Compute median sizeclass
median(rls$Sizeclass)

# Collapse some size bins
rls$new.sizeclass <- 
  ifelse(rls$Sizeclass < 10, "0-10",
         ifelse(rls$Sizeclass >= 10 & rls$Sizeclass < 30, "10-30",
                ifelse(rls$Sizeclass >= 30 & rls$Sizeclass < 50, "30-50",
                       ifelse(rls$Sizeclass >= 50 & rls$Sizeclass < 100, "50-100",
                              ifelse(rls$Sizeclass >= 100 & rls$Sizeclass < 200, "100-200", "200-400")))))

# Compute proportional difference biomasses in different size classes between baseline and comparison sites
sizeclass.prop.df <- do.call(rbind, lapply((1:nrow(price.km.list[["100"]])), function(i) {
  
  # Get biomass in baseline and comparison communities
  base <- rls[rls$SiteCode == price.km.list[["100"]][i, "baseline"], , drop = F]
  
  comp <- rls[rls$SiteCode == price.km.list[["100"]][i, "comparison"], , drop = F]
  
  # Get vector of species shared betwene baseline and comparison sites
  shared <- unique(base$SPECIES_NAME[which(base$SPECIES_NAME %in% comp$SPECIES_NAME)])
  
  # Get shared and unique biomass in the baseline site
  base.shared <- base %>% filter(SPECIES_NAME %in% shared) %>% group_by(new.sizeclass) %>% summarize(biomass = sum(Biomass)) %>%
    mutate(biomass = biomass / sum(base$Biomass))
  
  base.unique <- base %>% 
    filter(!SPECIES_NAME %in% shared) %>%
    group_by(new.sizeclass) %>% 
    summarize(biomass = sum(Biomass)) %>%
    mutate(biomass = biomass / sum(base$Biomass)) 
  
  # Get shared and unique biomass in the comparison site
  comp.shared <- comp %>% filter(SPECIES_NAME %in% shared) %>% group_by(new.sizeclass) %>% summarize(biomass = sum(Biomass)) %>%
    mutate(biomass = biomass / sum(comp$Biomass))
  
  comp.unique <- comp %>% 
    filter(!SPECIES_NAME %in% shared) %>%
    group_by(new.sizeclass) %>% 
    summarize(biomass = sum(Biomass)) %>%
    mutate(biomass = biomass / sum(comp$Biomass)) 
  
  # Bind and return
  if(nrow(base.unique) != 0 & nrow(comp.unique) != 0) {
    
    rbind(
      cbind.data.frame(base.shared, site = "baseline", type = "shared"),
      cbind.data.frame(base.unique, site = "baseline", type = "unique"),
      cbind.data.frame(comp.shared, site = "comparison", type = "shared"),
      cbind.data.frame(comp.unique, site = "comparison", type = "unique")
    )
    
  } else data.frame()
  
} ) )

sizeclass.prop.df$new.sizeclass <- factor(sizeclass.prop.df$new.sizeclass, levels = c("0-10", "10-30", "30-50", "50-100", "100-200", "200-400"))

# Get unique plots
(p.lefta <- ggplot(subset(sizeclass.prop.df, type == "unique" & site == "comparison"), 
                   aes(y = biomass, x = new.sizeclass)) +
    stat_eye(fill = "grey80", col = "dodgerblue4", lwd = 0.5) +
    theme_bw(base_size = 8) +
    xlab("Size class (cm)") +
    ggtitle("a. Unique species") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(), 
      plot.title = element_text(size = 10, face = "bold"),
      plot.margin = unit(c(1, 0, 1, 0), "mm")) +
    scale_y_reverse(breaks = c(0.2, 0.5, 0.8)) + #limits = c(0.6, 0), expand = c(0, 0), breaks = c(0.3, 0.5)) +
    coord_flip()
)

(p.righta <- ggplot(subset(sizeclass.prop.df, type == "unique" & site == "baseline"), 
                    aes(y = biomass, x = new.sizeclass)) +
    stat_eye(fill = "grey80", col = "dodgerblue2", lwd = 0.5) +
    scale_y_continuous(breaks = c(0.2, 0.5, 0.8)) + #limits = c(0, 0.6), expand = c(0, 0), breaks = c(0, 0.3, 0.5)) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 8),
      plot.margin = unit(c(1, 1, 1, -1), "mm")) +
    coord_flip() 
)

fig5a <- cowplot::plot_grid(p.lefta, p.righta,  align = "h", rel_widths = c(1, 0.70))

fig5a <- ggdraw(add_sub(fig5a, "Proportional biomass", vpadding = grid::unit(0, "lines"), size = 8, 
                        y = 0, x = 0.6, vjust = -0.5))

# Create legend
legend.a <- data.frame(y = 1:2, fill = c("Reference", "Comparison"))

legend.a$fill <- factor(legend.a$fill, levels = c("Reference", "Comparison"))

plegend.a <- ggplot(legend.a,
                    aes(x = fill, y = y, fill = fill)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = c("dodgerblue2", "dodgerblue4"), name = "",
                                                  guide = guide_legend(reverse = TRUE)) +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"))

plegend.a <- ggdraw(get_legend(plegend.a))

fig5a <- plot_grid(fig5a, NULL, plegend.a, ncol = 1, rel_heights = c(1, 0.01, .1))

# Get shared plots

(p.leftc <- ggplot(subset(sizeclass.prop.df, type == "shared" & site == "comparison"), 
                   aes(y = biomass, x = new.sizeclass)) +
    stat_eye(fill = "grey80", col = "firebrick4", lwd = 0.5) +
    theme_bw(base_size = 8) +
    xlab("Size class (cm)") +
    ggtitle("b. Shared species") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(), 
      plot.title = element_text(size = 10, face = "bold"),
      plot.margin = unit(c(1, 0, 1, 0), "mm")) +
    scale_y_reverse(breaks = c(0.2, 0.5, 0.8)) + #limits = c(0.6, 0), expand = c(0, 0), breaks = c(0.3, 0.5)) +
    coord_flip()
)

(p.rightd <- ggplot(subset(sizeclass.prop.df, type == "shared" & site == "baseline"), 
                    aes(y = biomass, x = new.sizeclass)) +
    stat_eye(fill = "grey80", col = "firebrick2", lwd = 0.5) +
    scale_y_continuous(breaks = c(0.2, 0.5, 0.8)) + #limits = c(0, 0.6), expand = c(0, 0), breaks = c(0, 0.3, 0.5)) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_blank(), 
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 8),
      plot.margin = unit(c(1, 1, 1, -1), "mm")) +
    coord_flip() 
)

fig5b <- cowplot::plot_grid(p.leftc, p.rightd,  align = "h", rel_widths = c(1, 0.70))

fig5b <- ggdraw(add_sub(fig5b, "Proportional biomass", vpadding = grid::unit(0, "lines"), size = 8, 
                        y = 0, x = 0.6, vjust = -0.5))

# Create legend
legend.b <- data.frame(y = 1:2, fill = c("Reference", "Comparison"))

legend.b$fill <- factor(legend.b$fill, levels = c("Reference", "Comparison"))

plegend.b <- ggplot(legend.b, aes(x = fill, y = y, fill = fill)) + 
  geom_bar(stat = "identity") + scale_fill_manual(values = c("firebrick2", "firebrick4"), name = "",
                                                  guide = guide_legend(reverse = TRUE)) +
  theme_bw(base_size = 8) + theme(legend.position = "bottom", legend.key.size = unit(0.3, "cm"))

plegend.b <- ggdraw(get_legend(plegend.b))

fig5b <- plot_grid(fig5b, NULL, plegend.b, ncol = 1, rel_heights = c(1, 0.01, .1))

# Merge into single figure
fig5 <- plot_grid(fig5a, fig5b, nrow = 1)

ggsave("Figure 5.PDF", fig5, width = 5.75, height = 3.75)

# remove y-axis on panel b in photoshop

#######################################################################################################

# Analysis of environmental drivers
price.env.df <- price.km.list[["100"]]

# Bring in environmental covariates
names(price.env.df)[3] <- "SiteCode"

price.env.df <- left_join(price.env.df, env[, -1]) 

# Bring in depth
price.env.df <- cbind(price.env.df, rls[match(price.env.df$SiteCode, rls$SiteCode), c("Depth")])

f <- paste(c("chlomean", "phos", "pop500", "sstmean",  
             "sstrange", "parmean", "salinity" ,"dissox" ,"nitrate", "sstmin",  
             "sstmax", "chlomax", "chlomin", "Depth"), collapse = " + ")

comps <- c("RICH_L", "RICH_G", "COMP_L", "COMP_G", "CDE", "TOTAL_DIV")

rf.list <- lapply(comps, function(i) {
  
  form <- formula(paste(i, " ~ ", f))
  
  randomForest(form, data = na.omit(price.env.df), importance = TRUE)
  
} )

names(rf.list) <- comps

# Write function to extract variable importance and graph using ggplot
varImpPlot.ggplot2 <- function(rf) {
  
  imp <- as.data.frame(importance(rf))
  
  imp$pred <- rownames(imp)
  
  imp <- imp[order(imp$pred), ]
  
  imp$pred <- recode(imp$pred, 
                     "sstmax" = "Max SST",
                     "SiteLong" = "Longitude",
                     "chlomin" = "Min chl-a",
                     "sstmean" = "Mean SST",
                     "parmean" = "Mean PAR",
                     "sstrange" = "Range SST",
                     "dissox" = "DO",
                     "sstmin" = "Min SST",
                     "chlomax" = "Max chl-a",
                     "nitrate" = "Nitrate",
                     "chlomean" = "Mean chl-a",
                     "salinity" = "Salinity",
                     "phos" = "Phosphorus",
                     "SiteLat" = "Latitude",
                     "pop500" = "Human pop index" )
  
  imp <- imp[order(imp$`%IncMSE`), ]
  
  imp$pred <- factor(imp$pred, levels = imp$pred)
  
  ggplot(imp, aes(x = `%IncMSE`, y = pred)) +
    # geom_hline(yintercept = 1:length(imp$pred), col = "grey90") +
    geom_point(size = 2) +
    theme_bw(base_size = 18) +
    labs(x = "% Inc. MSE", y = "") +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
}  

# Plot results
rf.plots.list <- lapply(comps, function(i) {
  
  p <- varImpPlot.ggplot2(rf.list[[i]])
  
  p + labs(title = i)
  
} )

pdf("Figure S6.PDF",
    width = 10, height = 12)
bquiet = print(cowplot::plot_grid(plotlist = rf.plots.list, ncol = 2))
dev.off() 

#######################################################################################################

# Model population size against biomass and species richness

sitecodes <- c(as.character(price.km.list[["100"]]$baseline), as.character(price.km.list[["100"]]$comparison))

sitecodes <- sitecodes[!duplicated(sitecodes)]

# Summarize RLS data by site
rls.summary <- rls %>%
  
  group_by(Ecoregion, SiteCode) %>%
  
  filter(SiteCode %in% sitecodes) %>%
  
  summarize(richness = length(unique(SPECIES_NAME)), 
            biomass = sum(Biomass),
            biomass100 = sum(Biomass[Sizeclass >= 100]),
            size.class = max(Sizeclass))

# Bring in pop'n data
rls.summary <- left_join(rls.summary, env) 

# Run models
biomass.mod <- lme(log10(biomass + 0.1) ~ log10(pop500 + 0.1) + log10(richness + 0.1), random = ~ 1 | Ecoregion, na.action = na.omit, data = rls.summary)

plot(biomass.mod)

coefs(biomass.mod)

rsquared(biomass.mod)

sizeclass.mod <- lme(log10(size.class + 0.1) ~ log10(pop500 + 0.1) + log10(richness + 0.1), random = ~ 1 | Ecoregion, na.action = na.omit, data = rls.summary)

plot(sizeclass.mod)

coefs(sizeclass.mod)

rsquared(sizeclass.mod)

richness.mod <- lme(log10(richness + 0.1) ~ log10(pop500 + 0.1) + log10(biomass + 0.1), random = ~ 1 | Ecoregion, na.action = na.omit, data = rls.summary)

plot(richness.mod)

coefs(richness.mod)

rsquared(richness.mod)

# Gather for plotting
rls.pop500 <- rls.summary %>%
  
  filter(!is.na(pop500)) %>%
  
  group_by(SiteCode) %>%
  
  # mutate(biomass = log10(biomass), biomass100 = log10(biomass100 + 1e3)) %>%
  
  # mutate(biomass = biomass/1000, biomass100 = biomass100/1000) %>%
  
  select(SiteCode, pop500, biomass, biomass100, size.class, richness) %>%
  
  gather(variable, value, biomass, biomass100, size.class, richness)

rls.pop500$variable <- factor(rls.pop500$variable, levels = c("biomass", "biomass100", "size.class", "richness"))

# Plot results
(fig3a <- ggplot(subset(rls.pop500, variable == "biomass"), aes(x = log10(pop500 + 0.1), y = log10(value + 0.1))) +
    geom_point(size = 1, col = "grey75", shape = 1) +
    stat_smooth(method = "lm", col = "black", lwd = 1, se = F) +
    geom_text(x = -Inf, y = Inf, label = "A", hjust = -1, vjust = 1.5, size = 4, fontface = "bold") +
    labs(y = bquote(log[10](Biomass)), x = "") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(size = 10),
      plot.margin = unit(c(0.5, 0.5, 0, 0), "mm")
    )
)

(fig3b <- ggplot(subset(rls.pop500, variable == "size.class"), aes(x = log10(pop500 + 0.1), y = log10(value + 0.1))) +
    geom_point(size = 1, col = "grey75", shape = 1) +
    stat_smooth(method = "lm", formula = y ~ x, col = "black", lwd = 1, se = F) +
    geom_text(x = -Inf, y = Inf, label = "B", hjust = -1, vjust = 1.5, size = 4, fontface = "bold") +
    labs(y = bquote(log[10](Size~class)), x = bquote(log[10](Human~"population"))) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(size = 10),
      plot.margin = unit(c(0.5, 0.5, 0, 0.5), "mm")
    )
)

(fig3c <- ggplot(subset(rls.pop500, variable == "richness"), aes(x = log10(pop500 + 0.1), y = log10(value + 0.1))) +
    geom_point(size = 1, col = "grey75", shape = 1) +
    stat_smooth(method = "lm", formula = y ~ x, col = "black", lwd = 1, se = F) +
    geom_text(x = -Inf, y = Inf, label = "C", hjust = -1, vjust = 1.5, size = 4, fontface = "bold") +
    labs(y = bquote(log[10](Richness)), x = "") +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text = element_text(size = 10),
      plot.margin = unit(c(0.5, 0, 0, 0.5), "mm")
    )
)

pdf("Figure 4.PDF",
    width = 6.5, height = 2)
bquiet = print(plot_grid(fig3a, fig3b, fig3c, align = "h", nrow = 1))
dev.off() 

#######################################################################################################

# Remove non-site attached pelagic species & re-run analyses
pelagicsp <- c("Aldrichetta forsteri", "Arctocephalus pusillus", "Arripis georgianus", "Arripis spp.", 
               "Arripis trutta", "Arripis truttaceus", "Atherinason hepsetoides", "Atule mate", 
               "Carangoides orthogrammus", "Caranx ignobilis", "Caranx sexfasciatus", 
               "Carcharhinus galapagensis", "Carcharhinus melanopterus", "Carcharhinus obscurus", 
               "Carcharodon carcharias", "Chelonia mydas", "Clupeid spp.", "Decapterus macarellus", 
               "Decapterus muroadsi", "Engraulis australis", "Gnathanodon speciosus", "Hyperlophus vittatus",
               "Hyporhamphus australis", "Hyporhamphus melanochir", "Leptatherina presbyteroides", 
               "Neophoca cinerea", "Pomatomus saltatrix", "Pseudocaranx georgianus", 
               "Pseudocaranx georgianus", "Sarda australis", "Sardinops neopilchardus", 
               "Scomber australasicus", "Scomberoides commersonnianus", "Seriola hippos", 
               "Seriola lalandi", "Seriolella brama", "Sphyraena novaehollandiae", "Sphyraena obtusata", 
               "Spratelloides robustus", "Thyrsites atun", "Trachurus declivis", "Trachurus novaezelandiae", 
               "Triaenodon obesus", "Tursiops truncatus", "Tylosurus crocodilus", "Phyllorhiza punctata", 
               "Decapterus koheru", "Pseudocaranx georgianus", "Arripis georgianus", "Pseudocaranx wrighti",
               "Scomberoides commersonnianus", "Trachinotus baillonii", "Trachinotus blochii", 
               "Galeocerdo cuvier", "Aurelia aurita", "Atherinid spp.", "Negaprion acutidens", 
               "Carcharhinus amblyrhynchos", "Herklotsichthys castelnaui", "Chanos chanos", 
               "Hyporhamphus affinis", "Tylosurus gavialoides", "Strongylura incisa", 
               "Atherinomorus vaigiensis", "Echeneis naucrates", "Carangoides chrysophrys", 
               "Selaroides leptolepis", "Carangoides gymnostethus", "Seriola dumerili",
               "Elagatis bipinnulata", "Carangoides fulvoguttatus", "Scomberoides lysan", 
               "Ulua mentalis", "Caranx melampygus", "Seriola rivoliana", "Caranx lugubris", 
               "Caranx papuensis", "Carangoides ferdau", "Carangoides plagiotaenia", 
               "Sphyraena barracuda", "Scomberomorus commerson", "Cybiosarda elegans", 
               "Grammatorcynus bilineatus", "Caretta caretta", "Eretmochelys imbricata", 
               "Aipysurus laevis", "Megaptera novaeangliae", "Tursiops aduncus", "Arctocephalus forsteri", 
               "Carangoides spp.", "Caranx spp.", "Carcharhinus spp.", "Seriola spp.", "Spratelloides spp.",
               "Chrysaora spp.", "Pseudocaranx spp.", "Belonid spp.", "Carangid spp.", "Cheloniid spp.", 
               "Delphinid spp.", "Scombrid spp.", "Sphyraenid spp.", "Cubozoa spp.", "Pinnipedia spp.",
               "Sphyraena flavicauda", "Carcharhinus albimarginatus", "Carcharhinus limbatus", 
               "Sphyrna lewini", "Mobula thurstoni", "Manta birostris", "Remora remora", "Alectis ciliaris",
               "Carangoides bajad", "Caranx crysos", "Sphyraena forsteri", "Thunnus albacares", 
               "Sarda orientalis", "Thunnus obesus", "Gymnosarda unicolor", "Mola mola", 
               "Laticauda laticaudata", "Sphyraena jello", "Euthynnus affinis", "Grammatorcynus spp.", 
               "Rhincodon typus", "Rastrelliger kanagurta", "Echeneis neucratoides", "Elops affinis", 
               "Sphyraena ensis", "Sphyraena idiastes", "Sphyraena qenie", "Sphyraena viridensis", 
               "Champsocephalus esox", "Scomber colias", "Scomberomorus regalis", "Scomberomorus sierra",
               "Belone belone", "Carangoides bartholomaei", "Carangoides ruber", "Caranx caballus",
               "Caranx caninus", "Hemiramphus depauperatus", "Hemiramphus saltator", "Euthynnus lineatus",
               "Trachinotus ovatus", "Trachinotus rhodopus", "Trachinotus stilbe", "Uraspis helvola",
               "Manta birostris", "Carangoides ruber", "Scomberomorus maculatus", "Inermia vittata", 
               "Phoca vitulina", "Atherina presbyter", "Decapterus spp.", "Sphyraena spp.", 
               "Spratelloides gracilis", "Arctocephalus spp.", "Exocoetid spp.", "Otariid spp.", 
               "Pelagia noctiluca", "Sardinella aurita", "Sardina pilchardus", "Sphyraena putnamae", 
               "Sarda sarda", "Megalops atlanticus", "Trachurus trachurus", "Liza richardsonii", 
               "Scomberomorus brasiliensis", "Dichistius capensis", "Otaria flavescens", 
               "Trachurus japonicus", "Carangoides bajad", "Caranx crysos", "Aequorea coerulescens",
               "Grammatorcynus bicarinatus", "Dugong dugon", "Carcharhinus plumbeus", 
               "Pseudocaranx cheilio", "Selar crumenophthalmus", "Scyphozoa spp.", "Aurelia spp.", 
               "Chrysaora pacifica", "Thunnus spp.", "Arctocephalus galapagoensis", "Liza aurata", 
               "Hyporhamphus dussumieri", "Atherina hepsetus", "Selene spp.", "Hydrophis elegans",
               "Selar boops", "Emydocephalus annulatus", "Belone spp.", "Labracoglossa nitida",
               "Arripis xylabion", "Scomberoides tol", "Scomberomorus queenslandicus", 
               "Dipterygonotus balteatus", "Dicentrarchus labrax", "Sphyraena helleri", 
               "Alepes vari", "Hypoatherina barnesi", "Iso rhothophilus", "Carangoides oblongus", 
               "Acanthocybium solandri")

# How many species in the surveys fall into the pelagics?
length(unique(subset(survey, SPECIES_NAME %in% pelagicsp)$SPECIES_NAME))/
  length(unique(survey$SPECIES_NAME))

# Loop over kms
nopelagic.km.mat.list <- lapply(100, function(km) {
  
  # Create vector of site codes to ignore
  ignore <- c()
  
  # Compute Price components for each site within km radius
  mat.list <- lapply(sitecodes, function(i) {
    
    if(i %in% ignore) data.frame() else {
      
      # Get latlong of site
      longlat <- unique(subset(survey, SiteCode == i)[, c("SiteLong", "SiteLat")])
      
      # Compute distance between baseline and all other sites
      d <- spDistsN1(longlats, as.numeric(longlat), longlat = T)
      
      # Get vector of comparison sites
      n <- unique(survey$SiteCode)[d <= km]
      
      n <- n[!duplicated(n)]
      
      # Get matrix of biomasses
      mat <- survey %>% 
        
        filter(SiteCode %in% n) %>% 
        
        filter(!SiteCode %in% ignore) %>%
        
        filter(!SPECIES_NAME %in% pelagicsp) %>%
        
        select(SiteCode, SPECIES_NAME, Biomass) %>%
        
        spread(SPECIES_NAME, Biomass, fill = 0) %>%
        
        as.data.frame()
      
      # Set rownames
      rownames(mat) <- mat$SiteCode
      
      mat <- mat[, -1, drop = F]
      
      # Remove any species with biomass == 0
      mat <- mat[, colSums(mat) != 0, drop = F]
      
      # Add SiteCode to ignore
      ignore <<- c(ignore, rownames(mat))
      
      return(mat)
      
    }
    
  } )
  
  # Remove empty data.frames from list
  mat.list[!sapply(mat.list, function(x) dim(x)[1] == 0)]
  
} )

names(nopelagic.km.mat.list) <- "100"

# Apply Price equation to biomass matrices
nopelagic.price.km.list <- lapply(names(nopelagic.km.mat.list), function(i) {
  
  do.call(rbind, lapply(nopelagic.km.mat.list[[i]], function(j) {
    
    # Remove all matrices with fewer than 5 communities
    if(dim(j)[1] < 5) data.frame() else {
      
      out <- cbind.data.frame(km = i, type = "observed", price2(j, standardize = TRUE))
      
      out <- out[!out$sharedS == 0, ]
      
      return(out)
      
    }
    
  } ) )
  
} )

names(nopelagic.price.km.list) <- "100"

# Summarize by baseline site
nopelagic.price.km.list <- lapply(names(nopelagic.price.km.list), function(i) {
  
  nopelagic.price.km.list[[i]] %>% filter(sharedS != 0) %>% 
    
    group_by(km, type, baseline) %>%
    
    summarize_at(vars(baselineS, comparisonS, sharedS, RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), 
                 funs(mean), na.rm = T) %>%
    
    ungroup()
  
} )

names(nopelagic.price.km.list) <- "100"

# Melt null data so it is shortways (raw components)
nopelagic.price.baseline.df.melt <- nopelagic.price.km.list[["100"]] %>% 
  
  # filter(sharedS != 0) %>%
  
  # select(-comparison, -baseline, -deltaFunc) %>%
  
  gather(variable, value, RICH_L:TOTAL_DIV, factor_key = T) %>%
  
  na.omit() 

# Create summary data for observed means +/- SEs
nopelagic.price.baseline.df.summary <- nopelagic.price.km.list[["100"]] %>%
  
  summarize_at(vars(RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), funs(mean, std.error)) %>%
  
  gather() %>%
  
  mutate(variable = unlist(strsplit(key, "_mean|_std.error")), type = c(rep("mean", 6), rep("se", 6))) %>%
  
  select(-key) %>%
  
  spread(type, value)

# Create summary data for means +/- SEs
price.km.df.summary <- price.km.list[["100"]] %>%
  
  summarize_at(vars(RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), funs(mean, std.error), na.rm = T) %>%
  
  gather() %>%
  
  mutate(variable = unlist(strsplit(key, "_mean|_std.error")), type = c(rep("mean", 6), rep("se", 6))) %>%
  
  select(-key) %>%
  
  spread(type, value)

# Re-level for better plotting
levels(nopelagic.price.baseline.df.melt$variable) <- c("RICH-L", "RICH-G", "COMP-L", "COMP-G", "CDE", "DIV")

nopelagic.price.baseline.df.melt$variable <- factor(nopelagic.price.baseline.df.melt$variable, levels = c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV"))

nopelagic.price.baseline.df.summary$variable <- c("CDE", "COMP-G", "COMP-L", "RICH-G", "RICH-L", "DIV")

nopelagic.price.baseline.df.summary$variable <- factor(nopelagic.price.baseline.df.summary$variable, levels = c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV"))

price.km.df.summary$variable <- factor(price.km.df.summary$variable, levels = c("RICH_L", "COMP_L", "RICH_G", "COMP_G", "CDE", "TOTAL_DIV"))

levels(price.km.df.summary$variable) <- c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV")

# Generate plots where x-axis is free
nopelagic.price.plot.list <- lapply(levels(nopelagic.price.baseline.df.melt$variable), function(j) {
  
  ylimit <- max(density(subset(nopelagic.price.baseline.df.melt, variable == j)$value)$y) * 0.2
  
  p <- ggplot(subset(nopelagic.price.baseline.df.melt, variable == j)) +
    geom_density(aes(value, fill = variable), col = NA) +
    geom_rug(aes(value), col = "grey80") +
    geom_vline(xintercept = 0, lwd = 0.5) +
    geom_point(
      data = subset(nopelagic.price.baseline.df.summary, variable == j),
      aes(y = ylimit, x = mean), pch = 24, fill = "white", size = 2) +
    geom_point(
      data = subset(price.km.df.summary, variable == j),
      aes(y = ylimit, x = mean), size = 2) +
    # geom_errorbarh(
    #   data = subset(price.km.df.summary, variable == j),
    #   aes(y = ylimit, x = mean, xmin = mean - 2 * se, xmax = mean + 2 * se),
    #   height = 0) +
    coord_flip() +
    scale_fill_manual(values = as.character(subset(colors, variable == j)$col)) +
    scale_x_continuous(limits = c(-1, 1)) +
    labs(y = j, x = "Standardized\nloss of biomass") +
    theme_bw(base_size = 11) + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.ticks = element_blank(), 
      axis.text.x = element_blank(),
      strip.text = element_blank(),
      panel.spacing.x = unit(0, "lines")
    )
  
  if(j != "RICH-L") {
    
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )
    
  }
  
  if(j %in% levels(nopelagic.price.baseline.df.melt$variable)[-1]) {
    
    if(j == "COMP-L")
      
      p <- p + theme(plot.margin = unit(c(1, 1, 1, -0.3), "cm")) else
        
        p <- p + theme(plot.margin = unit(c(1, 1, 1, -1.1), "cm"))
      
  }
  
  # Turn off clipping
  p <- ggplot_gtable(ggplot_build(p))
  p$layout$clip[p$layout$name == "panel"] <- "off"
  
  return(p)
  
} )

pdf(paste0("Figure S7.PDF"),
    width = 7.2, height = 2.5)
bquiet = print(grid.draw(do.call(cbind, c(nopelagic.price.plot.list, size = "last"))))
dev.off() 

#######################################################################################################

# Create vector of site codes ordered by total richness
sitecodes <- survey %>% group_by(SiteCode) %>% summarize(richness = length(unique(SPECIES_NAME))) %>%
  ungroup() %>% arrange(desc(richness)) %>% collect() %>% .[["SiteCode"]]

# Loop over kms
richness.km.mat.list <- lapply(100, function(km) {
  
  # Create vector of site codes to ignore
  ignore <<- c()
  
  # Compute Price components for each site within km radius
  mat.list <- lapply(sitecodes, function(i) {
    
    if(i %in% ignore) NULL else {
      
      # Get latlong of site
      longlat <- unique(subset(survey, SiteCode == i)[, c("SiteLong", "SiteLat")])
      
      # Compute distance between baseline and all other sites
      d <- spDistsN1(longlats, as.numeric(longlat), longlat = T)
      
      # Get vector of comparison sites
      n <- unique(survey$SiteCode)[d <= km]
      
      n <- n[!duplicated(n)]
      
      # Get matrix of biomasses
      bmat <- survey %>% 
        
        filter(SiteCode %in% n) %>% 
        
        filter(!SiteCode %in% ignore) %>%
        
        select(SiteCode, SPECIES_NAME, Biomass) %>%
        
        spread(SPECIES_NAME, Biomass, fill = 0) %>%
        
        as.data.frame()
      
      # Set rownames
      rownames(bmat) <- bmat$SiteCode
      
      bmat <- bmat[, -1, drop = F]
      
      # Get matrix of abundances
      amat <- survey %>% 
        
        filter(SiteCode %in% n) %>% 
        
        filter(!SiteCode %in% ignore) %>%
        
        select(SiteCode, SPECIES_NAME, Abundance) %>%
        
        spread(SPECIES_NAME, Abundance, fill = 0) %>%
        
        as.data.frame()
      
      # Set rownames
      rownames(amat) <- amat$SiteCode
      
      amat <- amat[, -1, drop = F]
      
      # Remove any species with biomass/abundance == 0
      bmat <- bmat[, colSums(bmat) != 0, drop = F]
      
      amat <- amat[, colSums(amat) != 0, drop = F]
      
      # Add SiteCode to ignore
      ignore <<- c(ignore, rownames(amat))
      
      if(nrow(amat) == 0 | nrow(bmat) == 0) NULL else
        
        list(
          amat = amat,
          bmat = bmat
        )
      
    }
    
  } )
  
  mat.list[!sapply(mat.list, is.null)]
  
} )

names(richness.km.mat.list) <- "100"

# Apply Price equation to biomass matrices
richness.price.km.list <- lapply(names(richness.km.mat.list), function(i) {
  
  do.call(rbind, lapply(richness.km.mat.list[[i]], function(j) {
    
    j <- j$bmat
    
    # Remove all matrices with fewer than 5 communities
    if(dim(j)[1] < 5) data.frame() else {
      
      out <- cbind.data.frame(km = i, type = "observed", price2(j, baseline = which.max(rowSums(j > 0)), standardize = TRUE))
      
      out <- out[!out$sharedS == 0, ]
      
      return(out)
      
    }
    
  } ) )
  
} )

names(richness.price.km.list) <- "100"

# # Summarize by baseline site
# richness.price.km.list <- lapply(names(richness.price.km.list), function(i) {
#   
#   richness.price.km.list[[i]] %>% 
#     
#     filter(sharedS != 0) %>% 
#     
#     group_by(km, type, baseline) %>%
#     
#     summarize_at(vars(baselineS, comparisonS, sharedS, RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), 
#                  funs(mean), na.rm = T) %>%
#     
#     ungroup()
#   
# } )
# 
# names(richness.price.km.list) <- "100"

# Get summary data
lapply(names(richness.price.km.list), function(i) {
  
  # Create summary data for means +/- SEs
  richness.price.km.list[[i]] %>%
    
    summarize_at(vars(RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), funs(mean, std.error), na.rm = T) %>%
    
    gather() %>%
    
    mutate(variable = unlist(strsplit(key, "_mean|_std.error")), type = c(rep("mean", 6), rep("se", 6))) %>%
    
    select(-key) %>%
    
    spread(type, value)
  
} )

# Get ratio of RICH to COMP, and DIV to CDE with bootstrapped CI's
lapply(names(richness.price.km.list), function(i) {
  
  dat <- richness.price.km.list[[i]] 
  
  # Get vector of ratios
  RICH_COMP_vec <- with(dat, COMP_L/RICH_L)
  
  DIV_CDE_vec <- with(dat, TOTAL_DIV/CDE)
  
  # Sample 5000 times (with replacement) to generate distribution of means
  RICH_COMP_means <- sapply(1:5000, function(i) {
    
    new_vec <- sample(RICH_COMP_vec, replace = T)
    
    # Take geometric mean
    exp(mean(log(new_vec[!is.infinite(new_vec)]), na.rm = T))
    
  } )
  
  RICH_COMP_means <- RICH_COMP_means[order(RICH_COMP_means)]
  
  # Sample 5000 times (with replacement) to generate distribution of means
  DIV_CDE_means <- sapply(1:5000, function(i) {
    
    new_vec <- sample(DIV_CDE_vec, replace = T)
    
    # Take geometric mean
    exp(mean(log(new_vec[!is.infinite(new_vec)]), na.rm = T))
    
  } )
  
  DIV_CDE_means <- DIV_CDE_means[order(DIV_CDE_means)]
  
  data.frame(
    ratio = c("COMP_RICH", "DIV_CDE"),
    mean = c(exp(mean(log(RICH_COMP_vec[!is.infinite(RICH_COMP_vec)]), na.rm = T)),
             exp(mean(log(DIV_CDE_vec[!is.infinite(DIV_CDE_vec)]), na.rm = T))),
    lower = c(RICH_COMP_means[125], DIV_CDE_means[125]),
    upper = c(RICH_COMP_means[4875], DIV_CDE_means[4875])
  )
  
} )

# Melt data so it is shortways (raw components)
richness.price.km.df.melt <- richness.price.km.list[["100"]] %>% 
  
  gather(variable, value, RICH_L:TOTAL_DIV, factor_key = T)

# Create summary data for means +/- SEs
richness.price.km.df.summary <- richness.price.km.list[["100"]] %>%
  
  filter(sharedS != 0) %>%
  
  summarize_at(vars(RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), funs(mean, std.error), na.rm = T) %>%
  
  gather() %>%
  
  mutate(variable = unlist(strsplit(key, "_mean|_std.error")), type = c(rep("mean", 6), rep("se", 6))) %>%
  
  select(-key) %>%
  
  spread(type, value)

# Re-level for better plotting
richness.price.km.df.melt$variable <- factor(richness.price.km.df.melt$variable, levels = c("RICH_L", "COMP_L", "RICH_G", "COMP_G", "CDE", "TOTAL_DIV"))

levels(richness.price.km.df.melt$variable) <- c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV")

richness.price.km.df.summary$variable <- factor(richness.price.km.df.summary$variable, levels = c("RICH_L", "COMP_L", "RICH_G", "COMP_G", "CDE", "TOTAL_DIV"))

levels(richness.price.km.df.summary$variable) <- c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV")

# Generate plots where x-axis is free
richness.price.km.plot.list <- lapply(levels(richness.price.km.df.melt$variable), function(j) {
  
  ylimit <- max(density(na.omit(subset(richness.price.km.df.melt, variable == j)$value))$y) * 0.2
  
  p <- ggplot(subset(richness.price.km.df.melt, variable == j)) +
    geom_density(aes(value, fill = variable), col = NA) +
    geom_rug(aes(value), col = "grey80") +
    geom_vline(xintercept = 0, lwd = 0.5) +
    geom_point(
      data = subset(richness.price.km.df.summary, variable == j),
      aes(y = ylimit, x = mean), pch = 22, fill = "white", size = 2) +
    geom_point(
      data = subset(price.km.df.summary, variable == j),
      aes(y = ylimit, x = mean), size = 2) +
    # geom_errorbarh(
    #   data = subset(richness.price.km.df.summary, variable == j),
    #   aes(y = ylimit, x = mean, xmin = mean - 2 * se, xmax = mean + 2 * se),
    #   height = 0) +
    coord_flip() +
    scale_fill_manual(values = as.character(subset(colors, variable == j)$col)) +
    scale_x_continuous(limits = c(-1, 1)) +
    labs(y = j, x = "Standardized\nloss of biomass") +
    theme_bw(base_size = 11) + 
    theme(
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      legend.position = "none",
      axis.ticks = element_blank(), 
      axis.text.x = element_blank(),
      strip.text = element_blank(),
      panel.spacing.x = unit(0, "lines")
    )
  
  if(j != "RICH-L") {
    
    p <- p + theme(
      axis.text.y = element_blank(),
      axis.title.y = element_blank()
    )
    
  }
  
  if(j %in% levels(richness.price.km.df.melt$variable)[-1]) {
    
    if(j == "COMP-L")
      
      p <- p + theme(plot.margin = unit(c(1, 1, 1, -0.3), "cm")) else
        
        p <- p + theme(plot.margin = unit(c(1, 1, 1, -1.1), "cm"))
      
  }
  
  # Turn off clipping
  p <- ggplot_gtable(ggplot_build(p))
  p$layout$clip[p$layout$name == "panel"] <- "off"
  
  return(p)
  
} )

pdf("Figure S7.PDF",
    width = 7.2, height = 2.5)
bquiet = print(grid.draw(do.call(cbind, c(richness.price.km.plot.list, size = "last"))))
dev.off() 

#######################################################################################################

# Repeat using original (Fox & Kerr) formulation of the Price equation

# Apply Price equation to biomass matrices
old.price.km.list <- lapply(names(km.mat.list)[4], function(i) {
  
  do.call(rbind, lapply(km.mat.list[[i]], function(j) {
    
    j <- j$bmat
    
    # Remove all matrices with fewer than 5 communities
    if(dim(j)[1] < 5) data.frame() else {
      
      out <- cbind.data.frame(km = i, type = "observed", price(j))
      
      out <- out[out$shared.S != 0, ]
      
      return(out)
      
    }
    
  } ) )
  
} )

names(old.price.km.list) <- "100"

# Get summary data
lapply(names(old.price.km.list), function(i) {
  
  # Create summary data for means +/- SEs
  old.price.km.list[[i]] %>%
    
    summarize_at(vars(RICH_LOSS, RICH_GAIN, COMP_LOSS, COMP_GAIN, CDE, TOTAL_DIV), lst(mean, std.error), na.rm = T) %>%
    
    gather() %>%
    
    mutate(variable = unlist(strsplit(key, "_mean|_std.error")), type = c(rep("mean", 6), rep("se", 6))) %>%
    
    select(-key) %>%
    
    spread(type, value)
  
} )

# Get ratio of RICH to COMP, and DIV to CDE with bootstrapped CI's
lapply(names(old.price.km.list), function(i) {
  
  dat <- old.price.km.list[[i]] 
  
  # Get vector of ratios
  RICH_COMP_vec <- with(dat, COMP_LOSS/RICH_LOSS)
  
  DIV_CDE_vec <- with(dat, TOTAL_DIV/CDE)
  
  # Sample 5000 times (with replacement) to generate distribution of means
  RICH_COMP_means <- sapply(1:5000, function(i) {
    
    new_vec <- sample(RICH_COMP_vec, replace = T)
    
    # Take geometric mean
    exp(mean(log(new_vec[!is.infinite(new_vec)]), na.rm = T))
    
  } )
  
  RICH_COMP_means <- RICH_COMP_means[order(RICH_COMP_means)]
  
  # Sample 5000 times (with replacement) to generate distribution of means
  DIV_CDE_means <- sapply(1:5000, function(i) {
    
    new_vec <- sample(DIV_CDE_vec, replace = T)
    
    # Take geometric mean
    exp(mean(log(new_vec[!is.infinite(new_vec)]), na.rm = T))
    
  } )
  
  DIV_CDE_means <- DIV_CDE_means[order(DIV_CDE_means)]
  
  data.frame(
    ratio = c("COMP_RICH", "DIV_CDE"),
    mean = c(exp(mean(log(RICH_COMP_vec[!is.infinite(RICH_COMP_vec)]), na.rm = T)),
             exp(mean(log(DIV_CDE_vec[!is.infinite(DIV_CDE_vec)]), na.rm = T))),
    lower = c(RICH_COMP_means[125], DIV_CDE_means[125]),
    upper = c(RICH_COMP_means[4875], DIV_CDE_means[4875])
  )
  
} )

# Plot results
lapply(as.character("100"), function(i) { # names(price.km.list)
  
  # Melt data so it is shortways (raw components)
  price.km.df.melt <- old.price.km.list[[i]] %>% 
    
    gather(variable, value, RICH_L:TOTAL_DIV, factor_key = T)
  
  # Create summary data for means +/- SEs
  price.km.df.summary <- old.price.km.list[[i]] %>%
    
    summarize_at(vars(RICH_L, RICH_G, COMP_L, COMP_G, CDE, TOTAL_DIV), funs(mean, std.error), na.rm = T) %>%
    
    gather() %>%
    
    mutate(variable = unlist(strsplit(key, "_mean|_std.error")), type = c(rep("mean", 6), rep("se", 6))) %>%
    
    select(-key) %>%
    
    spread(type, value)
  
  # Re-level for better plotting
  price.km.df.melt$variable <- factor(price.km.df.melt$variable, levels = c("RICH_L", "COMP_L", "RICH_G", "COMP_G", "CDE", "TOTAL_DIV"))
  
  levels(price.km.df.melt$variable) <- c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV")
  
  price.km.df.summary$variable <- factor(price.km.df.summary$variable, levels = c("RICH_L", "COMP_L", "RICH_G", "COMP_G", "CDE", "TOTAL_DIV"))
  
  levels(price.km.df.summary$variable) <- c("RICH-L", "COMP-L", "RICH-G", "COMP-G", "CDE", "DIV")
  
  # Generate plots where x-axis is free
  price.km.plot.list <- lapply(levels(price.km.df.melt$variable), function(j) {
    
    ylimit <- max(density(na.omit(subset(price.km.df.melt, variable == j)$value))$y) * 0.2
    
    p <- ggplot(subset(price.km.df.melt, variable == j)) +
      geom_density(aes(value, fill = variable), col = NA) +
      geom_rug(aes(value), col = "grey80") +
      geom_vline(xintercept = 0, lwd = 0.5) +
      geom_point(
        data = subset(price.km.df.summary, variable == j),
        aes(y = ylimit, x = mean), size = 2) +
      geom_errorbarh(
        data = subset(price.km.df.summary, variable == j),
        aes(y = ylimit, x = mean, xmin = mean - 2 * se, xmax = mean + 2 * se),
        height = 0) +
      coord_flip() +
      scale_fill_manual(values = as.character(subset(colors, variable == j)$col)) +
      scale_x_continuous(limits = c(-1, 1)) +
      labs(y = j, x = "Standardized\nloss of biomass") +
      theme_bw(base_size = 11) + 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text.x = element_blank(),
        strip.text = element_blank(),
        panel.spacing.x = unit(0, "lines")
      )
    
    if(j != "RICH-L") {
      
      p <- p + theme(
        axis.text.y = element_blank(),
        axis.title.y = element_blank()
      )
      
    }
    
    if(j %in% levels(price.km.df.melt$variable)[-1]) {
      
      if(j == "COMP-L")
        
        p <- p + theme(plot.margin = unit(c(1, 1, 1, -0.3), "cm")) else
          
          p <- p + theme(plot.margin = unit(c(1, 1, 1, -1.1), "cm"))
        
    }
    
    # Turn off clipping
    p <- ggplot_gtable(ggplot_build(p))
    p$layout$clip[p$layout$name == "panel"] <- "off"
    
    return(p)
    
  } )
  
  pdf(paste0("Figure 1_", i, "km-new.PDF"),
      width = 7.2, height = 2.5)
  bquiet = print(grid.draw(do.call(cbind, c(price.km.plot.list, size = "last"))))
  dev.off() 
  
} )

dev.off()