load('/Users/Mark/Downloads/bci.tree8.rdata')
a = bci.tree8

range(a$gx, na.rm=T) # 0 to 1000 is the range of x (10 ha)
range(a$gy, na.rm=T) # 0 to 500 is the range of y (5 ha)

x.seq = seq(from=0, to=1000, length.out=11) # 10 divisions on x, length 11 because we need endpoint
y.seq = seq(from=0, to= 500, length.out= 6) #  5 divisions on y, length  6 because we need endpoint

counter = 0 # this will help us combine unique i and j into one variable for the output list
bci.50 = vector('list', 50) # empty list, one for each of the 1-ha plots we want to make

for(i in 1:10){ # x axis
  for(j in 1:5){ # y axis
    counter = counter+1
    # "clunky but it works" way to select the right area
    bci.50[[counter]] = a[which(a$gx >= x.seq[i] & a$gx < x.seq[i+1] &
                                  a$gy >= y.seq[j] & a$gy < y.seq[j+1]),]
    bci.50[[counter]]$subplot = counter
  }
}

# checking the number of rows in each of the 50 new  1-ha plots, looks reasonable
plot(unlist(lapply(bci.50, nrow)), ylim=c(0,12000))

# checking what the r-binded df looks like, looks good
bci_50subplots = do.call(rbind, bci.50)
str(bci_50subplots)

# but need to make sure all the unique subplots showed up
unique(bci_50subplots$subplot)

# looks good
