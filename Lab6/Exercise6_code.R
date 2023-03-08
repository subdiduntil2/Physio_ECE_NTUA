### Exercise 6 - Practice part of analysis of microarray experiment data


### ZITOUMENO 2
 
setwd("C:/Users/manos/Lab06") # TODO complete by adding your working directory, with format e.g. C:/Users/.../

# 2-i

x <-read.table(file.path("GSE7117_series_matrix.txt"), skip=68, header=TRUE, sep="\t",row.names=1) # explain in your assignment the arguments used in this function

# TODO 2-v to 2-x



### ZITOUMENO 3

# 3-i 

round(apply(x,2, summary),digits=2)  # explain in your assignment the functionality/arguments in this line of code


# 3-ii 

# Boxplots to check the normalisation of data
boxplot(as.data.frame(x),col=c("red","green","blue","black","pink","yellow","orange","brown"))

# TODO write a function to plot boxplots for every sample of the microarray experiment described in x. Use a different colour for the two experimental conditions. Choose not to plot the outliers. 
boxplot(x, col=c("red","blue", "blue", "red", "blue", "red", "blue", "red"), outlier.colour=NA)
boxplot(x, col=c("yellow","red", "red", "yellow", "red", "yellow", "red", "yellow"), outlier.colour=NA,log="y")

### ZITOUMENO 4 

# perform t-test to estimate the statistical significance of genes' differential expression

# 4-i

xsamplelabels<-c(0,1,1,0,1,0,1,0) # TODO complete function c() according to the labels of samples (0 for control samples and 1 for diet intervention samples)

# 4-ii

# TODO perform t-test using function t.test() so that the p-value is calculated for every probe, for the two experimental conditions (control/diet intervention)
p_values <-  c()
c_control=c(x[,1],x[,4],x[,6],x[,8])
c_interv=c(x[,2],x[,3],x[,5],x[,7])
for (i in 1:54675){
  print(i)
  c_control=c(x[i,1],x[i,4],x[i,6],x[i,8])
  c_interv=c(x[i,2],x[i,3],x[i,5],x[i,7])
  temp <- t.test(c_control,c_interv)
  p_values[i] <- temp$p.value
}

# 4-iii

# TODO write code to find which probes have p-value less than 0.001
sign_values <- as.data.frame(p_values[p_values<0.001])
sign_indx <- (which(p_values<0.001))
sign_name <- c()
for (i in 1:38){
  temp_i<-sign_indx[i] 
  print(temp_i)
  sign_name[i]<-rownames(x)[temp_i]
}
sign_name <- as.data.frame(sign_name)
sign_total <- cbind(sign_indx,sign_name,sign_values)
# TODO write code to find the probe with the smallest p-value
min_indx=as.data.frame(which.min(p_values))
print(min_indx[1,1])
min_name=as.data.frame(rownames(x)[(min_indx[1,1])])
min_value=as.data.frame(p_values[(min_indx[1,1])])
min_total <- cbind(min_indx,min_name,min_value)






