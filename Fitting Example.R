library(lubridate) 
library(sf)
library(spacetime)
library(tidyverse)
library(crimedata)
library(zoo)
library(INLA)
library(INLAspacetime)
library(inlabru)
library(stelfi)
library(spatstat)


test2 <- homicides15 %>%
  filter(city_name=="Chicago")

df_time=data.frame(time=test2$date_single)


weekly_data <- df_time %>%
  mutate(week = cut(time,breaks = "week")) %>%
  group_by(week) %>%
  summarize(tot=n())

data_weekly <- weekly_data

acf_result <- acf(data_weekly$tot, plot = FALSE)

# Create a data frame with ACF values and lag
acf_data <- data.frame(lag = acf_result$lag, acf = acf_result$acf)

# Plot the ACF using ggplot2
ggplot(acf_data, aes(x = lag, y = acf)) +
  geom_segment(aes(xend = lag, yend = 0), color = "blue", size = 1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_point(color = "red", size = 3) +
  labs(title = "Autocorrelation Function (ACF)",
       x = "Lag",
       y = "ACF",
       caption = "Autocorrelation of weekly values") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_text(hjust = 0.5))



# Aggregate data by week
data_weekly <-  df_time %>%
  mutate(week = as.Date(cut(time, breaks = "week"))) %>%
  group_by(week) %>%
  summarise(tot = n())



# Compute a 3-week rolling average
data_weekly <- data_weekly %>%
  mutate(rolling_avg = rollmean(tot, 3, fill = NA, align = "right"))


# Plot the weekly values and rolling average
ggplot(data_weekly, aes(x = week)) +
  geom_line(aes(y = tot), color = "blue", size = 1, linetype = "dashed") +
  geom_line(aes(y = rolling_avg), color = "red", size = 1.5) +
  labs(title = "3-Week Rolling Average",
       x = "Week",
       y = "Value",
       caption = "Blue dashed line: Weekly values, Red line: 3-week rolling average") +
  theme_minimal(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.caption = element_text(hjust = 0.5))



time <- test2$date_single

time<-as.numeric(difftime(time,min(na.omit(time)),units="days"))

x<- test2$latitude-min(test2$latitude)

y<- test2$longitude-min(test2$longitude)

df<-data.frame(x = x,y=y,t=time)


#df <- df %>% mutate(t=(t-min(t))/(max(t)-min(t)))

region<-data.frame(x=c(min(x),max(x)),
                   y=c(min(y),max(y)),
                   t=c(min(df2$t),max(df2$t)))



##params<-data.frame(beta0=.5,theta=.3,lambda=.1,sigma2=.2)


#hawkes_est(df2,params,region)




#Remove Duplicates
df2 <- df %>% distinct()


#Spatially Visualize Data
df2 %>% ggplot(aes(x=x,y=y)) + 
  geom_point() + 
  theme_bw() + ggtitle("Homicide Locations in Chicago (2015)")
  



points_spat <- sf::st_as_sf(df2,coords=c("x","y"))


#Create Boundary
asdf <- points_spat %>% 
  summarise() %>% 
  concaveman::concaveman(concavity = 10)

#Create INLA Mesh
bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(asdf)[, 1:2]))
#smesh <- INLA::inla.mesh.2d(boundary = bnd,
#                            max.edge = 0.75, cutoff = 0.3)


domain=asdf
smesh <- fmesher::fm_mesh_2d_inla(
  loc = cbind(df$x,df$y), 
  max.edge = .2, 
  offset = .05)

plot(smesh)
#lines(asdf,lwd=2,col="green4")
#Visualize Points on MESH
points(cbind(df$x,df$y),pch = 19, col="blue")



df_spat <- data.frame(x=df$x,y=df$y)

w <- as.owin(list(xrange=c(0,0.365), yrange=c(0,.298814)))

df_ppp <- as.ppp(df_spat,W=w)

#Calculate Border corrected K function
K_est <- Kest(df_ppp, correction = "border")

#Calculate LGCP K function fit to empirical K function
LGCP<- lgcp.estK(K_est)
plot(LGCP)

#Calculate Matern Cluster Process K function fit to empirical K function
Mat_Clus <- matclust.estK(K_est)
plot(Mat_Clus)


locs <- df_spat


#Fit Spatial only LGCP to data (Not shown in Paper)
fit <- fit_lgcp(locs = locs, sf = domain, smesh = smesh,
                parameters = c(beta = 0, log_tau = log(1),
                               log_kappa = log(1)))

show_lambda(fit, smesh = smesh, sf = domain) + ggplot2::theme_void()


#bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[,1:2]))

#smesh <- INLA::inla.mesh.2d(boundary = bnd)


#locs_st <- df
df<- df %>% mutate(t=t+.5)

df<-df %>% 
  mutate(tmod=t/max(t))


df2 <- data.frame(x=df$x,y=df$y,t=df$t)

#bnd <- INLA::inla.mesh.segment(as.matrix(sf::st_coordinates(domain)[, 1:2]))

#smesh <- INLA::inla.mesh.2d(boundary = bnd,
#                            max.edge = 0.75, cutoff = 0.3)

w0 <- 30

tmesh <- fm_mesh_1d(seq(0,max(df$t),by=w0))

df2 <- df2[-499,]


#Fitting Spatio-Temporal LGCP
fit <- fit_lgcp(locs = df2, sf = domain, smesh = smesh, tmesh = tmesh,
                parameters = c(beta = 0, log_tau = log(1),
                               log_kappa = log(1), atanh_rho = 0.2))



get_coefs(fit)


df_stelfi <- df2 %>% unique()

df_stelfi <- df_stelfi[-240,]



locs = data.frame(x=df_stelfi$x,y=df_stelfi$y)


param <- list(mu = 1, alpha = 20, beta = 200, kappa = 2, tau = 0.1, 
              xsigma = 0.2, ysigma = 0.2, rho = 0)

#If receive: Error in log(parameters[["tau"]]) : 
#non-numeric argument to mathematical function, try rerunning

Hawkes_fit <- fit_stelfi(times=df_stelfi$t, locs=locs, 
                         sf=domain,smesh=smesh,parameters=param,
                         GMRF = T)


get_coefs(Hawkes_fit)

show_hawkes(list(times = df_stelfi$t, params = c(mu = 10.284218638, 
                                                 alpha = 0.213540795, 
                                                 beta = 19.825459470)))

