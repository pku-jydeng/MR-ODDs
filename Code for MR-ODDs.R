library(splines)
library(mgcv)
library(metafor)
library(ggplot2)
library(mixmeta)
library(tidyverse)
library(rootSolve)
library(meta)
rm(list = ls());gc()
dta <- readxl::read_xlsx("data.xlsx") #数据示例
dta <- filter(dta, dta$`Sex` =="All"|dta$`Sex`=="Both")

# 1. Clean data -----------------------------------------------------------------
dta$`Effect Indicator` <- replace(dta$`Effect Indicator`, dta$`Effect Indicator` == "IR", "ER")
dta$`Effect Indicator` <- replace(dta$`Effect Indicator`, dta$`Effect Indicator` == "Perc-Incr", "ER")
dta$`Effect Indicator` <- replace(dta$`Effect Indicator`, dta$`Effect Indicator` == "HR", "RR")
dta$`Unit of Measurement` <- replace(dta$`Unit of Measurement`, dta$`Unit of Measurement` == "mcg/m3", "μg/m3")

dta$`Associated Increment` <- gsub("μg/m3","",dta$`Associated Increment`)
dta$`Associated Increment` <- gsub("mcg/m3","",dta$`Associated Increment`)
dta$`Associated Increment` <- gsub("ppb","",dta$`Associated Increment`)
dta$`Associated Increment` <- as.numeric(dta$`Associated Increment`)
dta$`Associated Increment`[which(dta$`Unit of Measurement` == "ppb")]<- 47.997/24.47*dta$`Associated Increment`[which(dta$`Unit of Measurement` == "ppb")]
dta$Mean<- replace(dta$Mean,dta$Mean == "NA", NA)
dta$Mean<- replace(dta$Mean,dta$Mean == "", NA)
dta$Mean<- as.numeric(dta$Mean)
dta$median <- replace(dta$median,dta$median == "NA", NA)
dta$median<- replace(dta$median,dta$median == "", NA)
dta$median<- as.numeric(dta$median)

dta$`Unit of Measurement` <- replace(dta$`Unit of Measurement`, dta$`Unit of Measurement` == "ppb", "μg/m3")

dta$`Single pollutant - effect estimate` <- as.numeric(dta$`Single pollutant - effect estimate`)
dta$`Single pollutant - higher interval limit (95% CI)` <- as.numeric(dta$`Single pollutant - higher interval limit (95% CI)`)
dta$`Single pollutant - lower interval limit (95% CI)` <- as.numeric(dta$`Single pollutant - lower interval limit (95% CI)`)

dta <- dta %>% mutate(beta = if_else(`Effect Indicator` == "ER",
                                     log(`Single pollutant - effect estimate`/100+1)/`Associated Increment`,
                                     log(`Single pollutant - effect estimate`)/`Associated Increment`))
dta <- dta %>% mutate(beta_low = if_else(`Effect Indicator`== "ER",
                                         log(`Single pollutant - lower interval limit (95% CI)`/100+1)/`Associated Increment`,
                                         log(`Single pollutant - lower interval limit (95% CI)`)/`Associated Increment`))
dta <- dta %>% mutate(beta_up = if_else(`Effect Indicator`== "ER",
                                        log(`Single pollutant - higher interval limit (95% CI)`/100+1)/`Associated Increment`,
                                        log(`Single pollutant - higher interval limit (95% CI)`)/`Associated Increment`))
dta <- dta %>% mutate(se = (beta_up - beta_low)/(2*1.96))

dta_p <- filter(dta, Pollutant == "PM2.5", !is.na(beta))

# 2. Reconstruct the exposure distribution based on the gamma distribution ------------------------------------------------------
dta_newp <- dta_p
dta_newp <- dta_p[which(!is.na(dta_p$median)|(!is.na(dta_p$Mean))),]
dta_newp$Mean <-as.numeric(dta_newp$Mean)
dta_newp$median <-as.numeric(dta_newp$median)
dta_newp$`Interquartile range` <- as.numeric(dta_newp$`Interquartile range`)
dta_newp$`Standard deviation` <- as.numeric(dta_newp$`Standard deviation`)
dta_newp$P5 <- as.numeric(sub("-.*", "", dta_newp$`5-95th percentiles`))
dta_newp$P10 <- as.numeric(sub("-.*", "", dta_newp$`10-90th percentiles`))
dta_newp$P90 <- as.numeric(sub(".*-", "", dta_newp$`10-90th percentiles`))
dta_newp$P95 <- as.numeric(sub(".*-", "", dta_newp$`5-95th percentiles`))
dta_newp$P25 <- as.numeric(dta_newp$P25)
dta_newp$P75 <- as.numeric(dta_newp$P75)

library(nloptr)
# The function can contain most of the exposed information of the record, add penalty items, and have different weights for different statistics
estimate_gamma <- function(
    mean_val = NA, sd_val = NA, median_val = NA, iqr_val = NA,
    q5 = NA, q95 = NA, q10 = NA, q90 = NA, q25 = NA, q75 = NA
) {
  # Dynamically constructed based on available statistics
  objective <- function(params) {
    k <- params[1]
    theta <- params[2]
    
    # Calculate the theoretical statistics
    theor_mean <- k * theta
    theor_sd <- sqrt(k) * theta
    theor_median <- qgamma(0.5, shape = k, scale = theta)
    theor_iqr <- qgamma(0.75, shape = k, scale = theta) - qgamma(0.25, shape = k, scale = theta)
    
    # Calculate the theoretical quantile
    if (!is.na(q5)) theor_q5 <- qgamma(0.05, shape = k, scale = theta)
    if (!is.na(q95)) theor_q95 <- qgamma(0.95, shape = k, scale = theta)
    if (!is.na(q10)) theor_q10 <- qgamma(0.10, shape = k, scale = theta)
    if (!is.na(q90)) theor_q90 <- qgamma(0.90, shape = k, scale = theta)
    if (!is.na(q25)) theor_q25 <- qgamma(0.25, shape = k, scale = theta)
    if (!is.na(q75)) theor_q75 <- qgamma(0.75, shape = k, scale = theta)
    
    # Dynamic calculation error
    error <- 0
    if (!is.na(mean_val)) error <- error + (theor_mean - mean_val)^2*5 
    if (!is.na(sd_val)) error <- error + (theor_sd - sd_val)^2*5  
    if (!is.na(median_val)) error <- error + (theor_median - median_val)^2
    if (!is.na(iqr_val)) error <- error + (theor_iqr - iqr_val)^2  *2
    if (!is.na(q5)) error <- error + (theor_q5 - q5)^2 *0.5 
    if (!is.na(q95)) error <- error + (theor_q95 - q95)^2 *0.5
    if (!is.na(q10)) error <- error + (theor_q10 - q10)^2 *0.5
    if (!is.na(q90)) error <- error + (theor_q90 - q90)^2 *0.5
    if (!is.na(q25)) error <- error + (theor_q25 - q25)^2 *0.5
    if (!is.na(q75)) error <- error + (theor_q75 - q75)^2 *0.5
    
    return(error)
  }
  
  # Dynamic calculation error
  if (!is.na(mean_val) && !is.na(sd_val)) {
    k_init <- (mean_val / sd_val)^2
    theta_init <- sd_val^2 / mean_val
  } else if (!is.na(median_val) && !is.na(iqr_val)) {
    k_init <- (median_val / iqr_val)^2 * 1.5  # Empirical formula
    theta_init <- iqr_val / (qgamma(0.75, k_init) - qgamma(0.25, k_init))
  } else if(!is.na(mean_val) && !is.na(iqr_val)){
    k_init <- (mean_val / iqr_val)^2
    theta_init <- iqr_val / (qgamma(0.75, k_init) - qgamma(0.25, k_init))
  } else if(!is.na(mean_val)){
    k_init <- mean_val 
    theta_init <- 1
  } else {
    k_init <- median_val
    theta_init <- 1
  }
  
  # Constrained optimization
  result <- nloptr(
    x0 = c(k_init, theta_init),
    eval_f = objective,
    lb = c(1e-6, 1e-6),
    ub = c(Inf, Inf),
    opts = list(
      algorithm = "NLOPT_LN_COBYLA",
      xtol_rel = 1e-8,
      maxeval = 2000
    )
  )
  
  return(list(
    k = result$solution[1],
    theta = result$solution[2],
    goodness = objective(result$solution)
  ))
}

tmp_k<-numeric()
tmp_theta <-numeric()
for (i in 1:nrow(dta_newp)){
  tmp<-estimate_gamma(
    mean_val = dta_newp$Mean[i], 
    sd_val = dta_newp$`Standard deviation`[i], 
    median_val = dta_newp$median[i], 
    iqr_val = dta_newp$`Interquartile range`[i],
    q5 = dta_newp$P5[i], 
    q95 = dta_newp$P95[i], 
    q10 = dta_newp$P10[i], 
    q90 = dta_newp$P90[i], 
    q25 = dta_newp$P25[i], 
    q75 = dta_newp$P75[i]
  )
  tmp_k<-cbind(tmp_k,tmp$k)
  tmp_theta<-cbind(tmp_theta,tmp$theta)
}

dta_newp$k <- as.vector(tmp_k)
dta_newp$theta <-as.vector(tmp_theta)

# 3. Calculate the weights and the first-order derivative -----------------------------------------------------------
m1 <- rma.mv(
  yi = beta,           
  V = se^2,            
  random = ~ 1 | Study_ID,
  data = dta_newp,
  method = "ML" 
)
weights_raw <- weights(m1)
weights_normalized <- weights_raw / sum(weights_raw)

n_samples <- 100000
sample_counts <- round(weights_normalized * n_samples)

bs_sample <- vector()
for (i in 1:nrow(dta_newp)){
  shape <- dta_newp$k[i]
  scale <- dta_newp$theta[i]
  bs_sample <- c(bs_sample,rgamma(sample_counts[i], shape = shape, scale = scale))
}

tmp <-data.frame(x=bs_sample,y=1)

ndf=3
ndigits=10
dta_newp$x <- dta_newp$Mean
tmp_coef_new <- matrix(nrow = 0, ncol = (6+ndf))
bm <- gam(y~s(x, k=ndf), data = tmp) 

tmp_prop <- numeric()
P2.5<-numeric()
P97.5<-numeric()
for (i in 1:nrow(dta_newp)){
  # Set the range of integration
  max_x <-ceiling(qgamma(0.9999, shape = dta_newp$k[i] , scale = dta_newp$theta[i]))
  max_x <- max(max_x,200)
  max_x <- min(max_x,1000)# Set the range of integration
  # Obtain the exposed discrete distribution
  tmp_weight <- data.frame(x=seq(0,max_x,1/ndigits),weight=NA)
  tmp_p <- dgamma(tmp_weight$x, shape = dta_newp$k[i] , scale = dta_newp$theta[i])/ndigits
  tmp_p[1]<-0 
  P2.5[i] <- qgamma(0.025, shape = dta_newp$k[i] , scale = dta_newp$theta[i])
  P97.5[i] <- qgamma(0.975, shape = dta_newp$k[i] , scale = dta_newp$theta[i])
  
  bs_co <-  predict(bm, newdata = data.frame(x=seq(0,max_x,1/ndigits)), type ="lpmatrix")
  for (nrow in 1:nrow(tmp_weight)) 
  {
    tmp_weight$weight[nrow] <- (sum(tmp_weight$x[nrow:nrow(tmp_weight)]*tmp_p[nrow:nrow(tmp_weight)]/sum(tmp_p[nrow:nrow(tmp_weight)]),na.rm=T)-dta_newp$k[i]*dta_newp$theta[i])*sum(tmp_p[nrow:nrow(tmp_weight)])/(dta_newp$k[i]*dta_newp$theta[i]^2)/ndigits
  }
  tmp_weight$weight[is.na(tmp_weight$weight)]<-0
  
  tmp_w <- numeric(ndf)
  for (n in 1:ndf) {
    tmp_w[n] <- sum(tmp_weight$weight * bs_co[, n])
  }
  tmp_coef_new <- rbind(tmp_coef_new, cbind(dta_newp$Mean[i], dta_newp$k[i]*dta_newp$theta[i], dta_newp$beta[i], dta_newp$se[i],P2.5[i] ,P97.5[i],
                                            t(tmp_w)))
  if (sum(tmp_weight$weight)<0.99){
    tmp_prop <- cbind(tmp_prop ,i)
  }
  print(i)
}

# 4. Obtain first-order derivatives by meta-regression----------------------------------------------------------
dta_reg<-as.data.frame(tmp_coef_new) 
colnames(dta_reg)<-c("x","x_sim","beta","se","P2.5","P97.5",paste0("bs_gamma",1:ndf))
bs_gamma<-as.matrix(dta_reg[,7:(6+ndf)])
bs <-  predict(bm, newdata = data.frame(x=dta_newp$x[which(!is.na(dta_newp$Mean))]), type ="lpmatrix")
colnames(bs) <- paste0("bs",1:ndf)
dta_newp$bs_gamma1 <- bs_gamma[,1]
dta_newp$bs_gamma2 <- bs_gamma[,2]
dta_newp$bs_gamma3 <- bs_gamma[,3]


mr_gamma <- rma.mv(
  yi = beta,
  V = se^2,
  mods = ~ bs_gamma - 1,   
  random = ~ 1| Study_ID, 
  data = dta_newp,  
  method = "ML"
)


ERF<-data.frame(x=seq(0,300,1/ndigits))
bs <- predict(bm, newdata = ERF, type ="lpmatrix")

# Obtain the first-order derivative and the confidence interval
ERF$fit_gamma=bs%*%coef(mr_gamma)
ERF$se_gamma=sqrt(diag(bs%*%vcov(mr_gamma)%*%t(bs)))
scl<-max(ERF$fit_gamma)-min(ERF$fit_gamma)


ggplot()+
  geom_histogram(data=dta_reg,aes(x=x,y=(..ncount..)*scl*10),fill="grey80")+
  geom_ribbon(data=ERF,aes(x=x,ymin=(fit_gamma-se_gamma*1.96)*10,ymax=(fit_gamma+se_gamma*1.96)*10),alpha=0.5)+
  geom_path(data=ERF,aes(x=x,y=fit_gamma*10))+
  theme_bw()+
  xlab(expression("PM"[2.5]))+ylab(expression("ln(RR) for 10" ~mu~g/m^3~"increase of PM"[2.5]))+
  ggtitle("First-order derive of exposure-response function")

# 5. The exposure response function is obtained by integration. --------------------------------------------------------------
library(MASS)
coef_sim_gamma<-mvrnorm(500,coef(mr_gamma),vcov(mr_gamma))
coef_sim_gamma[1,]=coef(mr_gamma)
id1<-which.min(abs(ERF$x-15))
tmp1<-apply(bs[-(1:id1),]%*%t(coef_sim_gamma)*1/ndigits,2,cumsum)
tmp2<-apply(bs[c((id1-1):1),]%*%t(coef_sim_gamma)*-1/ndigits,2,cumsum)[c((id1-1):1),]
ERF_sim_gamma<-rbind(tmp2,0,tmp1)

ERF$ERF_fit_gamma=ERF_sim_gamma[,1]
ERF$ERF_lo_gamma=apply(ERF_sim_gamma,1,quantile,0.025)
ERF$ERF_up_gamma=apply(ERF_sim_gamma,1,quantile,1-0.025)
scl1<-max(ERF$ERF_fit_gamma)-min(ERF$ERF_fit_gamma)

coef<-coef(mr_gamma)
vc<-vcov(mr_gamma)
save(ERF,coef,vc,mr_gamma,bm,bs,file="final_de.RData")

left_y_range <- c(0.92, 1.16)
left_y_breaks <- seq(0.92, 1.16, by = 0.02)
right_y_range <- c(0.8, 1.5)
right_y_breaks <- seq(0.8, 1.5, by = 0.06)
scale_factor <- diff(left_y_range)/diff(right_y_range)
library(mgcv)

dta_samples <- data.frame(x = bs_sample)
dta_samples_binned <- dta_samples %>%
  mutate(x_bin = cut(x, breaks = seq(0, 300, by = 1), include.lowest = TRUE)) %>%  
  group_by(x_bin) %>%
  summarize(
    x = mean(x, na.rm = TRUE),  
    count = n(),                
    .groups = "drop"
  ) %>%
  mutate(alpha = scales::rescale(count, to = c(0, 1)))  

dens <- density(dta_samples$x, from = 0)
dual_plot <- ggplot() +
  geom_line(
    data = data.frame(x = dens$x, y = dens$y*1+0.96),
    aes(x = x, y = y,color = "The density function of PM2.5 distribution"), 
    linewidth = 0.5,
    alpha = 0.6)+
  geom_ribbon(data = data.frame(x = dens$x, y = dens$y*1+0.96)
              , aes(x = x, ymin = 0.96, ymax = y),
              alpha = 0.30, fill = "gray") +
  geom_point(data = dta_reg,
             aes(x = x, y = exp(beta * 10) ,color = "RR for 10 μg/m³ increase of PM2.5"),
             size = 0.1,
             alpha = 0.5) +
  geom_errorbar(data = dta_reg,
                aes(x = x,
                    ymin = exp((beta - 1.96 * se) * 10),
                    ymax = exp((beta + 1.96 * se) * 10),
                    color = "RR for 10 μg/m³ increase of PM2.5"),
                linewidth = 0.05, alpha = 0.20) +
  geom_segment(
    data = dta_reg,
    aes(y = exp(beta*10),color =  "PM2.5 95% CI", yend = exp(beta*10),
        x = P2.5 , xend = P97.5),
    linewidth = 0.05,
    alpha = 0.20) +
  geom_ribbon(data = ERF, aes(x = x, 
                              ymin = exp((fit_gamma - se_gamma*1.96)*10), 
                              ymax = exp((fit_gamma + se_gamma*1.96)*10)),
              fill = "navy", alpha = 0.15) +
  geom_path(data = ERF, aes(x = x, y = exp(fit_gamma*10),
                            color = "Marginal effect"), 
            linewidth = 0.8,alpha = 0.7) +
  geom_ribbon(data = ERF, aes(x = x, 
                              ymin = (exp(ERF_lo_gamma)-1)*scale_factor+1, 
                              ymax = (exp(ERF_up_gamma)-1)*scale_factor+1),
              alpha = 0.10, fill = "maroon") +
  geom_path(data = ERF, aes(x = x, y = (exp(ERF_fit_gamma)-1)*scale_factor+1,
                            color = "Exposure-response function"),
            linewidth = 0.8,alpha = 0.7) +
  scale_y_continuous(
    name = expression("RR for 10" ~"μg/m"^3~"increase of PM"[2.5]),
    breaks = left_y_breaks,
    limits = left_y_range,
    expand = c(0, 0),
    sec.axis = sec_axis(
      trans = ~ ((. - 1) / scale_factor+1),
      name = expression("Risk Ratio (RR)"),
      breaks = right_y_breaks)) +
  scale_x_continuous(
    breaks = c(0, 15, 100, 75, 200,300, 50, 37.5,25),
    labels = c("0", "15", "100", "75", "200","300", "50","37.5","25"),
    sec.axis = sec_axis(
      trans = ~ .,  
      name = NULL,  
      breaks = c(0, 15, 100,75,200,300,50,37.5,25),  
      labels = c("", "AQG", "","IT-1","","","IT-2","IT-3","IT-4"))   
  ) +
  labs(x = expression("Average"~"PM"[2.5]~"exposure of centers"~"("~"μg/m"^3~")")) +
  # geom_vline(xintercept = 100, linetype = "dotted", color = "black", linewidth = 0.5) + 
  geom_vline(xintercept = 75, linetype = "dotted", color = "black", linewidth = 0.5) + 
  geom_vline(xintercept = 15, linetype = "dashed", color = "black", linewidth = 0.5) + 
  geom_vline(xintercept = 50, linetype = "dotted", color = "black", linewidth = 0.5) + 
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.5) +
  theme_bw() +  
  theme(
    panel.grid = element_blank(), 
    axis.title.y.right = element_text(color = "maroon"),
    axis.title.y.left  = element_text(color = "navy"),
    axis.text.y.left = element_text(color = "navy"),
    axis.text.y.right = element_text(color = "maroon"),
    axis.text.x.top = element_text(size=7), 
    axis.ticks.x.top = element_line(),
    legend.position = c(0.90, 0.05),               
    legend.justification = c(1, 0.05),          
    legend.box.margin = margin(t = -10, r = -10),
    legend.background = element_rect(
      fill = alpha("white", 0),              
      color = NA)
  ) +
  ggtitle(expression("PM"[2.5]~" exposure and all cause mortality")) +
  scale_color_manual(name = NULL,values = c(
    "Exposure-response function" = "maroon",
    "Marginal effect" = "navy",
    "PM2.5 95% CI" = "gray10",
    "RR for 10 μg/m³ increase of PM2.5" = "navy",
    "The density function of PM2.5 distribution" = "grey5"),
    labels = c(
      "Exposure-response function",
      "Marginal effect", 
      expression("PM"[2.5]~"95% CI"),
      expression("RR for 10" ~"μg/m"^3~"increase of PM"[2.5]),
      expression("The density function of PM"[2.5]~"distribution")
    ))+
  # scale_alpha_identity() +
  coord_cartesian(ylim = c(0.96, 1.05),xlim= c(0,300))
ggsave("combined_plotall.png", plot = dual_plot, width = 8, height = 6, dpi = 300)

