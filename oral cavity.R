if (!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, survival, ggfortify, survminer, plotly, gridExtra, 
               Epi, KMsurv, gnm, cmprsk, mstate, flexsurv, splines, epitools, 
               eha, shiny, ctqr, scales)
library(Rcpp)


orca <- read.table("orca.txt", stringsAsFactors = TRUE)
head(orca,n=10L)


ggplotly(
  orca %>%
    mutate(
      text = paste("Subject ID = ", id, "<br>", "Time = ", time, "<br>", "Event = ",  
                   event, "<br>", "Age = ", round(age, 2), "<br>", "Stage = ", stage)
    ) %>%
    ggplot(aes(x = id, y = time, text = text)) +
    geom_linerange(aes(ymin = 0, ymax = time)) +
    geom_point(aes(shape = event, color = event), stroke = 1, cex = 2) +
    scale_shape_manual(values = c(1, 3, 4)) +
    labs(y = "alivetime(年)", x = "patientId") + coord_flip() + theme_classic(),
  tooltip = "text"
)


ggplotly(
  orca %>%
    mutate(age_orig = age,
           age_end = age + time) %>%
    ggplot(aes(x = id, y = age_end)) +
    geom_linerange(aes(ymin = age_orig, ymax = age_end)) +
    geom_point(aes(shape = event, color = event), stroke = 1, cex = 2) +
    scale_shape_manual(values = c(1, 3, 4)) + guides(fill = FALSE) +
    labs(y = "age）", x = "patientId") + coord_flip() + theme_classic(),
  tooltip = "text"
)



orca <- mutate(orca, all = event != "Alive")
fit_km <- survfit(Surv(time, all) ~ 1, data = orca)
print(fit_km, print.rmean = TRUE)
dat_km <- fortify(fit_km)
head(dat_km,10)
ggsurvplot(fit_km,xlab = "Time (years)",surv.median.line = "hv", censor = F)


glist <- list(
  ggsurvplot(fit_km, fun = "event", main = "Cumulative proportion"),
  ggsurvplot(fit_km, fun = "cumhaz",  main = "Cumulative Hazard"),
  ggsurvplot(fit_km, fun = "cloglog", main = "Complementary log−log")
)
arrange_ggsurvplots(glist, print = TRUE, ncol = 3, nrow = 1)


fit_km2 <- survfit(Surv(time, all) ~ sex, data = orca)

ggsurvplot(fit_km2,xlab = "Time (years)",surv.median.line = "hv", censor = F)



#cox
m1 <- coxph(Surv(time, all) ~ sex + I((age-65)/10) + stage, data = orca)
summary(m1)

cox.zph.m1 <- cox.zph(m1)
cox.zph.m1

#residual
ggcoxdiagnostics(m1, type = "dfbeta", linear.predictions = FALSE)

#forest
ggforest(m1, data = orca)


#flexsurvreg
m2w <- flexsurvreg(Surv(time, all) ~ sex + I((age-65)/10) + stage,
                   data = orca, dist = "weibull")
m2w



#compare
fit_exp <- flexsurvreg(Surv(time, all) ~ 1, data = orca, dist = "exponential")
fit_w <- flexsurvreg(Surv(time, all) ~ 1, data = orca, dist = "weibull")
fit_gp <- flexsurvreg(Surv(time, all) ~ 1, data = orca, dist = "gompertz")
fit_sp <- flexsurvspline(Surv(time, all) ~ 1, data = orca, k = 1, scale = "odds")
grid.arrange(
  ggplot(data.frame(summary(fit_exp)), aes(x = time)) + 
    geom_line(aes(y = est, col = "Exponential")) +
    geom_line(data = data.frame(summary(fit_w)), aes(y = est, col = "Weibull")) +
    geom_line(data = data.frame(summary(fit_gp)), aes(y = est, col = "Gompertz")) +
    geom_line(data = data.frame(summary(fit_sp)), aes(y = est, col = "Flex splines")) +
    labs(x = "Time (years)", y = "Survival", col = "Distributions") + theme_classic(),
  ggplot(data.frame(summary(fit_exp, type = "hazard")), aes(x = time)) + 
    geom_line(aes(y = est, col = "Exponential")) +
    geom_line(data = data.frame(summary(fit_w, type = "hazard")), aes(y = est, col = "Weibull")) +
    geom_line(data = data.frame(summary(fit_gp, type = "hazard")), aes(y = est, col = "Gompertz")) +
    geom_line(data = data.frame(summary(fit_sp, type = "hazard")), aes(y = est, col = "Flex splines")) +
    labs(x = "Time (years)", y = "Hazard", col = "Distributions") + theme_classic(),
  ncol = 2
)

















