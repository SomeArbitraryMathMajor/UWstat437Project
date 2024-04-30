library(survival)
library(randomForestSRC)
library(survAUC)
library(survminer)
library(dplyr)
library(finalfit)
library(patchwork)

data(cancer, package = 'survival')
df <- na.omit(mgus2)

## splitting data for training
X <- df[, 2:7] # (2:7)[-3]
died <- df$death == 1
times <- df$futime
split <- iai::split_data('survival', X, died, times, seed=437)
train_X <- split$train$X
train_died <- split$train$deaths
train_times <- split$train$times
test_X <- split$test$X
test_died <- split$test$deaths
test_times <- split$test$times

## corresponding split for coxph and rsf
train.df <- cbind(train_times, train_died, train_X)
names(train.df) <- c('time', 'status', names(train.df)[-c(1,2)])
test.df <- cbind(test_times, test_died, test_X)
names(test.df) <- c('time', 'status', names(test.df)[-c(1,2)])

###########################
# EDA
###########################

## check for NA
head(mgus2)
sum(is.na(mgus2))

## life table
survival_object <- mgus2 %$% 
  Surv(futime, death)
fit <- survfit(survival_object ~ 1, data = mgus2)
summary(fit, times = seq(0, 50, by=10))

## Kaplan-Meier for sex
dependent <- 'Surv(futime, death)'
explanatory <- 'sex'
mgus2 %>%
  surv_plot(
    dependent, explanatory, pval=F, size=0.7,
    censor.size=2, risk.table.fontsize=3,
    risk.table=T)

## training and testing sets for Cox and RSF
head(train.df)
head(test.df)

###########################
# models
###########################

## fit Cox models and visualize
fit.cox <- coxph(
  Surv(time, status) ~ .,
  data=train.df); summary(fit.cox)
fit.cox.test <- coxph(
  Surv(time, status) ~ .,
  data=test.df); summary(fit.cox.test)
p1 <- ggadjustedcurves(fit.cox, data=train.df, variable='sex')
p2 <- ggadjustedcurves(fit.cox.test, data=test.df, variable = 'sex')
# survival curve for sex
p1 + p2

## fit and visualize RSF models
fit.tree.rfsrc <- rfsrc(
  Surv(time, status) ~ .,
  data=train.df, ntree=1, nodedepth=3); fit.tree.rfsrc
plot(get.tree(fit.tree.rfsrc, 1))
fit.tree.rfsrc.test <- rfsrc(
  Surv(time, status) ~ .,
  data=test.df, ntree=1, nodedepth=3); fit.tree.rfsrc.test
plot(get.tree(fit.tree.rfsrc.test, 1))

## optimal survival tree
grid <- iai::grid_search(
  iai::optimal_tree_survival_learner(
    random_seed = 1,
    missingdatamode = 'separate_class',
    minbucket = 1,
  ),
  max_depth = 1:3,
)
fit.ost <- iai::fit(grid, train_X, train_died, train_times,
                    validation_criterion = 'harrell_c_statistic')
iai::get_learner(grid)
iai::get_grid_result_summary(grid)

###########################
# comparison
###########################

## concordance value
cox.c.stat <- c(concordance(fit.cox)$concordance,
                concordance(fit.cox.test)$concordance)
tree.c.fn <- function(x, died, times)
  iai::score(grid, x, died, times, criterion='harrell_c_statistic')
ost.c.stat <- c(
  tree.c.fn(train_X, train_died, train_times),
  tree.c.fn(test_X, test_died, test_times))
tree.c.stat <- c(1-fit.tree.rfsrc$err.rate,
                 1-fit.tree.rfsrc.test$err.rate)
# table
c.stat <- rbind(ost.c.stat, tree.c.stat, cox.c.stat)
colnames(c.stat) <- c('training', 'testing')

## OST AUC
ost.auc <- c()
treeaucmap <- function(t) {
  tryCatch(
    expr = {
      iai::score(grid, test_X, test_died, test_times,
                 criterion = 'auc', evaluation_time=t)
    },
    error = function(e){
      message(paste('error at:',t))
    }
  )
}
for (i in test_times) {
  ost.auc <- c(ost.auc, treeaucmap(i))
}
## tree AUC
pred.struct <- predict(fit.tree.rfsrc)$yvar
lp <- pred.struct[,2]
pred.struct.test <- predict(fit.tree.rfsrc, newdata=test.df)$yvar
lpnew <- pred.struct.test[,2]
Surv.rsp <- Surv(train.df$time, train.df$status)
Surv.rsp.new <- Surv(test.df$time, test.df$status)
tree.auc <- AUC.cd(
  Surv.rsp, Surv.rsp.new, lp, lpnew, test_times)$auc
## cox AUC
lp <- predict(fit.cox)
lpnew <- predict(fit.cox, newdata=test.df)
cox.auc <- AUC.cd(Surv.rsp, Surv.rsp.new, lp, lpnew, test_times)$auc
# table
auc.stat <- rbind(ost.auc, tree.auc, cox.auc)
colnames(auc.stat) <- test_times

## comparison results
round(c.stat, 4)
round(rowMeans(auc.stat), 4)
