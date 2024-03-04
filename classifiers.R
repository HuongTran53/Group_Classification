##-- Filename: classifiers.R
## Classify Crabs data set 
Crabs <- read.table("http://www.stat.ufl.edu/~aa/cat/data/Crabs.dat", header=TRUE)
head(Crabs)
set.seeed(1)
split <- rbinom(nrow(Crabs), 1, .3)
test <- Crabs[split == 1, ]
train <- Crabs[split == 0, ]
library(ggplot2)
train$y <- as.factor(train$y)
ggplot(train, aes(width, weight, color = y)) + geom_point()

#Logistic 
mod.glm <- glm(y ~ width + weight, train, family = binomial)
fit.glm <- predict(mod.glm, train, type = "response")
fit.glm <- ifelse(fit.glm > .5, 1, 0)
test.glm <- ifelse(predict(mod.glm, test, type = "response") > .5, 1, 0)
table(train$y, fit.glm)
table(test$y,  test.glm)
library(MASS)
# LDA
mod.lda <- MASS::lda(as.factor(y) ~ width + weight, train, prior = c(.5, .5))
fit.lda <- predict(mod.lda, train)
test.lda <- predict(mod.lda, test)
table(train$y,  fit.lda$class)
table(~test$y,   test.lda$class)
# qda
mod.qda <- MASS::qda(as.factor(y) ~ width + weight, train, prior = c(.5, .5))
fit.qda <- predict(mod.qda, train)
test.qda <- predict(mod.qda, test)
table(train$y,  fit.qda$class)
table(~test$y,   test.qda$class)
# SVM
library("e1071")
mod.svm <- e1071::svm(as.factor(y) ~ width + weight, data = train, kernel = "linear")
fit.svm <- predict(mod.svm, train)
table(train$y, fit.svm)

mod.svm <- e1071::svm(as.factor(y) ~ width + weight, data = train, kernel = "radial")
fit.svm <- predict(mod.svm, train)
table(train$y, fit.svm)

#Naive Bayes
library(naivebayes)
mod.nb <- naivebayes::naive_bayes(y ~ width + weight, data = train, usekernel = T) 
fit.nb <- predict(mod.nb, train)
table(train$y, fit.nb)

# Decision Tree
mod.tree <- rpart::rpart(y ~ width + weight, data = train)
rattle::fancyRpartPlot(mod.tree)
fit.tree <- predict(mod.tree, data = train, type = "class")
table(train$y, fit.tree)

# ANN
install.packages(c('neuralnet','keras','tensorflow'),dependencies = T)
library(tidyverse)
library(neuralnet)
mod.ann <- neuralnet::neuralnet(y ~ width + weight, 
                                data = train, 
                                hidden = c(3), 
                                linear.output = FALSE)
plot(mod.ann, rep = "best")
fit.ann <- predict(mod.ann, train, type = "class") 
fit.ann <- apply(fit.ann, 1, which.max)
fit.ann <- fit.ann - 1
table(train$y, fit.ann)
  
