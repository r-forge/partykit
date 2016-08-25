
library("partyNG")
library("sandwich")
library("libcoin")
data("PimaIndiansDiabetes", package = "mlbench")


## partykit: recursive partitioning of a logistic regression model
pid_tree <- partykit::glmtree(diabetes ~ glucose | pregnant +
       pressure + triceps + insulin + mass + pedigree + age,
       data = PimaIndiansDiabetes, family = binomial)

pid_tree

nd <- factor(predict(pid_tree, type = "node"))
m <- glm(diabetes ~ glucose:nd, data = PimaIndiansDiabetes, family = binomial)
logLik(m)
AIC(m)

pid_tree1 <- partyNG::glmtree(diabetes ~ glucose | pregnant +
       pressure + triceps + insulin + mass + pedigree + age,
       data = PimaIndiansDiabetes, family = binomial)

pid_tree1

info_node(node_party(pid_tree1))

nd1 <- factor(predict(pid_tree1, type = "node"))
m1 <- glm(diabetes ~ glucose:nd1, data = PimaIndiansDiabetes, family = binomial)
logLik(m1)
AIC(m1)


table(nd, nd1)

pid_tree2 <- partyNG::glmtree(diabetes ~ glucose | pregnant +
       pressure + triceps + insulin + mass + pedigree + age,
       data = PimaIndiansDiabetes, family = binomial, 
       control = ctree_control(nmax = 20))

pid_tree2

info_node(node_party(pid_tree2))

nd2 <- factor(predict(pid_tree2, type = "node"))
m2 <- glm(diabetes ~ glucose:nd2, data = PimaIndiansDiabetes, family = binomial)
logLik(m2)
AIC(m2)


table(nd, nd2)
table(nd1, nd2)
