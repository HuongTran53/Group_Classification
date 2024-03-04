library(dplyr)
library(ggplot2)
library(GGally)
mytheme <-  ggplot2::theme_bw() + 
  ggplot2::theme(
    aspect.ratio = 1/1, 
    panel.grid = element_blank(), 
    axis.line = element_line(colour = "black"
    ), 
    axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"), 
    legend.background = element_rect( 
      size=0.5, linetype="solid"), 
    legend.text = element_text(size=12))
# Data processing
df <- 
  read.csv(file = "/Users/huongtran/OU /Dissertation/Datasets/ETF prices.csv")
tmp <- df %>% mutate_at(c("price_date"), function(x) as.Date(x)) %>%  
  filter(price_date >= as.Date("2016-01-01")) %>%
  group_by(fund_symbol) %>% 
  mutate(n_ = n_distinct(price_date)) %>% 
  ungroup() 
chosen_symbol <- 
  unique(tmp[tmp$n_ == max(tmp$n_) , "fund_symbol"])$fund_symbol
mysymbol <- c("IWF", "LRGF")
df_stock <- tmp %>%
  filter(fund_symbol %in% mysymbol) %>% 
  mutate_at(c("price_date"), function(x) as.Date(x))

fund1 <- df_stock %>% 
  filter(fund_symbol == mysymbol[1]) %>% 
  arrange(price_date)
fund1 <- (diff(fund1$close)/fund1$close[-nrow(fund1)])

fund2 <- df_stock %>% 
  filter(fund_symbol == mysymbol[2]) %>% 
  arrange(price_date)
fund2 <- (diff(fund2$close)/fund2$close[-nrow(fund2)])

mydf <- data.frame(fund1 = fund1, fund2 = fund2)

p <- ggplot2::ggplot(mydf, aes(fund1, fund2)) + 
  ggplot2::geom_point(shape = 1, color = "#4292c6", size = 2) + 
  mytheme + xlab(mysymbol[1]) + ylab(mysymbol[2])
ggExtra::ggMarginal(p, type = "histogram", fill = "#9ECAE1")
###################################
U <- copula::pobs(mydf);
p <- ncol(U); nx <- nrow(U)
ggplot2::ggplot(data.frame(U), aes(fund1, fund2)) + 
  ggplot2::geom_point(shape = 1, color = "#4292c6", size = 2) + 
  mytheme + xlab(mysymbol[1]) + ylab(mysymbol[2])
seed <- 1234
par(mfrow = c(2, 3))
set.seed(seed)
frank.param <- coef(
  copula::fitCopula(
    copula = copula::frankCopula(dim = p), data = U, method = "mpl"
  )
)
frank.cop <- copula::frankCopula(param = frank.param, dim = p)
frank.Y <- copula::rCopula(nx, copula = frank.cop)
frank.g <- ggplot2::ggplot(data.frame(frank.Y), aes(X1, X2)) +
  ggplot2::geom_point(shape = 1, color = "#4292c6", size = 1) + 
  mytheme + ggtitle("Frank") + 
  theme(plot.title = element_text(face = "bold"), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank())
set.seed(seed)
gumbel.param <-coef(
  copula::fitCopula(
    copula = copula::gumbelCopula(dim = p), data = U, method = "mpl"
  )
)
gumbel.cop <- copula::gumbelCopula(param = gumbel.param, dim = p)
gumbel.Y <- copula::rCopula(nx, copula = gumbel.cop)
gumbel.g <- ggplot2::ggplot(data.frame(gumbel.Y), aes(X1, X2)) +
  ggplot2::geom_point(shape = 1, color = "#4292c6", size = 1) + 
  mytheme + xlab(bquote(U[1])) + ylab(bquote(U[2])) + ggtitle("Gumbel") + 
  theme(plot.title = element_text(face = "bold"), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank())
# 
set.seed(seed)
clayton.param <- coef(
  copula::fitCopula(
    copula = copula::claytonCopula(dim = p), data = U, method = "mpl"
  )
)
clayton.cop <- copula::claytonCopula(param = clayton.param, dim = p)
clayton.Y <- copula::rCopula(nx, copula = clayton.cop)
clayton.g <- ggplot2::ggplot(data.frame(clayton.Y), aes(X1, X2)) +
  ggplot2::geom_point(shape = 1, color = "#4292c6", size = 1) + 
  mytheme + xlab(bquote(U[1])) + ylab(bquote(U[2])) + ggtitle("Clayton") + 
  theme(plot.title = element_text(face = "bold"), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank())
# 
set.seed(seed)
joe.param <- coef(
  copula::fitCopula(
    copula = copula::joeCopula(dim= p), data = U, method = "mpl"
  )
)
joe.cop <- copula::joeCopula(param = joe.param, dim = p)
joe.Y <-  copula::rCopula(nx, copula = joe.cop)
joe.g <- ggplot2::ggplot(data.frame(joe.Y), aes(X1, X2)) +
  ggplot2::geom_point(shape = 1, color = "#4292c6", size = 1) + 
  mytheme + xlab(bquote(U[1])) + ylab(bquote(U[2])) + ggtitle("Joe") + 
  theme(plot.title = element_text(face = "bold"), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank())
#
set.seed(seed)
normal.param <- coef(
  copula::fitCopula(
    copula = copula::normalCopula(dim= p), data = U, method = "mpl"
  )
)
normal.cop <- copula::normalCopula(param = normal.param, dim = p)
normal.Y <-  copula::rCopula(nx, copula = normal.cop)
normal.g <- ggplot2::ggplot(data.frame(normal.Y), aes(X1, X2)) +
  ggplot2::geom_point(shape = 1, color = "#4292c6", size = 1) + 
  mytheme + xlab(bquote(U[1])) + ylab(bquote(U[2])) + ggtitle("Normal") + 
  theme(plot.title = element_text(face = "bold"), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank())
# 
set.seed(seed)
t.param <- coef(
  copula::fitCopula(
    copula = copula::tCopula(dim= p), data = U, method = "mpl"
  )
)
t.cop <- copula::tCopula(
  param = c(t.param["rho.1"]), df = t.param["df"], dim = p
)
t.Y <-  copula::rCopula(nx, copula = t.cop)
t.g <- ggplot2::ggplot(data.frame(t.Y), aes(X1, X2)) +
  ggplot2::geom_point(shape = 1, color = "#4292c6", size = 1) + 
  mytheme + xlab(bquote(U[1])) + ylab(bquote(U[2])) + ggtitle("t") + 
  theme(plot.title = element_text(face = "bold"), 
        axis.title.x=element_blank(), 
        axis.title.y=element_blank())
# 
gridExtra::grid.arrange(grobs = list(frank.g, gumbel.g, clayton.g, 
                                     joe.g, normal.g, t.g), 
                        padding = unit(0, "line"),
                        layout_matrix = matrix(1:6, nrow = 3, byrow = T))

lstX <- list(frank.Y, gumbel.Y, clayton.Y, joe.Y, normal.Y, t.Y)
lstCop <- list(frank.cop, gumbel.cop, clayton.cop, 
               joe.cop, normal.cop, t.cop)

re.eigen <- mul.discr.eigen.reg(lstX, U, w = rep(1, 3))
re.NP <- mul.NP(lstCop, U)

re_df <- data.frame(Eigen_reg = re.eigen$value, 
                    NP = re.NP$value)
rownames(re_df) <- c("Frank", "Gumble", "Clayton", "Joe","Normal", "t")
re_df
c(
  Eigen_reg = re.eigen$pop, 
  NP = re.NP$pop
)

############# 4 funds #####################
mysymbol <- c("IWF", "LRGF", "ONEO", "OUSA")
df_stock <- tmp %>%
  filter(fund_symbol %in% mysymbol) %>% 
  mutate_at(c("price_date"), function(x) as.Date(x)) #%>% 
# select(c("price_date", "fund_symbol", "close"))

funds <- sapply(mysymbol, function(x) {
  y <- df_stock %>% filter(fund_symbol == x) %>% arrange(price_date)
  (diff(y$close)/y$close[-nrow(y)])
}, simplify = "array")

mydf <- as.data.frame(funds)
GGally::ggpairs(
  mydf,
  diag = list(continuous = wrap(ggally_barDiag, fill = "#9ECAE1")), 
  lower = list(
    continuous = wrap(ggally_points, shape = 1, color = "#4292c6", size = 1)
  ), 
  upper = list(
    continuous = wrap(ggally_cor, method = "spearman", 
                      fontface ="bold", colour = "black")
  )
) + 
  theme_bw() + theme(panel.grid = element_blank())
###################################
U <- copula::pobs(mydf)
p <- ncol(U); nx <- 16*nrow(U)
GGally::ggpairs(
  data.frame(U),
  upper = list(
    continuous = wrap(
      ggally_cor, method = "spearman", fontface ="bold", colour = "black")
  ), 
  diag = list(continuous = "blankDiag"), 
  lower = list(
    continuous = wrap(
      ggally_points, shape = 1, color = "#4292c6", size = 1)
  )
) + 
  theme_bw() + theme(panel.grid = element_blank())
###################################
par(mfrow = c(2, 3))
set.seed(seed)
frank.param <- coef(
  copula::fitCopula(
    copula = copula::frankCopula(dim = p), data = U, method = "mpl"
  )
)
frank.Y <- copula::rCopula(
  nx, copula = copula::frankCopula(param = frank.param, dim = p)
)
frank.cop <- copula::frankCopula(param = frank.param, dim = p)
set.seed(seed)
gumbel.param <-coef(
  copula::fitCopula(
    copula = copula::gumbelCopula(dim = p), data = U, method = "mpl"
  )
)
gumbel.Y <- copula::rCopula(
  nx, copula = copula::gumbelCopula(param = gumbel.param, dim = p)
)
gumbel.cop <- copula::gumbelCopula(param = gumbel.param, dim = p)

set.seed(seed)
clayton.param <- coef(
  copula::fitCopula(
    copula = copula::claytonCopula(dim = p), data = U, method = "mpl"
  )
)
clayton.Y <- copula::rCopula(
  nx, copula = copula::claytonCopula(param = clayton.param, dim = p)
)
clayton.cop <- copula::claytonCopula(param = clayton.param, dim = p)

set.seed(seed)
joe.param <- coef(
  copula::fitCopula(
    copula = copula::joeCopula(dim= p), data = U, method = "mpl"
  )
)
joe.Y <-  copula::rCopula(
  nx, copula = copula::joeCopula(param = joe.param, dim = p)
)
joe.cop <- copula::joeCopula(param = joe.param, dim = p)

set.seed(seed)
normal.param <- coef(
  copula::fitCopula(
    copula = copula::normalCopula(dim= p), data = U, method = "mpl"
  )
)
normal.Y <-  copula::rCopula(
  nx, copula = copula::normalCopula(param = normal.param, dim = p)
)
normal.cop <- copula::normalCopula(param = normal.param, dim = p)

set.seed(seed)
t.param <- coef(
  copula::fitCopula(
    copula = copula::tCopula(dim= p), data = U, method = "mpl"
  )
)
# t.param[ "df"] <- round(t.param[ "df"], 1)
t.cop <- copula::tCopula(
  param = c(t.param["rho.1"]), df = t.param["df"], dim = p
)
t.Y <-  copula::rCopula(
  nx, copula = copula::tCopula(
    param = c(t.param["rho.1"]), df = t.param["df"], dim = p
  )
)
# lstX <- list(frank.Y, gumbel.Y, clayton.Y, joe.Y, normal.Y, t.Y)
# lstCop <- list(frank.cop, gumbel.cop, 
#                clayton.cop, joe.cop, normal.cop, t.cop)
lstX <- list(frank.Y, gumbel.Y, clayton.Y, joe.Y, normal.Y, t.Y)
lstCop <- list(frank.cop, gumbel.cop, clayton.cop, 
               joe.cop, normal.cop, t.cop)

re.eigen <- mul.discr.eigen.reg(lstX, U, w = rep(1, 3))
re.NP <- mul.NP(lstCop, U)

re_df <- data.frame(Eigen_reg = re.eigen$value, 
                    NP = re.NP$value)
rownames(re_df) <- c("Frank","Gumble", "Clayton", "Joe","Normal", "t")
re_df
c(
  Eigen_reg = re.eigen$pop, 
  NP = re.NP$pop
)


