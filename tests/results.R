library(normalmultinomial)
library(dplyr)
library(ggplot2)

source('scenarios.R')

get_parameters = function(x) gsub(paste0(x, ' = '), '', grep(paste0(x, ' = '), readLines('Makefile'), value = TRUE))

L_N = as.numeric(get_parameters('L_N'))
L_n = scan(text = get_parameters('L_n'))
L_s = scan(text = get_parameters('L_s'))
L_seed = scan(text = get_parameters('L_seed'))
L_method = scan(text = get_parameters('L_method'), what = character())

names(L_N) = L_N
names(L_n) = L_n
names(L_s) = L_s
names(L_seed) = L_seed
names(L_method) = L_method

compare = function(H1, H2) mean(apply(H1-H2, 1, function(x) sqrt(sum(x^2))))
evaluate = function(N, n, s, seed, method){
  load(sprintf('datasets/replacement-N_%05d-n_%05d-s_%05d-seed_%05d-method_%s.RData',
               N, n, s, seed, method))
  
  p = params[[s]]$p
  P.gs = matrix(p, nrow = N, ncol = length(p), byrow = TRUE)
  H.gs = ilr_coordinates(P.gs)
  
  P = fit$expected
  H = ilr_coordinates(P)
  
  data_frame(
    'paired.dist' = compare(H, H.gs))
}


df = lapply(L_method, function(method){
  lapply(L_seed, function(seed, method){
    lapply(L_s, function(s, seed, method){
      lapply(L_n, function(n, s, seed, method){
        lapply(L_N, function(N, n, s, seed, method){
          evaluate(N, n, s, seed, method)
        }, n, s, seed, method) %>% bind_rows(.id = 'N')
      }, s, seed, method) %>% bind_rows(.id = 'n')
    }, seed, method) %>% bind_rows(.id = 's')
  }, method) %>% bind_rows(.id = 'seed')
}) %>% bind_rows(.id = 'method')

df = df %>% 
  group_by(method, s, n, N) %>%
  summarise(paired.dist = mean(paired.dist))

df$n = as.factor(as.numeric(df$n))
df$s = sprintf("Scenario %s", df$s)
g <- ggplot() + 
  geom_bar(data=df, aes(n, paired.dist, fill = method), col='black', stat='identity', position = 'dodge') +
  facet_wrap(~s, scales = 'free') +
  theme_classic() + 
  xlab('Counting total') + ylab('Mean differences') + 
  theme(legend.position = 'top') +
  scale_fill_manual(values=c("grey", "white"), 
                    name="Methods:",
                    breaks=c("dm", "nm"),
                    labels=c("Dirichlet-Multinomial", "Log-ratio-Normal-Multinomial"))


ggsave('figures/multinomial.pdf', g, width = 7, height = 5.25)
