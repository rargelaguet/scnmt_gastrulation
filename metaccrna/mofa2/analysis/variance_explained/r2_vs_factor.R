library(MOFA2)
library(data.table)
library(purrr)
library(ggpubr)

#####################
## Define settings ##
#####################

io <- list()
io$model <- "/Users/ricard/data/gastrulation/mofa2/hdf5/test_1.hdf5"
io$outdir <- "/Users/ricard/data/gastrulation/mofa2/pdf/variance_explained"

################
## Load model ##
################

model <- load_model(file = io$model)

plot_factor_cor(model)

# subset factors
# model <- subset_factors(model, factors=1:10)

# rename views
tmp <- c(
  "met_Promoters" = "Promoter methylation",
  "met_E7.5 enhancers" = "Enhancer methylation",
  "acc_Promoters" = "Promoter accessibility",
  "acc_E7.5 enhancers" = "Enhancer accessibility",
  "RNA" = "RNA expression"
)
views(model) = stringr::str_replace_all(views(model), tmp[views(model)])

#######################################
## Plot variance explained vs factor ##
#######################################

r2 <- model@cache$variance_explained$r2_per_factor
r2 <- Reduce("+", r2)

r2.dt <- r2 %>%
  as.data.table %>% .[,factor:=as.factor(1:model@dimensions$K)] %>%
  melt(id.vars=c("factor"),variable.name="view", value.name = "r2") %>%
  .[,r2:=r2*100] %>%
  .[,cum_r2:=cumsum(r2), by="view"]

p <- ggline(r2.dt, x="factor", y="cum_r2", color="view") +
  scale_color_brewer(palette = "Dark2") +
  labs(x="Factor number", y="Cumulative variance explained (%)") +
  # theme(legend.title = element_blank())
  theme(legend.title = element_blank(), legend.position = "top")

pdf(paste0(io$outdir,"/r2_vs_factor.pdf"), width=6, height=4, useDingbats = F)
print(p)
dev.off()
