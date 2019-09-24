# -------------------------------------------------------------------
# - NAME:   circforest_fit_models.R
# - AUTHOR: Moritz N. Lang, Lisa Schlosser
# - DATE:   2019-09-13
# -------------------------------------------------------------------
# - PURPOSE:
# -------------------------------------------------------------------
# - L@ST MODIFIED: 2019-09-24 on thinkmoritz
# -------------------------------------------------------------------

# -------------------------------------------------------------------
# Preliminaries
# -------------------------------------------------------------------
## Set time to utc (just in case)
Sys.setenv('TZ'='UTC')

## Load packages
library("optparse")
library("zoo")
library("circtree")
library("disttree")
library("verification") # Careful, version 1.35 needed!!
library("ggplot2")
library("reshape2")
library("ggpubr")

## Create folder for outputs
if (! dir.exists("results")) dir.create("results")


# -------------------------------------------------------------------
# NAMELIST PARAMETERS
# -------------------------------------------------------------------
option_list <- list(
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
    help = "Print extra output [default]"),
  make_option(c("-q", "--quietly"), action = "store_false",
    dest = "verbose", help = "Print little output"),
  make_option("--run_name", type = "character", default = "v7",
    help = "Run name or version of script used for output name [default \"%default\"]"),
  make_option("--plot", action = "store_true", default = FALSE,
    help = "Plot validation [default]"),
  make_option(c("--seed"), type = "integer", default = 123,
    help="Set seed for bootstrapping [default %default]")
  )

opt <- parse_args(OptionParser(option_list = option_list))


## Set theme
theme_set(theme_bw(base_size = 14.5) +
   theme(panel.grid.major = element_line(linetype = "dotted", colour = "grey80"),
         panel.grid.minor = element_blank(),
         plot.title = element_text(hjust = 0.5),
         strip.text.x = element_text(margin = margin(.15, 0, .15, 0, "cm")),
         strip.text.y = element_text(margin = margin(0, 0.15, 0, 0.15, "cm")),
         axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))))

## -------------------------------------------------------------------
## Loading data and plots per station and lag
## -------------------------------------------------------------------
## Loop over stations and lag
crps.mall <- crps.agg_skill.mall <- data.frame()

for (i_station in c("ibk", "vie")){
  for (i_lag in c("6", "18")){

    load(file = sprintf("results/circforest_validation_%s_lag%s_%s.rda", i_station, i_lag, opt$run_name))

    ## Plot raw crps
    crps.m <- melt(crps[, c("climatology", "persistence", "tree", "forest")], variable.name = "model")
    crps.m$model <- plyr::revalue(crps.m$model,
      c("climatology" = "Climatology",
      "persistence" = "Persistence",
      "tree" = "Tree",
      "forest" = "Forest"))
    crps.m$station <- factor(i_station, levels = c("ibk", "vie"), labels = c("Innsbruck", "Vienna"))
    crps.m$lag <- factor(i_lag, levels = c("6", "18"), labels = c("1-hour forecast", "3-hour forecast"))

    crps.agg_skill.m <- melt(crps.agg_skill, variable.name = "model")
    crps.agg_skill.m$model <- plyr::revalue(crps.agg_skill.m$model,
      c("climatology" = "Climatology",
      "persistence" = "Persistence",
      "tree" = "Tree",
      "forest" = "Forest"))
    crps.agg_skill.m$station <- factor(i_station, levels = c("ibk", "vie"), labels = c("Innsbruck", "Vienna"))
    crps.agg_skill.m$lag <- factor(i_lag, levels = c("6", "18"), labels = c("1-hour forecast", "3-hour forecast"))

    crps.mall <- rbind(crps.mall, crps.m)
    crps.agg_skill.mall <- rbind(crps.agg_skill.mall, crps.agg_skill.m)

    ## Plot raw crps    
    p1 <- ggplot(crps.m, aes(x = model, y = value)) +
          geom_hline(yintercept = 0, linetype ="solid", colour = "gray80") +
          stat_boxplot(geom = "errorbar", width = 0.2) +
          geom_boxplot(fill = "gray60")
    p1 <- p1 + labs(x = "", y = "CRPS [rad]") 

    dev.new(width=8, height=4.5)
    print(ggarrange(p1, legend = "none"))
    pdf_file <- sprintf("results/_plot_circforest_validation_crpsraw_%s_lag%s_%s.pdf",
      i_station, i_lag, opt$run_name)
    ggsave(pdf_file)
    
    ## Plot crps skill scores   
    p2 <- ggplot(crps.agg_skill.m, aes(x = model, y = value)) +
          geom_hline(yintercept = 0, linetype ="solid", colour = "gray80") +
          stat_boxplot(geom = "errorbar", width = 0.2) +
          geom_boxplot(fill = "gray60")
    p2 <- p2 + labs(x = "", y = "CRPS skill score [%]") 
    
    dev.new(width=8, height=4.5)
    print(ggarrange(p2, legend = "none"))
    pdf_file <- sprintf("results/_plot_circforest_validation_crpsskill_agg_%s_lag%s_%s.pdf",
      i_station, i_lag, opt$run_name)
    ggsave(pdf_file)
    
    ## Plot single tree
    m_ct.plot <- readRDS(file = sprintf("results/circforest_model_tree_%s_lag%s_%s_4plotting.rds",    
      i_station, i_lag, opt$run_name))
    
    pdf(file = sprintf("results/_plot_circforest_exampletree_%s_lag%s_%s.pdf", 
      i_station, i_lag, opt$run_name), width = 18, height = 10)
    par(mar = c(3.1, 4.1, 2.1, 2.1))
    plot(m_ct.plot, ep_args = list(justmin = 10), tp_args = list(type = "response", plot_type = "geographics"), 
      ip_args = list(pval = FALSE))
    dev.off()
  }
}

## -------------------------------------------------------------------
## Summary plots
## -------------------------------------------------------------------
## Make summary plots with four panels (cprs raw)
p3 <- ggplot(crps.mall, aes(x = model, y = value)) +
      geom_hline(yintercept = 0, linetype ="solid", colour = "gray80") +
      stat_boxplot(geom = "errorbar", width = 0.2) +
      geom_boxplot(fill = "gray60") + 
      facet_grid(lag ~ station, scales = "free")
p3 <- p3 + labs(x = "", y = "CRPS [rad]") 
p3 <- p3 + annotate("text", -Inf, Inf, label = paste0("(", letters[1:4], ")"), hjust = -0.2, vjust = 1.3)

dev.new(width=10, height=6.5)
print(ggarrange(p3, legend = "none"))
pdf_file <- sprintf("results/_plot_circforest_validation_crpsraw_comparison_%s.pdf", opt$run_name)
ggsave(pdf_file)

## Make summary plots with four panels (cprs skill score)
p4 <- ggplot(crps.agg_skill.mall, aes(x = model, y = value)) +
      geom_hline(yintercept = 0, linetype ="solid", colour = "gray80") +
      stat_boxplot(geom = "errorbar", width = 0.2) +
      geom_boxplot(fill = "gray60") + 
      facet_grid(lag ~ station, scales = "free")
p4 <- p4 + labs(x = "", y = "CRPS skill score [%]") 
p4 <- p4 + annotate("text", -Inf, Inf, label = paste0("(", letters[1:4], ")"), hjust = -0.2, vjust = 1.3)

dev.new(width=10, height=6.5)
print(ggarrange(p4, legend = "none"))
pdf_file <- sprintf("results/_plot_circforest_validation_crpsskill_agg_comparison_%s.pdf", opt$run_name)
ggsave(pdf_file)
