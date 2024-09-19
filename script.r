###Code accompanying Bianchi et al 2024###

#Libraries
library(rethinking)
library(ggplot2)
library(tidyverse)
library(gridExtra)
library(rstan)
library(grid)
library(gt)
library(gtsummary)
library(tidybayes)

#Stan optimization
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

#read file
df <- read.csv("df.csv", stringsAsFactors = FALSE)

#use GTSummary to make a demo table (Table 1 in manuscript)
t1 <- tbl_summary(df[c(9, 1, 4, 10, 6, 8, 7,11)], by = Group,
  missing = "no") %>%
  bold_labels() %>%
  as_gt() %>%
  gt::gtsave("t1.html")

##Modelling##

#create datalist
dat <- list(N = nrow(df),
            injury = as.integer(df1$`Future MSI`),
            K = length(unique(df$block_id)),
            block_id = as.integer(df$block_id),
            G = length(unique(df$Group)),
            group = as.integer(as.factor(df$Group)),
            X = length(unique(df$Sex)), #number of sexes
            sex = as.integer(as.factor(df$Sex))
)

str(dat)

#prior_dat list
dat_prior <- list(N = nrow(df),
            injury = as.integer(df$`Future MSI`)*0,
            K = length(unique(df$block_id)),
            block_id = as.integer(df$block_id),
            G = length(unique(df$Group)),
            group = as.integer(as.factor(df$Group)),
            X = length(unique(df$Sex)), #number of sexes
            sex = as.integer(as.factor(df$Sex))
)

str(dat_prior)

#compile models
logit_mod <- stan_model(file = "logit_mod.stan")
logit_mod_prior <- stan_model(file = "logit_mod_prior.stan")

#sample model
m <- sampling(logit_mod, data = dat, iter = 4000, chains = 4, cores = 4)

#prior model
m_prior <- sampling(logit_mod_prior, data = dat_prior, iter = 2000, chains = 4, cores = 4)

##posterior predictive analyses##

#extract posterior samples
m_post <- data.frame(extract.samples(m))
colnames(m_post)
precis(m_post, depth = 2, prob = 0.9)

#extract posterior prior samples
m_post_prior <- data.frame(extract.samples(m_prior))

#create a posterior predictive distribution at the group level

#extract simulated draws from the posterior
sim <- m_post[grepl("injury_sim", colnames(m_post))]

#calculate the probability for each column - this represents each data sample
sim2 <- as.data.frame(apply(sim, 2, mean))

#bind back to the original data 
sim_vec <- as.vector(sim2$`apply(sim, 2, mean)`)
df$sim <- sim_vec
df$sim <- df$sim * 100

#create a dataframe of raw probabilities to add to the PPD plot
raw_probs_df1 <- data.frame("Group" = c("SRC", "MSI", "Healthy Control"),
"prob" = c(mean(df$`Future MSI`[df$Group == "SRC" & df$Sex == "Male"]),
         mean(df$`Future MSI`[df$Group == "MSI" & df$Sex == "Male"]),
         mean(df$`Future MSI`[df$Group == "Healthy Control" & df$Sex == "Male"])))

raw_probs_df2 <- data.frame("Group" = c("SRC", "MSI", "Healthy Control"),
"prob" = c(mean(df$`Future MSI`[df$Group == "SRC" & df$Sex == "Female"]),
         mean(df$`Future MSI`[df$Group == "MSI" & df$Sex == "Female"]),
         mean(df$`Future MSI`[df$Group == "Healthy Control" & df$Sex == "Female"])))         

raw_probs_df <- rbind(raw_probs_df1, raw_probs_df2)
raw_probs_df$Sex <- rep(c("Male", "Female"), each = 3)
raw_probs_df$prob <- raw_probs_df$prob * 100

#plot posterior predictive distribution against real data
theme_set(theme_tidybayes() + panel_border())
sim_plot <- ggplot(df, aes(x = sim, y = Group, fill = Group)) +
   facet_grid( ~ Sex) +
   scale_fill_manual(values = c("green","blue","red")) +
  theme(strip.text = element_text(face = "bold",size = 12),
        axis.title = element_text(size = 10, face = "bold"),
        legend.position = "none",
        axis.text.y = element_text(face = "bold")) +
        ylab("") +
        xlab("Probability of Future MSI (%)") +
        stat_gradientinterval(.width = c(0.70,0.9)) +
        geom_point(data = raw_probs_df,
        aes(x = prob, y = Group, group = Sex), 
        color = "black",fill = "white", size = 4, shape = 21)
sim_plot

ggsave(sim_plot, file = "sim_fig.jpg", dpi = 600)

##Evaluate prior simulation##
colnames(m_post_prior)
#female
f_prior <- m_post_prior[, c(1200, 1201, 1202)]
colnames(f_prior)

#male col
m_prior <- m_post_prior[, c(1203, 1204, 1205)]
colnames(m_prior)

#change colnames to match
colnames(m_prior) <- c("Healthy", "MSI", "SRC")
colnames(f_prior) <- colnames(m_prior)
#add sex columns
f_prior$Sex <- "Female"
m_prior$Sex <- "Male"

#combine dataframes
posterior_prior_df <- rbind(f_prior, m_prior)

#gather for figure
posterior_prior_df2 <- gather(posterior_prior_df, key = "Group", 
                              value = "Prior Probability (%)", 
                              Healthy:SRC)

#convert to prob
posterior_prior_df2$`Prior Probability (%)` <-
posterior_prior_df2$`Prior Probability (%)` * 100

#plot
theme_set(theme_tidybayes() + panel_border())
plot_prior <- ggplot(posterior_prior_df2, aes(x = `Prior Probability (%)`, 
y = Group, fill = Group)) +
   facet_grid(~ Sex) +
   scale_fill_manual(values = c("green","blue","red")) +
  theme(strip.text = element_text(face = "bold",size = 12),
        axis.title = element_text(size = 10, face = "bold"),
        legend.position = "none",
        axis.text.y = element_text(face = "bold")) +
        ylab("") +
        xlab("Probability of Future MSI (%)") +
        stat_halfeye(.width = c(0.70,0.9)) 
plot_prior

ggsave(plot_prior, file = "prior.jpg", dpi = 600)

##Evaluate posterior samples##
colnames(m_post)

#female
f <- m_post[, c(1200, 1201, 1202)]
colnames(f) <- c("Healthy", "MSI", "SRC")

#male col
ma <- m_post[, c(1203, 1204, 1205)]
colnames(ma) <- colnames(f)

#add sex columns
f$Sex <- "Female"
ma$Sex <- "Male"

#combine dataframes
posterior_df <- rbind(f, ma)

#gather for figure
posterior_df2 <- gather(posterior_df, key = "Group", 
                        value = "Posterior Probability (%)", 
                        Healthy:SRC)

#convert
posterior_df2$`Posterior Probability (%)` <-
posterior_df2$`Posterior Probability (%)` * 100

by(posterior_df2$`Posterior Probability (%)`[posterior_df2$Sex =="Male"],
posterior_df2$Group[posterior_df2$Sex=="Male"], median)

###Manuscript Figure###
###Male posterior density plots
#save first
pdf("fig1.pdf", width = 12, height = 8, pointsize = 12)
#set multiplot
par(mfrow = c(1, 2))
#Males
dens(x = posterior_df2$`Posterior Probability (%)`[posterior_df2$Sex == "Male" &
 posterior_df2$Group == "Healthy"],
     xlim = c(0, 100), ylim = c(0, 0.1),
     lwd = 3, col = 3, xlab = "",
     font.lab = 2, cex.lab = 0.9, main = "Males",
     frame = FALSE)
dens(x = posterior_df2$`Posterior Probability (%)`[posterior_df2$Sex == "Male" &
      posterior_df2$Group == "MSI"],
     lwd = 3, col = 4, add = TRUE)
dens(x = posterior_df2$`Posterior Probability (%)`[posterior_df2$Sex == "Male" &
     posterior_df2$Group == "SRC"],
     lwd = 3, col = 2, add = TRUE)
dens(posterior_prior_df2$`Prior Probability (%)`[posterior_prior_df2$Sex == "Male" &
     posterior_prior_df2$Group == "Healthy"],
     lwd = 1, lty = 6, col = 'gray30', add = TRUE)
dens(posterior_prior_df2$`Prior Probability (%)`[posterior_prior_df2$Sex == "Male" &
      posterior_prior_df2$Group == "MSI"],
     lwd = 1, lty = 6, col = 'gray60', add = TRUE)
dens(posterior_prior_df2$`Prior Probability (%)`[posterior_prior_df2$Sex == "Male" &
     posterior_prior_df2$Group == "SRC"],
     lwd = 1, lty = 6, col = 'gray90', add = TRUE)

legend(0, 0.1, box.lty = 0, cex = 1,
       legend = c("SRC", "MSI", "Healthy"),
       fill = c("red", "blue", "green"))
text(20, 0.008, 'SRC Prior', cex = 0.7, col = 'gray30')
text(20, 0.006, 'MSI Prior', cex = 0.7, col = 'gray60')
text(20, 0.004, 'Healthy Prior', cex = 0.7, col = 'gray90')
text(-1, 0.1, "A", cex = 0.9, font = 2)

dens(x = posterior_df2$`Posterior Probability (%)`[posterior_df2$Sex == "Female" &
 posterior_df2$Group == "Healthy"],
     xlim = c(0, 100), ylim = c(0, 0.1),
     lwd = 3, col = 3, xlab = "", ylab = "",
     font.lab = 2, cex.lab = 0.9, main = "Females",
     frame = FALSE)
text(-1, 0.1, "B", cex = 0.9, font = 2)     
dens(x = posterior_df2$`Posterior Probability (%)`[posterior_df2$Sex == "Female" &
      posterior_df2$Group == "MSI"],
     lwd = 3, col = 4, add = TRUE)
dens(x = posterior_df2$`Posterior Probability (%)`[posterior_df2$Sex == "Female" &
     posterior_df2$Group == "SRC"],
     lwd = 3, col = 2, add = TRUE)
dens(posterior_prior_df2$`Prior Probability (%)`[posterior_prior_df2$Sex == "Female" &
     posterior_prior_df2$Group == "Healthy"],
     lwd = 1, lty = 6, col = 'gray30', add = TRUE)
dens(posterior_prior_df2$`Prior Probability (%)`[posterior_prior_df2$Sex == "Female" &
      posterior_prior_df2$Group == "MSI"],
     lwd = 1, lty = 6, col = 'gray60', add = TRUE)
dens(posterior_prior_df2$`Prior Probability (%)`[posterior_prior_df2$Sex == "Female" &
     posterior_prior_df2$Group == "SRC"],
     lwd = 1, lty = 6, col = 'gray90', add = TRUE)

mtext("Probability of Future MSI (%)", side = 3, line = -39, outer = TRUE, font = 2, cex = 1)     

fig1 <- recordPlot()
dev.off()

#male summary
precis(posterior_sum_df[-c(4)][posterior_sum_df$Sex == "Male", ], probs = 0.9)

#female summary
precis(posterior_sum_df[-c(4)][posterior_sum_df$Sex == "Female", ], probs = 0.9)

#pprob masses
#male
apply(posterior_sum_df[c(5:7)][posterior_sum_df$Sex == "Male", ],
2, function(x) mean(x > 0))

#female
apply(posterior_sum_df[c(5:7)][posterior_sum_df$Sex == "Female", ],
2, function(x) mean(x > 0))

##calculating proportional changes in the outcomes by converting to log odds and then finding
##the exponent of the difference

#create new df
posterior_sum_df2 <- posterior_sum_df[c(1:4)]
#/100 to get to a probability
posterior_sum_df2[c(1:3)] <- posterior_sum_df2[c(1:3)] / 100
#convert to log odds
posterior_sum_df2[c(1:3)] <- apply(posterior_sum_df2[c(1:3)], 2, function(x) log(x / (1 - x)))

#compare log odds
posterior_sum_df2$SRC_MSI <- exp(posterior_sum_df2$SRC - posterior_sum_df2$MSI)
posterior_sum_df2$SRC_Healthy <- exp(posterior_sum_df2$SRC - posterior_sum_df2$Healthy)
posterior_sum_df2$MSI_Healthy <- exp(posterior_sum_df2$MSI - posterior_sum_df2$Healthy)

#summarize
#male summary
precis(posterior_sum_df2[-c(4)][posterior_sum_df2$Sex == "Male", ], probs = 0.9)

#female summary
precis(posterior_sum_df2[-c(4)][posterior_sum_df2$Sex == "Female", ], probs = 0.9)

#pprob masses
#male
apply(posterior_sum_df2[c(5:7)][posterior_sum_df$Sex == "Male", ],
2, function(x) mean(x > 1)) #compared to 1, because its relative risk

#female
apply(posterior_sum_df2[c(5:7)][posterior_sum_df$Sex == "Female", ],
2, function(x) mean(x > 1))

