# Script to plot community-level responses for PA and Normal models

comm_level <- readRDS(here::here("output/Comm_responses.rds"))
length(comm_level)

pred.Temp <- comm_level[[1]]
pred.Forest <-comm_level[[2]]
pred.Grass <- comm_level[[3]]
pred.Crops <- comm_level[[4]]

# PA model - Focal: temperature
r.all <- list()
for(i in 1:length(pred.Temp$PA)){
  r.all[[i]] <- as.character(pred.Temp$PA[[i]]$labels$caption)
}
PA.cor.values <- matrix(as.numeric(do.call(rbind, r.all)[,3]), ncol = 2, byrow = F)
colnames(PA.cor.values) <- c("Richness.PA", "CWM.PA")
rownames(PA.cor.values) <- c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2")
PA.cor.values

df.all <- cbind(rbind(pred.Temp$PA[[1]]$data, pred.Temp$PA[[2]]$data,
                      pred.Temp$PA[[3]]$data, pred.Temp$PA[[4]]$data,
                      pred.Temp$PA[[5]]$data, pred.Temp$PA[[6]]$data, 
                      pred.Temp$PA[[7]]$data),
                Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", PA.cor.values[1,1], sep = ""),
             paste("r = ", PA.cor.values[2,1], sep = ""),
             paste("r = ", PA.cor.values[3,1], sep = ""),
             paste("r = ", PA.cor.values[4,1], sep = ""),
             paste("r = ", PA.cor.values[5,1], sep = ""),
             paste("r = ", PA.cor.values[6,1], sep = ""),
             paste("r = ", PA.cor.values[7,1], sep = ""))

PA.p1 <- ggplot(df.all, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Comm)) + 
  labs(tag = "a)", x = "Temperature", y = "Species Richness") +
  theme(legend.position = c(.15, .25), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

df.all.WTR <- cbind(rbind(pred.Temp$PA[[1+7]]$data, pred.Temp$PA[[2+7]]$data,
                          pred.Temp$PA[[3+7]]$data, pred.Temp$PA[[4+7]]$data,
                          pred.Temp$PA[[5+7]]$data, pred.Temp$PA[[6+7]]$data, 
                          pred.Temp$PA[[7+7]]$data),
                    Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", PA.cor.values[1,2], sep = ""),
             paste("r = ", PA.cor.values[2,2], sep = ""),
             paste("r = ", PA.cor.values[3,2], sep = ""),
             paste("r = ", PA.cor.values[4,2], sep = ""),
             paste("r = ", PA.cor.values[5,2], sep = ""),
             paste("r = ", PA.cor.values[6,2], sep = ""),
             paste("r = ", PA.cor.values[7,2], sep = ""))

PA.p2 <- ggplot(df.all.WTR, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = .2, aes(fill = Comm)) +
  labs(tag = "b)", x = "Temperature", y = "CWM-WTR") +
  theme(legend.position = c(.15, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

# Normal model - Focal: Temperature
r.all <- list()
for(i in 1:length(pred.Temp$Normal)){
  r.all[[i]] <- as.character(pred.Temp$Normal[[i]]$labels$caption)
}
Nor.cor.values <- matrix(as.numeric(do.call(rbind, r.all)[,3]), ncol = 2, byrow = F)
colnames(Nor.cor.values) <- c("Richness.PA", "CWM.PA")
rownames(Nor.cor.values) <- c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2")
Nor.cor.values

df.all <- cbind(rbind(pred.Temp$Normal[[1]]$data, pred.Temp$Normal[[2]]$data,
                      pred.Temp$Normal[[3]]$data, pred.Temp$Normal[[4]]$data,
                      pred.Temp$Normal[[5]]$data, pred.Temp$Normal[[6]]$data, 
                      pred.Temp$Normal[[7]]$data),
                Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", Nor.cor.values[1,1], sep = ""),
             paste("r = ", Nor.cor.values[2,1], sep = ""),
             paste("r = ", Nor.cor.values[3,1], sep = ""),
             paste("r = ", Nor.cor.values[4,1], sep = ""),
             paste("r = ", Nor.cor.values[5,1], sep = ""),
             paste("r = ", Nor.cor.values[6,1], sep = ""),
             paste("r = ", Nor.cor.values[7,1], sep = ""))

Nor.p1 <- ggplot(df.all, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Comm)) + 
  labs(tag = "a)", x = "Temperature", y = "Summed response") +
  theme(legend.position = c(.15, .25), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

df.all.WTR <- cbind(rbind(pred.Temp$Normal[[1+7]]$data, pred.Temp$Normal[[2+7]]$data,
                          pred.Temp$Normal[[3+7]]$data, pred.Temp$Normal[[4+7]]$data,
                          pred.Temp$Normal[[5+7]]$data, pred.Temp$Normal[[6+7]]$data, 
                          pred.Temp$Normal[[7+7]]$data),
                    Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", Nor.cor.values[1,2], sep = ""),
             paste("r = ", Nor.cor.values[2,2], sep = ""),
             paste("r = ", Nor.cor.values[3,2], sep = ""),
             paste("r = ", Nor.cor.values[4,2], sep = ""),
             paste("r = ", Nor.cor.values[5,2], sep = ""),
             paste("r = ", Nor.cor.values[6,2], sep = ""),
             paste("r = ", Nor.cor.values[7,2], sep = ""))

Nor.p2 <- ggplot(df.all.WTR, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = .2, aes(fill = Comm)) +
  labs(tag = "b)", x = "Temperature", y = "CWM-WTR") +
  theme(legend.position = c(.15, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

PA.temp <- cowplot::plot_grid(PA.p1, PA.p2, ncol = 1)

Nor.temp <- cowplot::plot_grid(Nor.p1, Nor.p2, ncol = 1)

# PA model - Focal: Forest
r.all <- list()
for(i in 1:length(pred.Forest$PA)){
  r.all[[i]] <- as.character(pred.Forest$PA[[i]]$labels$caption)
}
PA.cor.values <- matrix(as.numeric(do.call(rbind, r.all)[,3]), ncol = 2, byrow = F)
colnames(PA.cor.values) <- c("Richness.PA", "CWM.PA")
rownames(PA.cor.values) <- c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2")
PA.cor.values

df.all <- cbind(rbind(pred.Forest$PA[[1]]$data, pred.Forest$PA[[2]]$data,
                      pred.Forest$PA[[3]]$data, pred.Forest$PA[[4]]$data,
                      pred.Forest$PA[[5]]$data, pred.Forest$PA[[6]]$data, 
                      pred.Forest$PA[[7]]$data),
                Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", PA.cor.values[1,1], sep = ""),
             paste("r = ", PA.cor.values[2,1], sep = ""),
             paste("r = ", PA.cor.values[3,1], sep = ""),
             paste("r = ", PA.cor.values[4,1], sep = ""),
             paste("r = ", PA.cor.values[5,1], sep = ""),
             paste("r = ", PA.cor.values[6,1], sep = ""),
             paste("r = ", PA.cor.values[7,1], sep = ""))

PA.p3 <- ggplot(df.all, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Comm)) + 
  labs(tag = "c)", x = "Forest", y = "Species Richness") +
  theme(legend.position = c(.8, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

df.all.WTR <- cbind(rbind(pred.Forest$PA[[1+7]]$data, pred.Forest$PA[[2+7]]$data,
                          pred.Forest$PA[[3+7]]$data, pred.Forest$PA[[4+7]]$data,
                          pred.Forest$PA[[5+7]]$data, pred.Forest$PA[[6+7]]$data, 
                          pred.Forest$PA[[7+7]]$data),
                    Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", PA.cor.values[1,2], sep = ""),
             paste("r = ", PA.cor.values[2,2], sep = ""),
             paste("r = ", PA.cor.values[3,2], sep = ""),
             paste("r = ", PA.cor.values[4,2], sep = ""),
             paste("r = ", PA.cor.values[5,2], sep = ""),
             paste("r = ", PA.cor.values[6,2], sep = ""),
             paste("r = ", PA.cor.values[7,2], sep = ""))

PA.p4 <- ggplot(df.all.WTR, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = .2, aes(fill = Comm)) +
  labs(tag = "d)", x = "Forest", y = "CWM-WTR") +
  theme(legend.position = c(.8, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

# Normal model - Focal: Forest
r.all <- list()
for(i in 1:length(pred.Forest$Normal)){
  r.all[[i]] <- as.character(pred.Forest$Normal[[i]]$labels$caption)
}
Nor.cor.values <- matrix(as.numeric(do.call(rbind, r.all)[,3]), ncol = 2, byrow = F)
colnames(Nor.cor.values) <- c("Richness.PA", "CWM.PA")
rownames(Nor.cor.values) <- c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2")
Nor.cor.values

df.all <- cbind(rbind(pred.Forest$Normal[[1]]$data, pred.Forest$Normal[[2]]$data,
                      pred.Forest$Normal[[3]]$data, pred.Forest$Normal[[4]]$data,
                      pred.Forest$Normal[[5]]$data, pred.Forest$Normal[[6]]$data, 
                      pred.Forest$Normal[[7]]$data),
                Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", Nor.cor.values[1,1], sep = ""),
             paste("r = ", Nor.cor.values[2,1], sep = ""),
             paste("r = ", Nor.cor.values[3,1], sep = ""),
             paste("r = ", Nor.cor.values[4,1], sep = ""),
             paste("r = ", Nor.cor.values[5,1], sep = ""),
             paste("r = ", Nor.cor.values[6,1], sep = ""),
             paste("r = ", Nor.cor.values[7,1], sep = ""))

Nor.p3 <- ggplot(df.all, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Comm)) + 
  labs(tag = "c)", x = "Forest", y = "Summed response") +
  theme(legend.position = c(.8, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

df.all.WTR <- cbind(rbind(pred.Forest$Normal[[1+7]]$data, pred.Forest$Normal[[2+7]]$data,
                          pred.Forest$Normal[[3+7]]$data, pred.Forest$Normal[[4+7]]$data,
                          pred.Forest$Normal[[5+7]]$data, pred.Forest$Normal[[6+7]]$data, 
                          pred.Forest$Normal[[7+7]]$data),
                    Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", Nor.cor.values[1,2], sep = ""),
             paste("r = ", Nor.cor.values[2,2], sep = ""),
             paste("r = ", Nor.cor.values[3,2], sep = ""),
             paste("r = ", Nor.cor.values[4,2], sep = ""),
             paste("r = ", Nor.cor.values[5,2], sep = ""),
             paste("r = ", Nor.cor.values[6,2], sep = ""),
             paste("r = ", Nor.cor.values[7,2], sep = ""))

Nor.p4 <- ggplot(df.all.WTR, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = .2, aes(fill = Comm)) +
  labs(tag = "d)", x = "Forest", y = "CWM-WTR") +
  theme(legend.position = c(.8, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

PA.forest <- cowplot::plot_grid(PA.p3, PA.p4, ncol = 1)

Nor.forest <- cowplot::plot_grid(Nor.p3, Nor.p4, ncol = 1)

# PA model - Focal: Grasslands
r.all <- list()
for(i in 1:length(pred.Grass$PA)){
  r.all[[i]] <- as.character(pred.Grass$PA[[i]]$labels$caption)
}
PA.cor.values <- matrix(as.numeric(do.call(rbind, r.all)[,3]), ncol = 2, byrow = F)
colnames(PA.cor.values) <- c("Richness.PA", "CWM.PA")
rownames(PA.cor.values) <- c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2")
PA.cor.values

df.all <- cbind(rbind(pred.Grass$PA[[1]]$data, pred.Grass$PA[[2]]$data,
                      pred.Grass$PA[[3]]$data, pred.Grass$PA[[4]]$data,
                      pred.Grass$PA[[5]]$data, pred.Grass$PA[[6]]$data, 
                      pred.Grass$PA[[7]]$data),
                Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", PA.cor.values[1,1], sep = ""),
             paste("r = ", PA.cor.values[2,1], sep = ""),
             paste("r = ", PA.cor.values[3,1], sep = ""),
             paste("r = ", PA.cor.values[4,1], sep = ""),
             paste("r = ", PA.cor.values[5,1], sep = ""),
             paste("r = ", PA.cor.values[6,1], sep = ""),
             paste("r = ", PA.cor.values[7,1], sep = ""))

PA.p5 <- ggplot(df.all, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Comm)) + 
  labs(tag = "e)", x = "Grassland", y = "Species Richness") +
  theme(legend.position = c(.15, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

df.all.WTR <- cbind(rbind(pred.Grass$PA[[1+7]]$data, pred.Grass$PA[[2+7]]$data,
                          pred.Grass$PA[[3+7]]$data, pred.Grass$PA[[4+7]]$data,
                          pred.Grass$PA[[5+7]]$data, pred.Grass$PA[[6+7]]$data, 
                          pred.Grass$PA[[7+7]]$data),
                    Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", PA.cor.values[1,2], sep = ""),
             paste("r = ", PA.cor.values[2,2], sep = ""),
             paste("r = ", PA.cor.values[3,2], sep = ""),
             paste("r = ", PA.cor.values[4,2], sep = ""),
             paste("r = ", PA.cor.values[5,2], sep = ""),
             paste("r = ", PA.cor.values[6,2], sep = ""),
             paste("r = ", PA.cor.values[7,2], sep = ""))

PA.p6 <- ggplot(df.all.WTR, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = .2, aes(fill = Comm)) +
  labs(tag = "f)", x = "Grassland", y = "CWM-WTR") +
  theme(legend.position = c(.15, .25), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

# Normal model - Focal: Grassland
r.all <- list()
for(i in 1:length(pred.Grass$Normal)){
  r.all[[i]] <- as.character(pred.Grass$Normal[[i]]$labels$caption)
}
Nor.cor.values <- matrix(as.numeric(do.call(rbind, r.all)[,3]), ncol = 2, byrow = F)
colnames(Nor.cor.values) <- c("Richness.PA", "CWM.PA")
rownames(Nor.cor.values) <- c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2")
Nor.cor.values

df.all <- cbind(rbind(pred.Grass$Normal[[1]]$data, pred.Grass$Normal[[2]]$data,
                      pred.Grass$Normal[[3]]$data, pred.Grass$Normal[[4]]$data,
                      pred.Grass$Normal[[5]]$data, pred.Grass$Normal[[6]]$data, 
                      pred.Grass$Normal[[7]]$data),
                Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", Nor.cor.values[1,1], sep = ""),
             paste("r = ", Nor.cor.values[2,1], sep = ""),
             paste("r = ", Nor.cor.values[3,1], sep = ""),
             paste("r = ", Nor.cor.values[4,1], sep = ""),
             paste("r = ", Nor.cor.values[5,1], sep = ""),
             paste("r = ", Nor.cor.values[6,1], sep = ""),
             paste("r = ", Nor.cor.values[7,1], sep = ""))

Nor.p5 <- ggplot(df.all, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Comm)) + 
  labs(tag = "e)", x = "Grassland", y = "Summed response") +
  theme(legend.position = c(.2, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

df.all.WTR <- cbind(rbind(pred.Grass$Normal[[1+7]]$data, pred.Grass$Normal[[2+7]]$data,
                          pred.Grass$Normal[[3+7]]$data, pred.Grass$Normal[[4+7]]$data,
                          pred.Grass$Normal[[5+7]]$data, pred.Grass$Normal[[6+7]]$data, 
                          pred.Grass$Normal[[7+7]]$data),
                    Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", Nor.cor.values[1,2], sep = ""),
             paste("r = ", Nor.cor.values[2,2], sep = ""),
             paste("r = ", Nor.cor.values[3,2], sep = ""),
             paste("r = ", Nor.cor.values[4,2], sep = ""),
             paste("r = ", Nor.cor.values[5,2], sep = ""),
             paste("r = ", Nor.cor.values[6,2], sep = ""),
             paste("r = ", Nor.cor.values[7,2], sep = ""))

Nor.p6 <- ggplot(df.all.WTR, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = .2, aes(fill = Comm)) +
  labs(tag = "f)", x = "Grassland", y = "CWM-WTR") +
  theme(legend.position = c(.2, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

PA.grass <- cowplot::plot_grid(PA.p5, PA.p6, ncol = 1)

Nor.grass <- cowplot::plot_grid(Nor.p5, Nor.p6, ncol = 1)

# PA model - Focal: Crops
r.all <- list()
for(i in 1:length(pred.Crops$PA)){
  r.all[[i]] <- as.character(pred.Crops$PA[[i]]$labels$caption)
}
PA.cor.values <- matrix(as.numeric(do.call(rbind, r.all)[,3]), ncol = 2, byrow = F)
colnames(PA.cor.values) <- c("Richness.PA", "CWM.PA")
rownames(PA.cor.values) <- c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2")
PA.cor.values

df.all <- cbind(rbind(pred.Crops$PA[[1]]$data, pred.Crops$PA[[2]]$data,
                      pred.Crops$PA[[3]]$data, pred.Crops$PA[[4]]$data,
                      pred.Crops$PA[[5]]$data, pred.Crops$PA[[6]]$data, 
                      pred.Crops$PA[[7]]$data),
                Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", PA.cor.values[1,1], sep = ""),
             paste("r = ", PA.cor.values[2,1], sep = ""),
             paste("r = ", PA.cor.values[3,1], sep = ""),
             paste("r = ", PA.cor.values[4,1], sep = ""),
             paste("r = ", PA.cor.values[5,1], sep = ""),
             paste("r = ", PA.cor.values[6,1], sep = ""),
             paste("r = ", PA.cor.values[7,1], sep = ""))

PA.p7 <- ggplot(df.all, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Comm)) + 
  labs(tag = "g)", x = "Crops", y = "Species Richness") +
  theme(legend.position = c(.15, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

df.all.WTR <- cbind(rbind(pred.Crops$PA[[1+7]]$data, pred.Crops$PA[[2+7]]$data,
                          pred.Crops$PA[[3+7]]$data, pred.Crops$PA[[4+7]]$data,
                          pred.Crops$PA[[5+7]]$data, pred.Crops$PA[[6+7]]$data, 
                          pred.Crops$PA[[7+7]]$data),
                    Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", PA.cor.values[1,2], sep = ""),
             paste("r = ", PA.cor.values[2,2], sep = ""),
             paste("r = ", PA.cor.values[3,2], sep = ""),
             paste("r = ", PA.cor.values[4,2], sep = ""),
             paste("r = ", PA.cor.values[5,2], sep = ""),
             paste("r = ", PA.cor.values[6,2], sep = ""),
             paste("r = ", PA.cor.values[7,2], sep = ""))

PA.p8 <- ggplot(df.all.WTR, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = .2, aes(fill = Comm)) +
  labs(tag = "h)", x = "Crops", y = "CWM-WTR") +
  theme(legend.position = c(.8, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

# Normal model - Focal: Crops
r.all <- list()
for(i in 1:length(pred.Crops$Normal)){
  r.all[[i]] <- as.character(pred.Crops$Normal[[i]]$labels$caption)
}
Nor.cor.values <- matrix(as.numeric(do.call(rbind, r.all)[,3]), ncol = 2, byrow = F)
colnames(Nor.cor.values) <- c("Richness.PA", "CWM.PA")
rownames(Nor.cor.values) <- c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2")
Nor.cor.values

df.all <- cbind(rbind(pred.Crops$Normal[[1]]$data, pred.Crops$Normal[[2]]$data,
                      pred.Crops$Normal[[3]]$data, pred.Crops$Normal[[4]]$data,
                      pred.Crops$Normal[[5]]$data, pred.Crops$Normal[[6]]$data, 
                      pred.Crops$Normal[[7]]$data),
                Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", Nor.cor.values[1,1], sep = ""),
             paste("r = ", Nor.cor.values[2,1], sep = ""),
             paste("r = ", Nor.cor.values[3,1], sep = ""),
             paste("r = ", Nor.cor.values[4,1], sep = ""),
             paste("r = ", Nor.cor.values[5,1], sep = ""),
             paste("r = ", Nor.cor.values[6,1], sep = ""),
             paste("r = ", Nor.cor.values[7,1], sep = ""))

Nor.p7 <- ggplot(df.all, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = 0.2, aes(fill = Comm)) + 
  labs(tag = "g)", x = "Crops", y = "Summed response") +
  theme(legend.position = c(.85, .3), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

df.all.WTR <- cbind(rbind(pred.Crops$Normal[[1+7]]$data, pred.Crops$Normal[[2+7]]$data,
                          pred.Crops$Normal[[3+7]]$data, pred.Crops$Normal[[4+7]]$data,
                          pred.Crops$Normal[[5+7]]$data, pred.Crops$Normal[[6+7]]$data, 
                          pred.Crops$Normal[[7+7]]$data),
                    Comm = rep(c("EEA2", "JG1", "JG2", "QR1", "QR2", "SG1", "SG2"), each = 20))

new.lab <- c(paste("r = ", Nor.cor.values[1,2], sep = ""),
             paste("r = ", Nor.cor.values[2,2], sep = ""),
             paste("r = ", Nor.cor.values[3,2], sep = ""),
             paste("r = ", Nor.cor.values[4,2], sep = ""),
             paste("r = ", Nor.cor.values[5,2], sep = ""),
             paste("r = ", Nor.cor.values[6,2], sep = ""),
             paste("r = ", Nor.cor.values[7,2], sep = ""))

Nor.p8 <- ggplot(df.all.WTR, aes(x = xx, y = v1, colour =  Comm)) +
  geom_smooth(method = "lm", alpha = .2, aes(fill = Comm)) +
  labs(tag = "h)", x = "Crops", y = "CWM-WTR") +
  theme(legend.position = c(.8, .8), legend.text = element_text(size = 6),
        legend.background = element_blank(), legend.key.size = unit(1, "lines")) +
  scale_color_manual(values = cbp, labels = new.lab, name = NULL) + 
  scale_fill_manual(values = cbp, labels = new.lab, name = NULL)

PA.crop <- cowplot::plot_grid(PA.p7, PA.p8, ncol = 1)

Nor.crop <- cowplot::plot_grid(Nor.p7, Nor.p8, ncol = 1)

leg <- cowplot::get_legend(PA.p1 + theme(legend.position = "right") + 
                             scale_color_manual(values = cbp, labels = unique(df.all$Comm), name = "Comm") + 
                             scale_fill_manual(values = cbp, labels = unique(df.all$Comm), name = "Comm"))

pred.PA <- cowplot::plot_grid(PA.temp, PA.forest, PA.grass, PA.crop, leg, ncol = 5,
                              rel_widths = c(2,2,2,2,1))

save_plot(here::here("output/figures/Final_figures/Fig4_PA_Pred_response.png"), pred.PA, 
          base_height = 8, base_width = 16)

pred.Nor <- cowplot::plot_grid(Nor.temp, Nor.forest, Nor.grass, Nor.crop, leg, ncol = 5,
                               rel_widths = c(2,2,2,2,1))

save_plot(here::here("output/figures/Final_figures/Fig5_Normal_Pred_response.png"), pred.Nor, 
          base_height = 8, base_width = 14)
