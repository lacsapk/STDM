# 1. Datensatz ----

## 1.1 Import ----

library(data.table)

data <- fread("Breast_GSE45827.csv")

## 1.3 Beschreib ----

'Datensatz auf Kaggle gefunden:
- https://www.kaggle.com/datasets/brunogrisci/breast-cancer-gene-expression-cumida/data

Originale Quelle vom Datensatz:
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45827

Der Datensatz GSE45827 enthält 151 Zellbeobachtungen mit jeweils 54.676 Genen, die 
in sechs Klassen unterteilt sind und der Untersuchung von Brustkrebs-Subtypen dienen.

Klassen:
- Triple-Negative (TN) – 41 Proben (basal)
- HER2+ – 30 Proben (HER)
- Luminal B – 30 Proben (luminal_B)
- Luminal A – 29 Proben (luminal_A)
- Normales Gewebe – 7 Proben (normal)
- Zelllinien – 14 Proben (cell_line)

Brustkrebs-subtypen die sich in Therapieansprechen und prognostizierten Aggressivitätsgrad unterscheiden.'

table(data$type)
str(data)

## 1.3 Optimierung ----

'Datensatz hat zu viele Variablen, weshalb wir nur die 2000 Gene mit der höchsten Varianz behalten'

data <- as.data.frame(data)

var <- apply(data[,-c(1,2)], 2, var) # Varianz der Spalten (Gene)
sorted.var <-  sort(var, decreasing = T)
keep <- names(sorted.var)

dat <- data[ ,c("type", keep)]

str(dat)

## 1.4 Skalieren? ----

# boxplot(dat[,-1]) # nichts erkennbar (Auslassen)

vars <- apply(dat[,-1], 2, var)

hist(vars, breaks = 100, xlim = c(0, 4),
     main = "Verteilung Varianzen",
     xlab = "Varianz")
abline(v = mean(vars), col = "red", lwd = 1.5, lty = 2)
legend("topright",
       legend = "Mittelwert",
       lty = 2,
       lwd = 1.5,
       col = "red",
       bty = "n")

### 1.4.1 nein ----
'Gene wurden in derselben Einheit gemessen, was gegen eine Skalierung sprechen würde.
Auch ist die Varianz der einzelnen Variablen ziemlich ähnlich'


### 1.4.2 ja ----
'Wenige Gene weisen eine viel grössere Varianz auf als andere.
Zudem sei es weit verbreitet Gen-Daten zu skalieren (standardisieren):
- https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/#:~:text=In%20principal%20component%20analysis%2C%20variables%20are%20often%20scaled,the%20PCA%20outputs%20obtained%20will%20be%20severely%20affected.
SW01 40/51: Bei unsicherheit ist es besser zu standardisieren'

# 2. PCA ----

## 2.1 prcomp() ----

pca <- prcomp(dat[,-1], scale. = T)
pca <- prcomp(dat[,-1], scale. = F) # Ohne Skalierung Übersichtlicher
                          
### 2.1.1 Anteil der Varianz ----

var <- pca$sdev^2
cum <- cumsum(var)/sum(var)

plot(1:length(cum), cum, type = "l",
     main = "Varianz erklärt durch PC",
     xlab = "PCx", ylab = "Anteil Varianz",
     col = "red", lwd = 1.5)
abline(h = 0.8, v = min(which(cum >= 0.8)), col = "grey", lty = 2, lwd = 1.5)
axis(1, at = min(which(cum >= 0.8)), labels = min(which(cum >= 0.8)))
grid()

'um 80% der Varianz zu erklären müsste man die ersten 58 Hauptkomponenten verwenden.'

### 2.1.2 Welches Gene trägt am meisten zu PC1 bei ----

rot <- abs(pca$rotation[,1]) / sum(abs(pca$rotation[,1])) # % Anteil rotation
rot <- sort(rot, decreasing = T)

toprot <- rot[1:10] # Top 10 Gene

barplot(toprot,
        las = 2,
        main = "Top 10 Gene in PC1",
        ylim = c(0, 0.00035))
grid()

### Plot base R

levels(as.factor(dat[,1]))
colours = c("red", "lightblue", "darkred", "yellow", "orange", "green")

plot(pca$x[,1], pca$x[,2], pch = 20, 
     main = "PCA", xlab = "PC1", ylab = "PC2",
     col = colours[as.factor(dat[,1])])
legend("topleft", 
       legend = levels(as.factor(dat[,1])),
       col = colours,
       pch = 20,
       cex = 0.8,
       pt.cex = 1,
       bty = "n")

### Plot ggplot

library(ggplot2)
library(ggthemes)

ggplot(pca$x, aes(PC1, PC2, colour = dat[,1])) +
  geom_point() +
  scale_color_manual(values = colours) +
  labs(colour = "Tumor type",
        title = "PCA") +
  theme_bw() 

## 2.2 PcaHuber() ----

library(MASS)
library(rrcov)

set.seed(2)
pca_h <- PcaHubert(dat[,-1], scale = T, k = 2)
pca_h <- PcaHubert(dat[,-1], scale = F, k = 2)

### Plot base R

plot(pca_h@scores[,1], pca_h$scores[,2], pch = 20, 
     main = "PCA_Hubert", xlab = "PC1", ylab = "PC2",
     col = colours[as.factor(dat[,1])])
legend("bottomright", 
       legend = levels(as.factor(dat[,1])),
       col = colours,
       pch = 20,
       cex = 0.8,
       pt.cex = 1,
       bty = "n")

### Plot ggplot

ggplot(pca_h@scores, aes(PC1, PC2, colour = dat[,1])) +
  geom_point() +
  scale_color_manual(values = colours) +
  labs(colour = "Tumor type",
        title = "PCA_Hubert") +
  theme_bw() 

'Robuste PCA zeigt, dass normale Proben eine homogene, von Tumoren klar 
unterscheidbare Gruppe bilden und dass diese Abgrenzung erst nach 
Eliminierung von Verzerrungen (durch Ausreisser/Varianzunterschiede) 
sichtbar wird.'

### 2.2.1 Orthogonale Distanz ----
n <- nrow(dat)
plot(x = 1:n, y = pca_h@od, pch = 20, col = colours[as.factor(dat[,1])])
abline(h = pca_h@cutoff.od)
legend("topright", 
       legend = levels(as.factor(dat[,1])),
       col = colours,
       pch = 20,
       cex = 0.8,
       pt.cex = 1,
       bty = "n")

'Die orthogonale Distanz zeigt, dass sich die cell_line Zellen stark von den
anderen Zelltypen unterscheiden. In der robusten PCA-Darstellung mit den ersten
beiden Hauptkomponenten wird dieser Unterschied jedoch nicht klar hervorgehoben
bzw. sie werden schlecht durch die PCA abgebildet. Ihre Variation könnte
hauptsächlich in anderen Dimensionen liegen.'


### 2.2.2 Welches Gene hat am meisten einfluss auf PC1----

options(scipen = 999) # Dezimal als ylim

anteilload <- abs(pca_h@loadings[,1]) / sum(abs(pca_h@loadings[,1]))
anteilload <- sort(anteilload, decreasing = T)

topload <- anteilload[1:10]

barplot(topload,
        las = 2,
        main = "Top 10 Genes in PC1",
        ylim = c(0, 0.00045))
grid()

### 2.3 Ergebnis ----

'Cell_line sind Zellen aus dem labor, welche zur Forschung verwendet werden.
Sie können nützliche Modelle sein, das klassischen PCA lässt vermuten, dass diese 
nicht gleich wie Patientenproben sind.

Die orthogonale Distanz zeigt, dass sich die cell_line Zellen stark von den
anderen Zelltypen unterscheiden. In der robusten PCA-Darstellung mit den ersten
beiden Hauptkomponenten wird dieser Unterschied jedoch nicht klar hervorgehoben
bzw. sie werden schlecht durch die PCA abgebildet. Ihre Variation könnte
hauptsächlich in anderen Dimensionen liegen.

Wir Vermuten, dass die Cell_line in den Genen viele extreme Werte beinhaltet,
was dazu geführt hat, dass der Unterschied zwischen normalen und Krebszellen beim
klassischen PCA untergeht. Bei PcaHubert haben diese Ausreisser durch eine robustere
Berechnung der Varianz weniger grossen Einfluss und die Abbildung ist somit representativer.

Bei allen Abbildungen sind Gruppen zu erkennen.'

# 3. UMAP ----

library(umap)

umap <- umap(dat[,-1], input = "data",
             n_neighbors = 7,         # Mit neighbors zwischen 5 und 10 i.O
             min_dist = 0.5,
             random_state = 29)

## Plot mit base R

plot(umap$layout[,1], umap$layout[,2], pch = 20,
     main = "UMAP", xlab = "Umap1", ylab = "Umap2",
     col = colours[as.factor(dat$type)])

legend("topright", legend = levels(as.factor(dat[,1])),
       pch = 20,
       col = colours,
       cex = 0.8,
       pt.cex = 1,
       bty = "n")

## Plot mit ggplot

colnames(umap$layout) <- c("Umap1", "Umap2")

ggplot(umap$layout, aes(Umap1, Umap2, colour = dat[,1])) +
  geom_point() +
  scale_color_manual(values = colours) +
  labs(colour = "tumor type",
       title = "UMAP") +
  theme_bw()

## 3.1 Ergebnis----
'
Es sind 3 klare Cluster zu erkennen, wobei man beim Cluster der Tumore
darüber diskutieren könnte ob es sich um einzele Gruppen handelt. 
'

# 4. T-Sne ----

library(Rtsne)

pca <- prcomp(dat[,-1], scale = F )
pca <- pca$x[,1:151]

set.seed(120)
tsne <- Rtsne(pca, dims = 2, is_distance = F,
              pca = F,  # Warum funktioniert es nicht bei pca = T?
              perplexity = 20, # Zwischen 10 und 20 ok
              max_iter = 1000)

tsne <- as.data.frame(tsne$Y)

# Plot

plot(tsne$V1, tsne$V2, pch = 20, 
     col = colours[as.factor(dat[,1])], main = "T-SNE")
legend("bottomleft", legend = levels(as.factor(dat[,1])),
       bty = "n",
       cex = 0.8,
       col = colours,
       pch = 20)

## GGplot

library(ggplot2)

ggplot(tsne, aes(V1, V2, colour = dat[,1])) +
  geom_point() +
  scale_color_manual(values = colours) +
  labs(title = "T-SNE",
       colour = "Tumor Types") +
  theme_bw()

