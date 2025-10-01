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
keep <- names(sorted.var)[1:2000]

dat <- data[ ,c("type", keep)]

str(dat)

## 1.4 Skalieren? ----

boxplot(dat[,-1]) # nichts erkennbar

vars <- apply(dat[,-1], 2, var)

hist(vars, breaks = 100, xlim = c(1, 15),
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
'Gene wurden in derselben Einheit gemessen, was gegen eine Skalierung sprechen würde.'

### 1.4.2 ja ----
'Wenige Gene weisen eine viel grössere Varianz auf als andere.
Zudem sei es weit verbreitet Gen-Daten zu skalieren (standardisieren):
- https://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/112-pca-principal-component-analysis-essentials/#:~:text=In%20principal%20component%20analysis%2C%20variables%20are%20often%20scaled,the%20PCA%20outputs%20obtained%20will%20be%20severely%20affected.
SW01 40/51: Bei unsicherheit ist es besser zu standardisieren'

# 2. PCA ----

## 2.1 prcomp() ----

pca <- prcomp(dat[,-1], scale. = T)
pca <- prcomp(dat[,-1], scale. = F) # Vermutlich besser 38 PC's bis 0.8 erklärt 
                                    # statt 40 und plot minimal klarer
                          
### Anteil der Varianz

var <- pca$sdev^2
cum <- cumsum(var)/sum(var)

plot(1:length(cum), cum, type = "l",
     main = "Varianz erklärt durch PC",
     xlab = "PCx", ylab = "Anteil Varianz",
     col = "red", lwd = 1.5)
abline(h = 0.8, v = min(which(cum >= 0.8)), col = "grey", lty = 2, lwd = 1.5)
axis(1, at = min(which(cum >= 0.8)), labels = min(which(cum >= 0.8)))

'um 80% der Varianz zu erklären müsste man die ersten 40 Hauptkomponenten verwenden.'

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

### 2.3 Ergebnis ----

'Cell_line sind Zellen aus dem labor, welche zur Forschung verwendet werden.
Sie können nützliche Modelle sein, sind aber laut dem klassischen PCA
nicht gleich wie Patientenproben.

Wenn wir PcaHubert verwenden (robustere PCA-Methode mit Median), sieht es jedoch
so aus als sein die Cell_line doch sehr ähnlich zu den basal Zellen, dafür unterscheiden
sich die normalen Zellen grob von den Krebszellen.

Wir Vermuten, dass die Cell_line in den Genen viele Ausreisser beinhaltet,
was dazu geführt hat, dass der Unterschied zwischen normalen und Krebszellen beim
klassischen PCA untergeht. Bei PcaHubert haben diese Ausreisser durch den Median
einen weniger grossen Einfluss und die Abbildung ist somit representativer.

Bei allen Abbildungen sind Gruppen zu erkennen.'

# 3. UMAP ----

library(umap)

umap <- umap(dat[,-1], input = "data",
             n_neighbors = 5,         # Mit neighbors zwischen 5 und 10 i.O
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
'Es sind 4 klare Unterteilungen erkennbar. Diese sind normale Zellen,
Zelllinien aus dem Labor, Basal (Triple negative (TN)) und luminal_A, B und HER.
Obwohl letztere nahe beieinander sind, sind grenzen durch die Farben erkennbar.'

