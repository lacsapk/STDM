# 1. Datensatz ----

library(data.table)

data <- fread("Breast_GSE45827.csv")
str(data)

dat <- data[,-c(1,2)] # in dat sind keine Labels
str(dat)

# 2. Dimensionsreduktion für Visualisierung ----

levels(as.factor(data$type))
colours = c("red", "lightblue", "darkred", "yellow", "orange", "green")

library(umap)

umap <- umap(dat, input = "data",
             n_neighbors = 7,         # Mit neighbors zwischen 5 und 10 i.O
             min_dist = 0.5,
             random_state = 29)

## Plot mit base R

plot(umap$layout[,1], umap$layout[,2], pch = 20,
     main = "UMAP", xlab = "Umap1", ylab = "Umap2",
     col = colours[as.factor(data$type)])

legend("topright", legend = levels(as.factor(data$type)),
       pch = 20,
       col = colours,
       cex = 0.8,
       pt.cex = 1,
       bty = "n")

# 3. Clustering ----

## 3.1 K-means ----
'Schwierigkeiten mit hochdimensionale Daten (curse of dimensionality) füht zu
 wenig aussagekräftigen Ergebnissen'

### 3.1.1 optimal k ----

library(factoextra)

fviz_nbclust(dat, kmeans, method = "wss", k.max = 20) # Kein klarer ellbow
fviz_nbclust(dat, kmeans, method = "silhouette", k.max = 20) # 2 cluster

### 3.1.2 clusteren ----

set.seed(17)
kmean <- kmeans(dat, centers = 2)

plot(umap$layout[,1], umap$layout[,2], pch = 20,
     main = "UMAP", xlab = "Umap1", ylab = "Umap2",
     col = kmean$cluster)

table(kmeans = kmean$cluster, # 13 fehler
      actual = data$type)

'Ohne Dimensionsreduktion würden 2 Cluster empfohlen werden, wobei es so scheint als ob
normal und luminal_A und B ein Cluster bilden würden (Weniger aggressive krebszellen)'

set.seed(17)
kmean <- kmeans(dat, centers = 6)

plot(umap$layout[,1], umap$layout[,2], pch = 20,
     main = "UMAP", xlab = "Umap1", ylab = "Umap2",
     col = kmean$cluster)

table(kmeans = kmean$cluster, # 13 fehler
      actual = data$type)

'mit den "korrekten" 6 clustern funktioniert es gut. normal und cell_line sind 100% korrekt'

#### 3.1.1.1 K-means mit dimensionsreduktion ----

#### 3.1.1.2 optimal k ----

library(factoextra)

fviz_nbclust(umap$layout, kmeans, method = "wss", k.max = 20) # Kein klarer ellbow
fviz_nbclust(umap$layout, kmeans, method = "silhouette", k.max = 20) # 5 cluster

'Mit den Dimensionsreduktion werden schonmal 5 cluster statt nur 2 empfohlen'

#### 3.1.1.3 clusteren ----

set.seed(17)
kmean <- kmeans(umap$layout, centers = 6)

plot(umap$layout[,1], umap$layout[,2], pch = 20,
     main = "UMAP", xlab = "Umap1", ylab = "Umap2",
     col = kmean$cluster)

table(kmeans = kmean$cluster, # 25 fehler
      actual = data$type)


#### 3.1.1.4 Mittlere Silhouettenbreite ----
library(cluster)

kmeans <- kmeans(dat, centers=2)
plot(silhouette(x=kmeans$cluster, dist=dist(dat)))

'Bei k=2 -> Mittlere Silhouettenbreite = 0.14'

kmeans <- kmeans(dat, centers=6)
plot(silhouette(x=kmeans$cluster, dist=dist(dat)))

'Bei k=6 -> Mittlere Silhouettenbreite = 0.11'

'Beide weisen eine sehr ungeeignete Struktur auf (< 0.25). Da selbst mit
Kenntnissen über die wahre Anzahl Gruppen keine guten Strukturen geschaffen
werden können, scheint K-Means ungeeignet für diesen Datensatz zu sein.'

#### Erkenntnisse ----

'Wenn k-means auf die Rohdaten angewendet wird, werden nur 2 cluster empfohlen.
Bei den Daten mit Dimensionsreduktion durch UMAP sind es zumindest schon 5.
Wenn die optimale Anzahl von 6 Cluster verwendet wird, wird der Grossteil der
Daten ohne Dimensionsreduktion korrekt geclustert (13 fehler), bei den Daten
mit Dimensionsreduktion sind es jedoch 25 fehler, und Cell_line und normal werden sind ein cluster.
K-means clustering funktioniert also relativ gut mit den Rohdaten. Problem
ist jedoch, dass das empfohlene k bei den Rohdaten nicht nahe am optimum liegt.'

## 3.2 Pam ----

library(cluster)

### 3.2.1 optimal k ----

library(factoextra)

# Rohdaten
fviz_nbclust(dat, pam, method = "wss", k.max = 10) # Kein klarer ellbow
fviz_nbclust(dat, pam, method = "silhouette", k.max = 10) # 4 cluster

# Daten Dimensionsreduktion
fviz_nbclust(umap$layout, pam, method = "wss", k.max = 10) # 4 cluster
fviz_nbclust(umap$layout, pam, method = "silhouette", k.max = 10) # 2 oder 4 cluster

### 3.2.2 clusteren ----

set.seed(17)
pam <- pam(dat, k = 6) # Optimale anzahl cluster

plot(umap$layout, pch = 20,
     main = "UMAP", xlab = "Umap1", ylab = "Umap2",
     col = pam$clustering)

table(pam = pam$clustering, # 29 fehler
      actual = data$type)

'--------------------------------------------------------'

set.seed(17)
pam <- pam(umap$layout, k = 6) # Optimale anzahl cluster

plot(umap$layout, pch = 20,
     main = "UMAP", xlab = "Umap1", ylab = "Umap2",
     col = pam$clustering)

table(pam = pam$clustering, # 13 fehler
      actual = data$type)

### 3.2.3 Mittlere Silhouettenbreite ----

pam <- pam(dat, k=4)
plot(silhouette(x=pam$clustering, dist=dist(dat)))

'Bei k=2 -> Mittlere Silhouettenbreite = 0.13'

pam <- pam(dat, k=6)
plot(silhouette(x=pam$clustering, dist=dist(dat)))

'Bei k=6 -> Mittlere Silhouettenbreite = 0.10'

'Ähnliche mittlere Silhouettenbreiten wie bei K-Means. Dennoch immernoch
ungeeignete Struktur.'

### Erkenntnisse ----

'Pam sollte robuster als kmeans sein, da echte Beobachtungen als Zentren dienen.
Bei der bestimmung des optimalen k ist dies bei den Rohdaten der fall.
Beim clustering mit den Rohdaten schneided pam jedoch schlechter ab als kmeans.
Dafür ist das Ergebniss mit den Daten der Dimensionsreduktion besser (umgekehrt wie bie kmeans)'

## 3.3 Hierarchie-cluster ----

### 3.3.1 optimal k ----

# Rohdaten
hi.ward <- hclust(dist(dat), method = "ward.D2") # Most Variance
hi.single <- hclust(dist(dat), method = "single") # kürzeste Distanz
hi.complete <- hclust(dist(dat), method = "complete") # längste Distanz
hi.avg <- hclust(dist(dat), method = "average")

plot(hi.ward)
rect.hclust(hi.ward, h = 400) # mit dieser Methode hätte ich mich für 5 cluster entschieden.

# ungeeignet
plot(hi.single)
plot(hi.complete)
plot(hi.avg)

# Daten dimensionsreduktion
hi_d <- hclust(dist(umap$layout), method = "ward.D2")

plot(hi_d)
rect.hclust(hi_d, h = 20) # Auch hier wären es 5 cluster

### 3.3.2 clusteren ----

hi.c <- cutree(hi, k = 5)

plot(umap$layout, pch = 20,
     main = "UMAP", xlab = "Umap1", ylab = "Umap2",
     col = hi.c)

table(hi = hi.c, # 14 fehler
      actual = data$type)

'--------------------------------------------------------'

hi_d.c <- cutree(hi_d, k = 5)

plot(umap$layout, pch = 20,
     main = "UMAP", xlab = "Umap1", ylab = "Umap2",
     col = hi_d.c)

table(hi = hi_d.c, # 15 fehler
      actual = data$type)

### 3.3.3 Mittlere Silhouettenbreite ----
library(cluster)

plot(silhouette(x=hi.c, dist=dist(dat)))

### Erkenntnisse ----
'Man erhält mit Rohdaten und Daten der Dimensionsreduktion zumindest 5 cluster empfohlen.
Und auch beim Clustern schneiden beide Daten etwa gleich gut ab, wobei luminal_A und B
zusammen in einem cluster sind. Die mittlere Silhouettenbreite sagt, dass es sich um
eine ungeeignete Struktur handelt.'

## 3.4 DBscan ----

library(dbscan)

### 3.4.1 optimal params ----

# Rohdaten
kNNdistplot(dat, minPts = 5) # Nur 2 cluster nach dieser Methode
abline(h = 160, col = "red")

# Daten Dimensionsreduktion
kNNdistplot(umap$layout, minPts = 5) # 6 cluster nach dieser Methode
abline(h = 0.8, col = "red")

### 3.4.2 clusteren ----

db <- dbscan(dat, eps = 160, minPts = 5)

plot(umap$layout, pch = 20,
     main = "UMAP", xlab = "Umap1", ylab = "Umap2",
     col = db$cluster + 1)

table(dbscan = db$cluster + 1, 
      actual = data$type)

'----------------------------------------------------'

db_d <- dbscan(umap$layout, eps = 0.8, minPts = 5)

plot(umap$layout, pch = 20,
     main = "UMAP", xlab = "Umap1", ylab = "Umap2",
     col = db_d$cluster + 1)

table(dbscan = db_d$cluster + 1, 
      actual = data$type)

### Erkenntnisse ----
'Hierbei ist die Methode mit den Daten der Dimensionsreduktion genauer.
Diese gibt 6 cluster und mit den rohdaten wären es nur 2 cluster. 
Bei dbscan gibt es viel spielraum mit den parametern, also ist bestimmt noch
eine bessere Lösung möglich.

Interpretation: Durch UMAP wurden die Punkte besser separatiert, weshalb mit
den Dimensionsreduzierten Daten mehr Cluster entstanden sind.'

## 3.5 MClust ----

library(mclust)
library(factoextra)

### 3.5.1  optimal k and Model ----

# Daten Dimensionsreduktion
set.seed(17)
mc <- Mclust(umap$layout) # Rohdaten stürzt es ab, zu viele Variablen

mc$modelName
fviz_mclust(mc, "BIC") # Optimal 5 clusters

### 3.5.2 clusteren ----

plot(mc, what = "classification")

# 2. Visualisierung
colnames(mc$data) <- c("x", "y")
mc$data

fviz_mclust(mc, "classification",
            geom = "point")

table(pam = mc$classification, # 8 fehler aber 5 cluster 
      actual = data$type)

# Unsicherheit
fviz_mclust(mc, "uncertainty")

### Erkenntnisse ----
'Hierbei musste der Dimensionsreduzierte Datensatz verwendet werden. 
mclust erstellte dabei 5 cluster, wobei luminal A und B in einem Cluster sind.
Nur bei basal, HER und limunial B kam es zu "fehler".'

