
# Daten ----

# Source:
# https://www.kaggle.com/datasets/fedesoriano/heart-failure-prediction?resource=download&select=heart.csv

dat <- read.csv("heart.csv", header = T, sep = ",")
str(dat)

# Es handelt sich um einen Datensatz zur Vorhersage von Herzerkrankungen.
# Der Datensatz wurde mit gemeinsamen Features aus 5 Datensätze zusammengestellt.
# Er beeinhaltet Kategorielle sowie numerische Werte.
# Beobachtungen sind jedoch nur 918.

'Attribute Information:
- Age: age of the patient [years]
- Sex: sex of the patient [M: Male, F: Female]
- ChestPainType: chest pain type [TA: Typical Angina, ATA: Atypical Angina, NAP: Non-Anginal Pain, ASY: Asymptomatic]
- RestingBP: resting blood pressure [mm Hg]
- Cholesterol: serum cholesterol [mm/dl]
- FastingBS: fasting blood sugar [1: if FastingBS > 120 mg/dl, 0: otherwise]
- RestingECG: resting electrocardiogram results [Normal: Normal, ST: having ST-T wave abnormality (T wave inversions and/or ST elevation or depression of > 0.05 mV), LVH: showing probable or definite left ventricular hypertrophy by Estes criteria]
- MaxHR: maximum heart rate achieved [Numeric value between 60 and 202]
- ExerciseAngina: exercise-induced angina [Y: Yes, N: No]
- Oldpeak: oldpeak = ST [Numeric value measured in depression]
- ST_Slope: the slope of the peak exercise ST segment [Up: upsloping, Flat: flat, Down: downsloping]
- HeartDisease: output class [1: heart disease, 0: Normal]'

# Cleaning ----

sum(is.na(dat)) # Keine Fehlerwerte
sum(duplicated(dat)) # Keine Dublikate

## Kategorielle Spalten umwandeln (Factor)
str(dat)
dat$Sex <- as.factor(dat$Sex)
dat$ChestPainType <- as.factor(dat$ChestPainType)
dat$RestingECG <- as.factor(dat$RestingECG)
dat$ExerciseAngina <-  as.factor(dat$ExerciseAngina)
dat$ST_Slope <- as.factor(dat$ST_Slope)
dat$HeartDisease <- as.factor(dat$HeartDisease)
dat$FastingBS <- as.factor(dat$FastingBS)

str(dat)

# Umkehren der Zielvariable, da sonst 0 für sensitivität und 1 für spec gilt
# Positive Klasse sollte Klasse von interesse sein, 1 = Herzerkrankung
dat$HeartDisease <- factor(dat$HeartDisease, levels = c("1", "0"))

str(dat)

# 1. Zielvariabel ----

# Binäre Zielvariable mit leichtem ungleichgewicht

table(dat$HeartDisease)
barplot(prop.table(table(dat$HeartDisease)))

'Denke kein Up-sampling nötig, da ungleichgewicht klein ist und die Klasse welche 
uns interessiert mehr Beobachtungen hat.'

# 2. Aufteilung ----

library(tidymodels)

split <- initial_split(dat, prop = 0.8, strata = HeartDisease) 
train_set <- training(split)
test_set <- testing(split) # Nur 184 Beobachtungen in test, vlt doch prop = 0.75?

# 3. Vorverarbeitung ----

str(dat)
boxplot(dat[c(1, 4, 5, 8, 10)], las = 2) # Skalierung notwendig

# Rezept
recipe <- recipe(HeartDisease ~., data = train_set) |> 
  step_normalize(all_numeric_predictors())

# 4. Setting ----

# Nicht viele Features, Fluch der Dimensionalität also kein Problem (KNN und Diskriminanzanalysen)
# 6 der 11 Features sind kategoriell, weshalb Diskriminanzanalysen und KNN weniger geeignet sind (Evtl. hot-encoding bei KNN).
# Evtl. naive bayes mit Gauss, da es sonst vermutlich zu informationsverlust kommt, wenn man stetige Daten diskretisiert und in Gruppen einteilt. 

# Random Forest
# Logistische Regression
# naive Bayes (mit Gauss wäre vermutlich optimal aber hatten wir nicht)

# Kreuzvalidierung
cv_fold <- vfold_cv(train_set, v = 10, repeats = 3, strata = HeartDisease)

## Naive-Bayes ----
library(discrim)

# Model
naiv <- naive_Bayes(smoothness = tune()) |> 
  set_engine("naivebayes") |> 
  set_mode("classification")

# Workflow
naiv_wf <- workflow() |> 
  add_model(naiv) |> 
  add_recipe(recipe)

## Logistische Regression ----
library(MASS)

# Model
log <- logistic_reg() |> 
  set_engine("glm") |> 
  set_mode("classification")

# Workflow
log_wf <- workflow() |> 
  add_model(log) |> 
  add_recipe(recipe)

## Random Forest ----
library(ranger)

# Model
rf <- rand_forest(mtry = tune(), trees = 500) |> 
  set_engine("ranger") |> 
  set_mode("classification")

# Workflow
rf_wf <- workflow() |> 
  add_model(rf) |> 
  add_recipe(recipe)

### Tune() ----

# Tune naive-Bayes Smoothness
naiv_tune <- tune_grid(naiv_wf,
                       resamples = cv_fold,
                       grid = grid_regular(smoothness(range = c(0.1, 3)), levels = 10),
                       metrics = metric_set(mn_log_loss),
                       control = control_grid(save_pred = T))

naiv_best <- select_best(naiv_tune, metric = "mn_log_loss")
naiv_final <- finalize_workflow(naiv_wf, naiv_best)

# Tune random forest mtry
rf_tune <- tune_grid(rf_wf,
                     resamples = cv_fold,
                     grid = grid_regular(mtry(range = c(1, 11)), levels = 11),
                     metrics = metric_set(mn_log_loss),
                     control = control_grid(save_pred = T))

rf_best <- select_best(rf_tune, metric = "mn_log_loss")
rf_final <- finalize_workflow(rf_wf, rf_best)

# 5. Vergleich ----

mymetrics <- metric_set(accuracy, mn_log_loss, roc_auc, precision, sens, spec)

# Resample Naive-Bayes
naiv_valid <- fit_resamples(naiv_final,
                          resamples = cv_fold,
                          metrics = mymetrics,
                          control = control_resamples(save_pred = T))

# Resample logistical Regression
log_valid <- fit_resamples(log_wf,
                        resamples = cv_fold,
                        metrics = mymetrics,
                        control = control_resamples(save_pred = T))

# Resample Random Forest
rf_valid <- fit_resamples(rf_final,
                        resamples = cv_fold,
                        metrics = mymetrics,
                        control = control_resamples(save_pred = T))

# Best Performance
collect_metrics(naiv_valid)
collect_metrics(log_valid)
collect_metrics(rf_valid) # Winner winner chicken dinner

# Confusion Matrix
collect_predictions(naiv_valid) |> 
  conf_mat(truth = HeartDisease, estimate = .pred_class)

collect_predictions(log_valid) |> 
  conf_mat(truth = HeartDisease, estimate = .pred_class)

collect_predictions(rf_valid) |> 
  conf_mat(truth = HeartDisease, estimate = .pred_class)

'Random Forest ist unsere wahl, es hat den tiefsten log-loss und den höchsten AUC, 
auch wenn nur knapp.
Zudem ist es besser darin Herzkrankheiten frühzeitig zu erkennen. Auch wenn es
im gegenzug minimal mehr risikoarme Patienten als Risikopatient einstuft.
'

# 6. Test ----

rf_fit <- last_fit(rf_final, split = split, metrics = mymetrics)

collect_predictions(rf_fit) |> 
  conf_mat(truth = HeartDisease, estimate = .pred_class)

collect_metrics(rf_fit)

'Unser Modell ist besser bei der Vorhersage der positiven Klasse als bei der negativen.
Die ist zum Glück was wir möchten, da es wichtig ist, alle risikopatienten frühzeitig zu erkennen.
Wir vermuten, dass dies mit dem Klassenungleichgewicht zugunsten der positiven Klasse zu tun hat.'

# Visuell

pred <- collect_predictions(rf_fit)
pred

# ROC
library(pROC)

p_roc <- roc(test_set$HeartDisease, pred$.pred_1, levels = c("0", "1"))

plot(p_roc, las = 2)
p_roc$auc

# Kalibrierungskurve
library(predtools)

str(pred)
pred2 <- data.frame(pred,
                    obs = ifelse(test_set$HeartDisease == "1", 1, 0)) 
str(pred2)

calibration_plot(pred2,
                 obs = "obs",
                 pred = ".pred_1",
                 nTiles = 10)

'In bezug auf die Kalibrierung, sehen wir, dass das Modell bei Wahrscheinlichkeiten unter
50% die Wahrscheinlichkeit unterschätzt, es liegt also öfters richtig als das Modell denkt.
Jedoch ist es sich dann bei höheren Wahrscheinlichkeiten etwas zu sicher (minimal)'

# Verbesserung ----

'Obwohl unser modell mit der positiven Klasse besser umgehen kann (vermutlich Klassenungleichgewicht)
könnte man diese Eigenschaft noch etwas verbessern.

z.B könnte man die positive Klasse mehr gewichten (GPT):
set_engine("ranger", class.weights = c("0" = 1, "1" = 3))

oder man könnte beim fertigen Modell wie folgt forgehen:
'
rf_fit_optim <- fit(rf_final, data = train_set)

prediction <- predict(rf_fit_optim, new_data = test_set, type = "prob")
prediction

pred_class <- ifelse(pred$.pred_1 > 0.3, "1", "0")

# Result
m <- table(prediction = pred_class, truth = test_set$HeartDisease)
m

m[2,1] / sum(m[1,1] + m[2,1]) # Sens um 3% erhöht
