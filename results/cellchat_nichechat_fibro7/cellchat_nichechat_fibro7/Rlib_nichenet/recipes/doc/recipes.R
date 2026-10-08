## -----------------------------------------------------------------------------
#| label: ex_setup
#| include: false
knitr::opts_chunk$set(
  message = FALSE,
  digits = 3,
  collapse = TRUE,
  comment = "#>",
  eval = requireNamespace("modeldata", quietly = TRUE) && requireNamespace("rsample", quietly = TRUE)
  )
options(digits = 3)


## -----------------------------------------------------------------------------
#| label: data
library(recipes)
library(rsample)
library(modeldata)

data("credit_data")

set.seed(55)
train_test_split <- initial_split(credit_data)

credit_train <- training(train_test_split)
credit_test <- testing(train_test_split)


## -----------------------------------------------------------------------------
#| label: missing
vapply(credit_train, function(x) mean(!is.na(x)), numeric(1))


## -----------------------------------------------------------------------------
#| label: first_rec
rec_obj <- recipe(Status ~ ., data = credit_train)
rec_obj


## -----------------------------------------------------------------------------
#| label: step_code
#| eval: false
# rec_obj <- step_{X}(rec_obj, arguments)    ## or
# rec_obj <- rec_obj |> step_{X}(arguments)


## -----------------------------------------------------------------------------
#| label: imp-steps
grep("impute_", ls("package:recipes"), value = TRUE)


## -----------------------------------------------------------------------------
#| label: dummy
imputed <- rec_obj |>
  step_impute_knn(all_predictors()) 
imputed


## -----------------------------------------------------------------------------
#| label: imputing
ind_vars <- imputed |>
  step_dummy(all_nominal_predictors()) 
ind_vars


## -----------------------------------------------------------------------------
#| label: center_scale
standardized <- ind_vars |>
  step_center(all_numeric_predictors()) |>
  step_scale(all_numeric_predictors()) 
standardized


## -----------------------------------------------------------------------------
#| label: trained
trained_rec <- prep(standardized, training = credit_train)
trained_rec


## -----------------------------------------------------------------------------
#| label: apply
train_data <- bake(trained_rec, new_data = credit_train)
test_data  <- bake(trained_rec, new_data = credit_test)


## -----------------------------------------------------------------------------
#| label: tibbles
class(test_data)
test_data
vapply(test_data, function(x) mean(!is.na(x)), numeric(1))


## -----------------------------------------------------------------------------
#| label: step_list
#| echo: false
grep("^step_", ls("package:recipes"), value = TRUE)


## -----------------------------------------------------------------------------
#| label: check
#| eval: false
# trained_rec <- trained_rec |>
#   check_missing(contains("Marital"))


## -----------------------------------------------------------------------------
#| label: check_list
#| echo: false
grep("^check_", ls("package:recipes"), value = TRUE)

