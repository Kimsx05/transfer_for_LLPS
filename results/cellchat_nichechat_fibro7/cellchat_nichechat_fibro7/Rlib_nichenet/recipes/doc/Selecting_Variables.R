## -----------------------------------------------------------------------------
#| label: ex_setup
#| include: false
knitr::opts_chunk$set(
  message = FALSE,
  digits = 3,
  collapse = TRUE,
  comment = "#>",
  eval = requireNamespace("modeldata", quietly = TRUE)
  )
options(digits = 3)


## -----------------------------------------------------------------------------
#| label: penguins
library(recipes)
library(modeldata)

data("penguins")
str(penguins)

rec <- recipe(body_mass_g ~ ., data = penguins)
rec


## -----------------------------------------------------------------------------
#| label: var_info_orig
summary(rec, original = TRUE)


## -----------------------------------------------------------------------------
#| label: var_info_orig_type
summary(rec, original = TRUE)$type


## -----------------------------------------------------------------------------
#| label: dummy_1
dummied <- rec |> step_normalize(all_numeric())


## -----------------------------------------------------------------------------
#| label: dummy_2
dummied <- rec |> step_normalize(bill_length_mm, bill_depth_mm, 
                                  flipper_length_mm) # or
dummied <- rec |> step_normalize(all_numeric(), - body_mass_g) # or
dummied <- rec |> step_normalize(all_numeric_predictors()) # recommended


## -----------------------------------------------------------------------------
rec |>
  step_dummy(sex) |>
  prep() |>
  juice()


## -----------------------------------------------------------------------------
#| label: dummy_3
dummied <- prep(dummied, training = penguins)
with_dummy <- bake(dummied, new_data = penguins)
with_dummy

