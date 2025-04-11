library(lubridate)
library(dplyr)
library(readr)

sqft_to_sqm <- 0.09290304

houses <- readr::read_csv("data-raw/kc_house_data.csv") |>
  mutate(
    date = year(date),
    # Convert square feet to square metres (metric system superiority)
    sq_m_living = sqft_living * sqft_to_sqm,
    sq_m_lot = sqft_lot * sqft_to_sqm,
    sq_m_above = sqft_above * sqft_to_sqm,
    sq_m_basement = sqft_basement * sqft_to_sqm,
    sq_m_living_15 = sqft_living15 * sqft_to_sqm,
    sq_m_lot_15 = sqft_lot15 * sqft_to_sqm,
    # Compute ages
    age = date - yr_built,
    yr_since_reno = date - pmax(yr_renovated, yr_built),
    # Convert to log values
    log_price = log(price),
    log_sq_m_living = log(sq_m_living),
    log_sq_m_lot = log(sq_m_lot),
    log_sq_m_living_15 = log(sq_m_living_15),
    log_sq_m_lot_15 = log(sq_m_lot_15)
  ) |>
  select(
    # Continuous variables
    age, yr_since_reno, log_sq_m_living, log_sq_m_lot, log_sq_m_living_15,
    log_sq_m_lot_15,
    # Discrete variables
    bedrooms, bathrooms, floors, waterfront, view, condition, grade,
    # Response variable
    log_price
  )

usethis::use_data(houses, overwrite = TRUE)
