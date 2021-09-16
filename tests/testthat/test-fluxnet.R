

test_that("get_fluxnet_metadata", {
  site <- "FI-Hyy"
  md <- get_fluxnet_metadata(site)
  num_names <- c(
    "Latitude","Longitude","Elevation..m."
    , "Mean.Annual.Temperature..degrees.C.", "Mean.Annual.Precipitation..mm."
  )
  expected_names <- c(
    num_names
    , "Site.ID","Site.Name","Network", "IGBP"
    , "Data.Products", "Data.Availability", "Data.Downloads.to.Date", "Data.DOIs"
    , "contact.name", "contact.email")
  expect_true(all(expected_names %in% names(md)))
  expect_equal(md$Site.ID, site)
  expect_true(all(sapply(md[num_names],is.numeric)))
})
