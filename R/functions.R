filter_to_ppg <- function(wfo_classification) {
  # Neither "source" nor "majorGroup" is guaranteed to hit everything we want
  # (ferns and lycophytes), so combine the two

  sources_keep <- c(
    "Ferns and Allies TEN",
    "Michael Hassler",
    "Pteridophyte Phylogeny Group"
  ) |>
    paste(collapse = "|")

  groups_keep <- c()

  keep_by_source <- wfo_classification |>
    filter(str_detect(source, sources_keep)) |>
    pull(taxonID)

  keep_by_group <- wfo_classification |>
    filter(majorGroup %in% c("Lycopodiophyta", "Polypodiophyta")) |>
    pull(taxonID)

  keep_list <- unique(c(keep_by_source, keep_by_group))

  wfo_classification |>
    filter(taxonID %in% keep_list)
}

# Filter taxonomic data to ranks above species.
# Changes taxonRank to a factor in order of rank hierarchy
filter_to_above_species <- function(wfo_ppg) {
  ranks_keep <- c(
    "subgenus",
    "genus",
    "subtribe",
    "tribe",
    "subfamily",
    "family",
    "suborder",
    "order",
    "subclass",
    "class",
    "phylum"
  )

  ranks_keep <- factor(ranks_keep, levels = ranks_keep)

  wfo_ppg |>
    filter(taxonRank %in% ranks_keep) %>%
    mutate(taxonRank = factor(taxonRank, levels = ranks_keep))
}

# Helper function for join_higher_taxa()
add_parent_info <- function(df, parent_df, number) {
  parent_name_col <- sym(glue("parent_{number}_name"))
  parent_rank_col <- sym(glue("parent_{number}_rank"))
  highest_parent_name_col <- sym(glue("parent_{number - 1}_name"))

  df %>%
    left_join(
      unique(select(
        parent_df,
        scientificName,
        !!parent_rank_col := parentNameRank,
        !!parent_name_col := parentNameUsage
      )),
      by = setNames("scientificName", as.character(highest_parent_name_col))
    )
}

get_tax <- function(tax, level, highest = FALSE) {
  if (!highest) {
    str_match(tax, paste0("\\b", level, "\\b ([^ ]+) ")) |>
      magrittr::extract(, 2)
  } else {
    str_match(tax, paste0("\\b", level, "\\b ([^ ]+)")) |>
      magrittr::extract(, 2)
  }
}

tax_to_col_single <- function(tax) {
  subtribe <- get_tax(tax, "subtribe")
  tribe <- get_tax(tax, "tribe")
  subfamily <- get_tax(tax, "subfamily")
  family <- get_tax(tax, "family")
  suborder <- get_tax(tax, "suborder")
  order <- get_tax(tax, "order")
  subclass <- get_tax(tax, "subclass")
  class <- get_tax(tax, "class")
  phylum <- get_tax(tax, "phylum", highest = TRUE)
  tibble(
    subtribe = subtribe,
    tribe = tribe,
    subfamily = subfamily,
    family = family,
    suborder = suborder,
    order = order,
    subclass = subclass,
    class = class,
    phylum = phylum
  )
}

tax_to_col <- function(taxonomy) {
  map_df(taxonomy, tax_to_col_single)
}

# Join higher taxa (genus, tribe, subfamily, family, order) to PPG dataframe
widen_dwct <- function(wfo_pterido_higher) {
  
  # Prepare initial dataframe for joining parent taxa
  # Adds non-DWC 'parentNameRank' column
  wf_dwc_p <- wfo_pterido_higher |>
    filter(taxonomicStatus == "Accepted") |>
    dct_fill_col(
      fill_to = "parentNameUsage",
      fill_from = "scientificName",
      match_to = "taxonID",
      match_from = "parentNameUsageID"
    ) |>
    select(
      taxonID,
      scientificName, taxonRank,
      parentNameUsage
    ) |>
    select(
      taxonID, scientificName, taxonRank,
      parentNameUsage
    )

  wf_dwc_p <-
    wf_dwc_p %>%
    left_join(
      unique(select(wf_dwc_p, scientificName, parentNameRank = taxonRank)),
      join_by(parentNameUsage == scientificName),
      relationship = "many-to-one"
    ) |>
    assert(not_na, taxonID) |>
    assert(is_uniq, taxonID)

  wf_dwc_p %>%
    # progressively map on higher taxa
    rename(parent_1_name = parentNameUsage, parent_1_rank = parentNameRank) %>%
    relocate(parent_1_rank, .before = parent_1_name) |>
    add_parent_info(wf_dwc_p, 2) %>%
    add_parent_info(wf_dwc_p, 3) %>%
    add_parent_info(wf_dwc_p, 4) %>%
    add_parent_info(wf_dwc_p, 5) %>%
    add_parent_info(wf_dwc_p, 6) %>%
    add_parent_info(wf_dwc_p, 7) %>%
    add_parent_info(wf_dwc_p, 8) %>%
    # done when there are no more names to add
    # for species, should be done after adding 6 levels.
    # confirm that 7th doesn't add new information.
    verify(!all(is.na(parent_7_name))) %>%
    verify(all(is.na(parent_8_name))) %>%
    select(-contains("_8_")) |>
    filter(taxonRank == "genus") |>
    select(-taxonID) |>
    unite("taxonomy", matches("_name|_rank"), sep = " ") |>
    select(genus = scientificName, taxonomy) |>
    mutate(taxonomy_df = tax_to_col(taxonomy)) |>
    select(-taxonomy) %>%
    unnest(cols = taxonomy_df)
}

clean_ppgi <- function(ppgi_raw) {
  ppgi_raw |>
    filter(notes == "original PPGI 2016") |>
    select(genus, subfamily, family, suborder, order, class) |>
    mutate(genus = str_replace_all(genus, "Isoëtes", "Isoetes")) |>
    mutate(family = str_replace_all(family, "Isoëtaceae", "Isoetaceae")) |>
    mutate(subfamily = str_replace_all(subfamily, "Isoëtaceae", "Isoetaceae")) |>
    mutate(order = str_replace_all(order, "Isoëtales", "Isoetales"))
}

compare_taxonomy <- function(wfo_pterido_table, ppgi_table) {
  # For whatever reason, WFO data do not include suborder.
  # Nothing will match since they are all NA, so exclude this
  wfo_pterido_table |>
    select(-suborder) |>
    # FIXME: need to determine which is correct, "Thelypteroideae" or "Thelypteridoideae"
    mutate(
      subfamily = str_replace_all(
        subfamily, "Thelypteroideae", "Thelypteridoideae")) |>
    left_join(
      mutate(ppgi_table, matches = TRUE),
      na_matches = "na") |>
    full_join(ppgi_table, by = join_by(genus)) |>
    rename_with(~ str_replace_all(., "\\.x", "_wfo")) |>
    rename_with(~ str_replace_all(., "\\.y", "_ppgi")) |>
    rename_with(~ paste0(., "_wfo"), c(subtribe, tribe, subclass, phylum)) |>
    mutate(matches = replace_na(matches, FALSE))
}

filter_to_mismatches <- function(comparison_full) {
  comparison_full |>
    filter(!matches) |>
    select(-matches) |>
    arrange(
      phylum_wfo, class_wfo, subclass_wfo, order_wfo, family_wfo, subtribe_wfo,
      tribe_wfo, subtribe_wfo, genus)
}