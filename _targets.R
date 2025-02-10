source("R/packages.R")
source("R/functions.R")

tar_plan(
  # Load pteridophyte taxonomy data from WFO ----
  # - full classification (all plants)
  tar_file_read(
    wfo_classification_full,
    "data/classification.csv",
    read_tsv(!!.x)
  ),
  # - filter to only pteridophytes
  wfo_pterido = filter_to_ppg(wfo_classification_full),
  # - filter to genus and higher
  wfo_pterido_higher = filter_to_above_species(wfo_pterido),
  # - convert to table (one column per taxonomic rank)
  wfo_pterido_table = widen_dwct(wfo_pterido_higher),
  # Load pteridophyte taxonomy data from PPG I ----
  tar_file_read(
    ppgi_raw,
    "data/ppgi.csv",
    read_csv(!!.x)
  ),
  ppgi_table = clean_ppgi(ppgi_raw),
  # Compare classifications ----
  comparison_full = compare_taxonomy(wfo_pterido_table, ppgi_table),
  comparison_mismatches = filter_to_mismatches(comparison_full)
)
