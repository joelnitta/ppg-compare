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
  # Load results of PPG voting
  tar_target(
    ppg_vote_results,
    fetch_ppg_vote_results(),
    cue = tar_cue("always")
  ),
  # Compare classifications ----
  comparison_full = compare_taxonomy(wfo_pterido_table, ppgi_table),
  comparison_mismatches = filter_to_mismatches(comparison_full),
  comparison_mismatches_with_issues = join_issues(
    comparison_mismatches, ppg_vote_results),
  tar_file(
    comparison_xlsx,
    write_xlsx_tar(
      comparison_mismatches_with_issues,
      "data/comparison_raw.xlsx"),
  )
)
