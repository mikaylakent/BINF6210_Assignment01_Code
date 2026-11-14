##load all libraries ----
rm(list = ls())
library(tidyverse)
library(Biostrings)

##download T.nativa into script and create object ----
df_te_nat = read.csv("../Assign_01/Data/GCA_001447565.2.csv")

##seperate TEs from non-TEs for T. nativa ----
df_te_nat$matching_class |> unique()

ls_no_te_nat = c("Simple_repeat", "Low_complexity", "Satellite", "rRNA", "snRNA", "ARTEFACT", "Unknown")
ls_no_te_nat

df_te_nat$matching_class %notin% ls_no_te_nat |>
  sum()

df_te_nat |>
  filter(matching_class %notin% ls_no_te_nat) |> View()

df_te_nat$qry_id |> unique() |> length()

##seperate DNA and RNA TEs for T. nativa ----
df_dna_te_nat = df_te_nat |> 
  filter(str_starts(matching_class, "DNA/")) |>
  filter(perc_del > 0) |>
  group_by(matching_class) |>
  count() |>
  arrange(desc(n))
df_dna_te_nat

df_rna_te_nat = df_te_nat |> 
  filter(str_starts(matching_class, "PLE|PLE/|LTR/|LTR|LINE|LINE/|SINE/")) |>
  filter(perc_del > 0) |>
  group_by(matching_class) |>
  count() |>
  arrange(desc(n))

df_te_nat |> 
  filter(matching_class %notin% ls_no_te_nat) |>
  filter(qry_id %notin% df_dna_te_nat) |> View()

df_te_nat |> 
  filter(matching_class %notin% ls_no_te_nat) |> 
  filter(qry_id %notin% df_rna_te_nat) |> View()

##repeat the same seperating process for T. spiralis ----
df_trich_te = read.csv("../Assign_01/Data/trichinella_spiralis_te.csv")

df_trich_te$matching_class |> unique()

ls_no_te_trich = c("Simple_repeat", "simple_repeat", "Low_complexity", "Satellite", "rRNA", "snRNA", "ARTEFACT", "Unknown")
ls_no_te_trich

df_trich_te$matching_class %notin% ls_no_te_trich |>
  sum()

df_trich_te |>
  filter(matching_class %notin% ls_no_te_trich) |> View()

df_trich_te$qry_id |> unique() |> length()

df_dna_te_trich = df_trich_te |> 
  filter(str_starts(matching_class, "DNA/|RC")) |>
  filter(perc_del > 0) |>
  group_by(matching_class) |>
  count() |>
  arrange(desc(n))

df_rna_te_trich = df_trich_te |> 
  filter(str_starts(matching_class, "PLE|PLE/|LTR/|LTR|LINE|LINE/|SINE/")) |>
  filter(perc_del > 0) |>
  group_by(matching_class) |>
  count() |>
  arrange(desc(n))

##determine the timings of the TEs ----
