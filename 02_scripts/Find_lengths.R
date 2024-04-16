length_extremes <- ROD_v0.4_genome_stats %>%
  mutate(assembly_id = rassembly_id) %>%
  separate_rows(all_lengths, sep = ";") %>%
  mutate(all_lengths = as.numeric(all_lengths)) %>%
  group_by(row) %>%
  summarise(
    length_max = max(all_lengths),
    length_min = min(all_lengths)
  ) %>%
  ungroup()

length_extremes <- ROD_v0.4_genome_stats %>%
  separate_rows(all_lengths, sep = ";") %>%
  mutate(all_lengths = as.numeric(all_lengths)) %>%
  group_by(assembly_id) %>%
  summarise(
    length_min = min(all_lengths),
    length_max = max(all_lengths),
    .groups = 'drop'
  )



ROD_v0.4_genome_stats <- left_join(ROD_v0.4_genome_stats,length_extremes)

ROD_v0.4_genome_stats$length_mean
