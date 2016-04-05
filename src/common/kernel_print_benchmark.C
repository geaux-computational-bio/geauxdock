printf ("Benchmark,\t\t");
printf ("%s,", ch->lig_file);
printf ("%s,", ch->prt_file);
printf ("%d,", ch->mcpara.steps_total);
printf ("%d,", ch->mcpara.steps_per_dump);
printf ("%d,", ch->size.n_lig);
printf ("%d,", ch->size.n_prt);
printf ("%d,", ch->size.n_tmp);
printf ("%d,", ch->size.n_rep);
printf ("%d,", ch->size.lig_natom);
printf ("%d,", ch->size.prt_npoint);
printf ("%d,", ch->size.kde_npoint);
printf ("%d,", ch->size.mcs_nrow);
printf ("%.3f", e[4].Span ()); // milli-seconds
printf ("\n");





#if IS_PAPI == 1
printf ("Benchmark_papi,\t\t");
printf ("%s,", ch->lig_file);
printf ("%s,", ch->prt_file);
printf ("%d,", ch->mcpara.steps_total);
printf ("%d,", ch->mcpara.steps_per_dump);
printf ("%d,", ch->size.n_lig);
printf ("%d,", ch->size.n_prt);
printf ("%d,", ch->size.n_tmp);
printf ("%d,", ch->size.n_rep);
printf ("%d,", ch->size.lig_natom);
printf ("%d,", ch->size.prt_npoint);
printf ("%d,", ch->size.kde_npoint);
printf ("%d,", ch->size.mcs_nrow);
printf ("%.3f,", e[4].Span () / 1000.0f); // seconds
for (int i = 0; i < m0.papi_event_n; ++i)
  printf ("%.3f,", (float) m0.papi_event_val[i] / 1000000.0f);
printf ("\n");
#endif


