

  t[5].Stop ();


  // clean up
  delete[]mcpara;
  delete[]exchgpara;
  delete[]inputfiles->lig_files;
  delete[]inputfiles->lhm_files;
  delete[]inputfiles;


  for (int i = 0; i < conf_lig; ++i) {
    delete[]lig0s[i];
    delete[]psp0s[i];
    delete[]kde0s[i];
    delete[]mcs0s[i];
  }
  delete[]lig0s;
  delete[]psp0s;
  delete[]kde0s;
  delete[]mcs0s;
  delete[]prt0;
  delete[]enepara0;


  for (int i = 0; i < conf_lig; ++i) {
    delete[]ligs[i];
    delete[]psps[i];
    delete[]kdes[i];
    delete[]mcss[i];
    delete[]mcss_r[i];
    delete[]mcss_ell[i];
    delete[]mcss_csr[i];
  }
  delete[]ligs;
  delete[]psps;
  delete[]kdes;
  delete[]mcss;
  delete[]mcss_r;
  delete[]mcss_ell;
  delete[]mcss_csr;

  delete[]prt;
  delete[]enepara;

  delete[]complex;


  printf ("server:\n");
  printf ("server wall time\t\t%8.3f\n", t[5].Span ());
  printf ("init time\t\t\t%8.3f\n", t[6].Span ());
  printf ("load misc time\t\t\t%8.3f\n", t[7].Span ());
  printf ("load enepara time\t\t%8.3f\n", t[8].Span ());
  printf ("load prt time\t\t\t%8.3f\n", t[9].Span ());
  printf ("construct complex time\t\t%8.3f\n", t[10].Span ());
