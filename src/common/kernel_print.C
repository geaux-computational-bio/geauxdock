
  // read outputfiles and print for result verification
  const int iter_begin = 0;
  const int iter_end = 0;
  const int arg = 2;

//for (int r = 0; r < ph.complexsize->n_rep; ++r)
  for (int r = 0; r <= 0; ++r)
    PrintRecord (record, steps_per_dump, r, iter_begin, iter_end, arg);

  printf ("  0 0 0.6434 -0.0367 -0.2079 -0.1845 0.8521 -0.8880 0.0523 0.1744 -1.0000 0.7735 Ref Result\n"); 





#if 0
#if IS_H4DUMP == 1
  char myoutputfile[MAXSTRINGLENG];
  sprintf(myoutputfile, "%s/%s_%04d.h5", ph.mcpara->outputdir, ph.mcpara->outputfile, 0);
  Record *record;
  record = (Record *) malloc (ph.record_sz);
  ReadRecord (record, ph.complexsize->n_rep, myoutputfile);

  const int myreplica = 0;
  const int iter_begin = 0;
  const int iter_end = minimal_int (steps_per_dump, 1) - 1;
  const int arg = 2;

  PrintRecord (record, steps_per_dump, myreplica, iter_begin, iter_end, arg);
  //PrintMoveRecord (record, steps_per_dump, myreplica, iter_begin, iter_end, arg);

  free (record);
#endif
#endif

