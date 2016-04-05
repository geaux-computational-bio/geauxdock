#if IS_H5DUMP == 1
    // dump ligand record from CPU memory to disk
    char myoutputfile[MAXSTRINGLENG];
    sprintf(myoutputfile, "%s/%s_%04d.h5", outputdir, outputfile, s1 / steps_per_dump);
    DumpRecord (ch->record, ch->size.n_rep, myoutputfile);
#endif
