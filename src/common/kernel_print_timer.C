

  const float mcpersec =
    ch->mcpara.steps_total * ch->size.n_rep / (e[4].Span () / 1000);
  printf ("mc kernel time\t\t\t%8.3f s\n", e[4].Span () / 1000);
  printf ("time per MC sweep per replica\t%8.3f * 1e-6 s\n", 1e6 / mcpersec);
  printf ("MC sweeps per second\t\t%8.3f\n", mcpersec);
  printf ("speedup over 1023.5\t\t%8.3f X\n", mcpersec / 1023.5);
  putchar ('\n');
  putchar ('\n');
  printf ("kernel:\n");
  printf ("kernel wall time\t\t%8.3f\n", e[10].Span ());
  printf ("alloc host\t\t\t%8.3f\n", e[0].Span ());
  printf ("H2Dcpy\t\t\t\t%8.3f\n", e[1].Span ());
  printf ("init 1\t\t\t\t%8.3f\n", e[2].Span ());
  printf ("init 2\t\t\t\t%8.3f\n", e[3].Span ());
  printf ("monte carlo\t\t\t%8.3f\n", e[4].Span ());
  printf ("\n");
