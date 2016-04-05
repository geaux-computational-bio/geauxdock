
  yeah::Timer t[16];
  t[5].Start ();

  //Banner ();

  t[6].Start ();
  McPara *mcpara = new McPara;
  ExchgPara *exchgpara = new ExchgPara;
  InputFiles *inputfiles = new InputFiles[1];
  ParseArguments (argc, argv, mcpara, exchgpara, inputfiles);


  // number of complexes
  const int conf_lig = inputfiles->lig_files_num;
  const int conf_prt = 1;
  const int ncomplex = conf_lig * conf_prt;
  Complex *complex = new Complex[ncomplex];

  printf ("task pool size = %d\n", ncomplex);



  Ligand0 **lig0s = new Ligand0*[conf_lig];
  Psp0 **psp0s = new Psp0*[conf_lig];
  Kde0 **kde0s = new Kde0*[conf_lig];
  Mcs0 **mcs0s = new Mcs0*[conf_lig];
  Ligand **ligs = new Ligand*[conf_lig];
  Psp **psps = new Psp*[conf_lig];
  Kde **kdes = new Kde*[conf_lig];
  Mcs **mcss = new Mcs*[conf_lig];
  Mcs_R **mcss_r = new Mcs_R*[conf_lig];
  Mcs_ELL **mcss_ell = new Mcs_ELL*[conf_lig];
  Mcs_CSR **mcss_csr = new Mcs_CSR*[conf_lig];

  for (int i = 0; i < conf_lig; ++i) {
    lig0s[i] = new Ligand0[MAX_CONF_LIG];
    psp0s[i] = new Psp0;
    kde0s[i] = new Kde0;
    mcs0s[i] = new Mcs0[MAX_MCS_ROW];
    ligs[i] = new Ligand[MAX_CONF_LIG];
    psps[i] = new Psp;
    kdes[i] = new Kde;
    mcss[i] = new Mcs[MAX_MCS_ROW];
    mcss_r[i] = new Mcs_R;
    mcss_ell[i] = new Mcs_ELL;
    mcss_csr[i] = new Mcs_CSR;
  }
  t[6].Stop ();


  t[7].Start ();
//#pragma omp parallel for
  for (int i = 0; i < conf_lig; ++i) {
    //printf ("Loading & pre-computing lig,kde,mcs %3d\n", i);
    //std::cout << "lhm_file: " << inputfiles->lig_files[i].path << std::endl;

    loadLigand (&inputfiles->lig_files[i], lig0s[i]);
    loadLHM (&inputfiles->lhm_files[i], psp0s[i], kde0s[i], mcs0s[i]);

    OptimizeKde (kde0s[i], kdes[i]);
    OptimizeLigand (lig0s[i], kdes[i], ligs[i], inputfiles->lig_files[i].conf_total);
    OptimizePsp (psp0s[i], psps[i]);
    OptimizeMcs (mcs0s[i], mcss[i], mcss_r[i], mcss_ell[i], mcss_csr[i], inputfiles->lhm_files[i].mcs_nrow);
  }
  t[7].Stop ();

  
  //printf ("Loading enepara\n");

  t[8].Start ();
  EnePara0 *enepara0 = new EnePara0;
  EnePara *enepara = new EnePara;
  loadEnePara (&inputfiles->enepara_file, enepara0);
  loadWeight(&inputfiles->weight_file, enepara0);
  loadNorPara(&inputfiles->norpara_file, enepara0);
  OptimizeEnepara (enepara0, enepara);
  t[8].Stop ();


  //printf ("Loading protein\n");

  t[9].Start ();
  Protein0 *prt0 = new Protein0[MAX_CONF_PRT];
  Protein *prt = new Protein[MAX_CONF_PRT]; // complexsize.n_prt
  loadProtein (&inputfiles->prt_file, prt0);
  t[9].Stop ();




  t[10].Start ();

// sould dynamically construct
//#pragma omp parallel for
  for (int i = 0; i < ncomplex; ++i) {

    //printf ("Comstructing complex %3d\n", i);

    const float cx = lig0s[i][0].pocket_center[0];
    const float cy = lig0s[i][0].pocket_center[1];
    const float cz = lig0s[i][0].pocket_center[2];
    SetPocketCenter (prt, cx, cy, cz, inputfiles->prt_file.conf_total);
    OptimizeProtein (prt0, prt, enepara0, inputfiles->prt_file.conf_total); // should optimize on the fly


    ComplexSize sz;
    sz.n_prt = inputfiles->prt_file.conf_total;	// number of protein conf
    sz.n_lig = inputfiles->lig_files[i].conf_total;	// number of ligand conf
    sz.n_tmp = exchgpara->num_temp;	// number of temperature
    sz.n_rep = sz.n_lig * sz.n_prt * sz.n_tmp;
    sz.lig_natom = inputfiles->lig_files[i].lig_natom;
    sz.prt_npoint = inputfiles->prt_file.prt_npoint;
    sz.kde_npoint = kde0s[i]->kde_npoint;
    sz.mcs_nrow = inputfiles->lhm_files[i].mcs_nrow;	// number of MCS positions


    if (sz.n_prt > MAX_CONF_PRT) {
      std::cout << inputfiles->prt_file.id << "conformation number exceeds the limit ";
      std::cout << sz.n_prt << " > " << MAX_CONF_PRT << std::endl;
      exit (1);
    }
    if (sz.n_lig > MAX_CONF_LIG) {
      std::cout << inputfiles->lig_files[i].id << "conformation number exceeds the limit ";
      std::cout << sz.n_lig << " > " << MAX_CONF_LIG << std::endl;
      exit (1);
    }
    if (sz.n_tmp > MAX_TMP) {
      std::cout << "temperature number exceeds the limit ";
      std::cout << sz.n_tmp << " > " << MAX_TMP << std::endl;
      exit (1);
    }
    if (sz.n_rep > MAX_REP) {
      std::cout << "ensemble number execeeds the limit ";
      std::cout << sz.n_rep << " > " << MAX_REP << std::endl;
      exit (1);
    }
    if (sz.lig_natom > MAXLIG) {
      std::cout << inputfiles->lig_files[i].id << "point number execeeds the limit ";
      std::cout << sz.lig_natom << " > " << MAXLIG << std::endl;
      exit (1);
    }
    if (sz.prt_npoint > MAXPRO) {
      std::cout << inputfiles->prt_file.id << "point number execeeds the limit ";
      std::cout << sz.prt_npoint << " > " << MAXPRO << std::endl;
      exit (1);
    }
    if (sz.kde_npoint > MAXKDE) {
      std::cout << "KDE(" << i << ")point number execeeds the limit ";
      std::cout << sz.kde_npoint << " > " << MAXKDE << std::endl;
      exit (1);
    }
    if (sz.mcs_nrow > MAX_MCS_ROW) {
      std::cout << "MCS(" << i << ")point number execeeds the limit ";
      std::cout << sz.mcs_nrow << " > " << MAX_MCS_ROW << std::endl;
      exit (1);
    }




    Temp *temp = new Temp[MAX_TMP]; // complexsize.n_tmp
    ReplicaMC *replica = new ReplicaMC[MAX_REP]; // complexsize.n_rep
    SetTemperature (temp, exchgpara);
    SetReplica (replica, sz);


    Complex *cc = &complex[i];

    cc->size = sz;
    cc->mcpara = *mcpara;
    cc->psp = *psps[i];
    cc->kde = *kdes[i];
    cc->enepara = *enepara;
    for (int j = 0; j < sz.n_lig; ++j) cc->lig[j] = ligs[i][j];
    for (int j = 0; j < sz.n_prt; ++j) cc->prt[j] = prt[j];
    for (int j = 0; j < sz.mcs_nrow; ++j)   cc->mcs[j] = mcss[i][j];
    cc->mcs_r =  mcss_r[i][0];
    cc->mcs_ell =  mcss_ell[i][0];
    cc->mcs_csr =  mcss_csr[i][0];
    for (int j = 0; j < sz.n_tmp; ++j) cc->temp[j] = temp[j];
    for (int j = 0; j < sz.n_rep; ++j) cc->replica[j] = replica[j];




    //std::strcpy (cc->lig_file, inputfiles->lig_files[i].path.c_str ());
    //std::strcpy (cc->prt_file, inputfiles->prt_file.path.c_str ());
    //std::strcpy (cc->lig_file, inputfiles->lig_files[i].id.c_str ());
    //std::strcpy (cc->prt_file, inputfiles->prt_file.id.c_str ());

    // get the file name from a path
    std::string aa;
    size_t last_slash_idx;

    aa = inputfiles->lig_files[i].path;
    last_slash_idx = aa.find_last_of ("/");
    if (std::string::npos != last_slash_idx) aa.erase (0, last_slash_idx + 1);
    std::strcpy (cc->lig_file, aa.c_str ());

    aa = inputfiles->prt_file.path;
    last_slash_idx = aa.find_last_of ("/");
    if (std::string::npos != last_slash_idx) aa.erase (0, last_slash_idx + 1);
    std::strcpy (cc->prt_file, aa.c_str ());



    std::strcpy (cc->lhm_file, inputfiles->lhm_files[i].path.c_str ());
    std::strcpy (cc->enepara_file, inputfiles->enepara_file.path.c_str ());
    std::strcpy (cc->weight_file, inputfiles->weight_file.path.c_str ());
    cc->rep_begin = 0;
    cc->rep_end = sz.n_rep - 1;
    cc->signal = i;

    delete[]temp;
    delete[]replica;
  }
  t[10].Stop ();





