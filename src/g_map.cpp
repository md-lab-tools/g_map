/*
  Copyright (C) 2015-2017
      Zidan Zhang, Jakub Krajniak

  This is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo++ is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "g_map.h"
#include <sstream>

#define NFILE asize(fnm)


/**
 *
 * @param fn
 * @param map_r
 */
void load_outer_rule(std::string fn, mapping_rule *map_r) {
  FILE *fp;
  int i, j;
  char temp[256];
  fp = fopen(fn.c_str(), "r");
  fgets(temp, 255, fp);
  strcpy(map_r->comment[0], temp);
  fgets(temp, 255, fp);
  sscanf(temp, "%d", &map_r->cg_res_num);
  for (i = 0; i < map_r->cg_res_num; i++) {
    fgets(temp, 255, fp);
    sscanf(temp, "%s%d", map_r->cg_name[i], &map_r->cg_num[i]);

    for (j = 0; j < map_r->cg_num[i]; j++) {
      fgets(temp, 255, fp);
      sscanf(temp, "%d", &map_r->residue_num[i][j]);
    }
  }
  fclose(fp);
}

/**
 *
 * @param fn
 * @param map_r
 */
void load_inner_rule(std::string fn, mapping_rule *map_r) {
  FILE *fp;
  int i, j, k;
  char temp[256];
  fp = fopen(fn.c_str(), "r");

  fgets(temp, 100, fp);

  strcpy(map_r->comment[1], temp);

  fgets(temp, 100, fp);

  int temp_num;
  sscanf(temp, "%d", &temp_num);

  if (temp_num != map_r->cg_res_num) {
    printf("Check your rule file! The residue number is not equal\n");
    exit(0);
  }

  for (i = 0; i < map_r->cg_res_num; i++) {
    fgets(temp, 255, fp);
    sscanf(temp, "%s%d", map_r->cg_name[i], &map_r->cg_particle_num[i]);

    for (j = 0; j < map_r->cg_particle_num[i]; j++) {
      fgets(temp, 255, fp);

      sscanf(temp, "%s%d", map_r->cg_particle_name[i][j], &map_r->atom_num_each_particle[i][j]);

      for (k = 0; k < map_r->atom_num_each_particle[i][j]; k++) {
        fgets(temp, 100, fp);
        sscanf(temp, "%s", map_r->atom_name_each_particle[i][j][k]);
      }
    }
  }
}

/**
 *
 * @param map_r maping rule
 * @param atoms atoms struct
 * @param x coordinate of original configuration
 * @param ih current head id
 * @param index_head head index
 * @param bGeo //whether using geometric center
 */
void map_traj(mapping_rule *map_r, t_atoms *atoms, rvec x[], int ih, int index_head[], bool bGeo) {
  int i;
  for (i = 0; i < ih; i++) {
    map_res(map_r, atoms, x, index_head[i], NULL, bGeo, FALSE);
  }
}

/**
 *
 * @param map_r
 * @param atoms
 * @param x
 * @param ih
 * @param index_head
 * @param bGeo
 */
void map_conf(mapping_rule *map_r, t_atoms *atoms, rvec x[], int ih, int index_head[], bool bGeo) {
  static t_symtab *symtab = NULL;
  snew(symtab, 1);
  open_symtab(symtab);

  int i;
  for (i = 0; i < ih; i++) {
    map_res(map_r, atoms, x, index_head[i], symtab, bGeo, TRUE);
  }
  sfree(symtab);

}

/*
 * these three procedures accomplish the mapping algorithm
the basic idea is to divide molecule configuration into
many residue, and in each residue we applied mapping method
to map atoms into coarse-grained particles. For the sake of
correctness, we asumme that all atoms in each residue must be
adjecent, and should be made a whole through rPBC procedue.
*/
/**
 *
 * @param map_r
 * @param atoms
 * @param x
 * @param head
 * @param symtab
 * @param bGeo
 * @param bMapName
 */
void map_res(mapping_rule *map_r, t_atoms *atoms, rvec x[], int head, t_symtab *symtab, bool bGeo, bool bMapName) {
  rvec xcom[CG_PART_TYPE];
  //paticle type number of specified residue
  int nres = resnm_to_resnr(map_r, *atoms->resinfo[atoms->atom[head].resind].name);

  for (int j = 0; j < map_r->cg_particle_num[nres]; j++) {
    clear_rvec(xcom[j]);
    real mtot = 0.0;
    for (int k = 0; k < map_r->atom_num_each_particle[nres][j]; k++) {
      int n = find_part_name_forward(atoms, head, map_r->atom_name_each_particle[nres][j][k]);
      int m = atoms->atom[n].m;
      if (bGeo)
        m = 1;
      for (int l = 0; l < DIM; l++)
        xcom[j][l] += m * x[n][l];
      mtot += m;
    }

    svmul(1 / mtot, xcom[j], xcom[j]);
    copy_rvec(xcom[j], x[head + j]);

    if (bMapName) {
      atoms->atomname[head + j] = put_symtab(symtab, map_r->cg_particle_name[nres][j]);
    }
  }
}

/**
 *
 * @param map_r
 * @param atoms
 * @param isize
 * @param index
 */
void proc_res_name(mapping_rule *map_r, t_atoms *atoms, int isize, int index[]) {
  static t_symtab *symtab = NULL;
  printf("Number of atoms: %d\n",isize);
  snew(symtab, 1);
  open_symtab(symtab);

  for (int i = 0; i < map_r->cg_res_num; i++) {
    for (int j = 0; j < map_r->cg_num[i]; j++) {
      atoms->resinfo[map_r->residue_num[i][j] - 1].
          name = put_symtab(symtab, map_r->cg_name[i]);
    }
  }
}

/**
 *
 * @param atoms
 * @param isize
 * @param index
 * @param ih
 * @param index_head
 */
void get_head(t_atoms *atoms, int isize, int index[], int *ih, int index_head[]) {
  index_head[0] = index[0];
  *ih = 1;  //number of the head index
  int current = 0;  //current position
  for (int i = 1; i < isize; i++) {
    int ia = index[i];
    if (atoms->atom[ia].resind != current) {
      *ih = *ih + 1;
      index_head[*ih - 1] = ia;
      current = atoms->atom[ia].resind;
    }
  }
}

/**
 *
 * @param map_r
 * @param res
 * @return
 */
int resnm_to_resnr(mapping_rule *map_r, char res[]) {
  for (int i = 0; i < map_r->cg_res_num; i++) {
    if (strcmp(map_r->cg_name[i], res) == 0)
      return i;
  }
  printf("resnm_to_resnr: no proper particle nr, check your rule file!");
  exit(0);
}

/**
 *
 * @param map_r
 * @param atoms
 * @param ih
 * @param index_head
 * @param io
 * @param index_output
 */
void make_output_index(mapping_rule *map_r, t_atoms *atoms, int ih, int index_head[], int *io, int index_output[]) {
  *io = 0;
  for (int i = 0; i < ih; i++) {
    int head = index_head[i];
    int t = resnm_to_resnr(map_r, *atoms->resinfo[atoms->atom[head].resind].name);
    int m = map_r->cg_particle_num[t];
    for (int j = 0; j < m; j++) {
      int n = head + j;
      index_output[*io] = n;
      *io = *io + 1;
    }
  }
}

/**
 *
 * @param map_r
 * @param atoms
 * @param from
 * @param atomnm
 * @return
 */
int find_part_name_forward(t_atoms *atoms, int from, char atomnm[]) {
  for (int i = from; i < atoms->nr; i++) {
    if (strcmp(*atoms->atomname[i], atomnm) == 0) {
      return i;
    }
  }

  printf("%s no such atom name, please check your rule file!\n", atomnm);
  exit(0);
}


/*! \brief
 * The main function for the analysis template.
 */
int main(int argc, char *argv[]) {
  static char *desc[] = {

  };

  static int n = 1;
  static bool bGeo = FALSE;
  int ngrps = 1;
  char title[STRLEN];
  t_trxframe fr;
  rvec *xtop;
  matrix box;
  t_trxstatus *status;
  int flags = TRX_READ_X;
  int isize_old, isize_new, isize_head, isize_output;
  char *grpname_old;
  int *index_old, *index_new, *index_head, *index_output;
  t_trxstatus *fp_confout;
  bool bIndex, bOutRule, bInnerRule;
  gmx_output_env_t *oenv;

  /*parameter input*/
  t_pargs pa[] = {
      {"-geo", FALSE, etBOOL, {&bGeo},
       "mapping with geometric center (default is center of mass)"}
  };
  /*file I/O part*/
  t_filenm fnm[] = {
      {efTPS, "-s", "topol", ffREAD},         /* this is the input topology */
      {efTRX, "-f", "traj", ffREAD},         /* this is the input trajectory */
      {efNDX, "-n", "index", ffREAD},    /* this is the grp to be coarse grained*/
      {efTRX, "-o", "trajout", ffWRITE},        /* this is the output topology */
      {efSTO, "-c", "confout", ffWRITE},        /* this is the output trajectory */
      {efDAT, "-or", "outer", ffREAD},
      {efDAT, "-ir", "inner", ffREAD},
  };


  if (!parse_common_args(&argc,
                         argv,
                         PCA_CAN_TIME | PCA_CAN_VIEW,
                         NFILE,
                         fnm,
                         asize(pa),
                         pa,
                         asize(desc),
                         (const char **) desc,
                         0,
                         NULL,
                         &oenv)) {
    return 0;
  }


  bIndex = opt2bSet("-n", NFILE, fnm);

  const char *infile_tps = ftp2fn_null(efTPS, NFILE, fnm);
  const char *infile_index = ftp2fn_null(efNDX, NFILE, fnm);

  int ePBC;
  t_topology top;
  bool bTop = (bool) read_tps_conf(infile_tps, &top, &ePBC, &xtop, NULL, box, FALSE);

  if (!bTop) {
    gmx_fatal(FARGS, "Could not read topology file from %s", ftp2fn(efTPS, NFILE, fnm));
  }

  get_index(&top.atoms, infile_index, 1, &isize_old, &index_old, &grpname_old);

  snew(index_head, isize_old);
  snew(index_output, isize_old);

  /* The first time we read data is a little special */
  read_first_frame(oenv, &status, ftp2fn(efTRX, NFILE, fnm), &fr, flags);

  //pre proc
  fp_confout = open_trx(opt2fn("-o", NFILE, fnm), "w");

  mapping_rule map_r;

  load_outer_rule(opt2fn("-or", NFILE, fnm), &map_r);    /*load the outer mapping rule*/
  load_inner_rule(opt2fn("-ir", NFILE, fnm), &map_r);    /*load the inner mapping rule*/

  proc_res_name(&map_r, &top.atoms, isize_old, index_old);
  get_head(&top.atoms, isize_old, index_old, &isize_head, index_head);
  make_output_index(&map_r, &top.atoms, isize_head, index_head, &isize_output, index_output);

  std::cout << "Loading rule file finished!" << std::endl;
  std::cout << "Processing trajectory" << std::endl;

  do {
    /* coordinates are available in the vector fr.x
    * you can find this and all other structures in
    * the types directory under the gromacs include dir.
    * Note how flags determines wheter to read x/v/f!
    */
    map_traj(&map_r, &top.atoms, fr.x, isize_head, index_head, bGeo);
    write_trxframe_indexed(fp_confout, &fr, isize_output, index_output, NULL);

  } while (read_next_frame(oenv, status, &fr));

  //close_trj(status);

  map_conf(&map_r, &top.atoms, xtop, isize_head, index_head, bGeo);

  std::stringstream ss;
  ss << "Input: " << infile_tps;

  write_sto_conf_indexed(ftp2fn(efSTO, NFILE, fnm),
                         ss.str().c_str(),
                         &top.atoms,
                         (rvec const *) xtop,
                         NULL,
                         ePBC,
                         box,
                         isize_output,
                         index_output);//need to be changed into new
  view_all(oenv, NFILE, fnm);
  return 0;
}
