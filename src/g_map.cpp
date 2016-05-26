/*
  Copyright (C) 2016
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

#include <string>
#include <vector>

#include "gromacs/trajectoryanalysis.h"
#include <gromacs/options/ioptionscontainer.h>
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/pbcutil/pbc.h"
#include <cstring>
#include <iostream>
#include <gromacs/fileio/trxio.h>

using namespace gmx;

#define CG_TYPE 10
#define CG_PART_TYPE 20
#define AA_RES_NUM 9999
#define AA_PART_NUM 999

typedef struct {
  char comment[2][255];                            //comment [0]:outer [1]:inner
  int cg_res_num;                                //total number of CG residues
  char cg_name[CG_TYPE][5];                        //CG residue names

  int cg_num[CG_TYPE];                            //corresponding AA residue numbers
  int residue_num[CG_TYPE][AA_RES_NUM];                    //AA residue numbers in each CG residue

  int cg_particle_num[CG_TYPE];                        //CG particle numbers for each CG residue
  char cg_particle_name[CG_TYPE][CG_PART_TYPE][5];            //CG particle names
  int atom_num_each_particle[CG_TYPE][CG_PART_TYPE];            //atom number in each CG particle
  char atom_name_each_particle[CG_TYPE][CG_PART_TYPE][AA_PART_NUM][5];    //atom name in each CG particle

} mapping_rule;

/*! \brief
 * Template class to serve as a basis for user analysis tools.
 */
class MappingAT2CG: public TrajectoryAnalysisModule {
 public:
  MappingAT2CG();

  virtual void initOptions(IOptionsContainer *options,
                           TrajectoryAnalysisSettings *settings);
  virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                            const TopologyInformation &top);

  virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                            TrajectoryAnalysisModuleData *pdata);

  virtual void finishAnalysis(int nframes);
  virtual void writeOutput();

 private:
  class ModuleData;

  void load_outer_rule(std::string fn, mapping_rule *map_r);
  void load_inner_rule(std::string fn, mapping_rule *map_r);

  t_trxstatus *f_confout;

  std::string fn_topol_;
  double cutoff_;
  Selection refsel_;
  SelectionList sel_;

  std::string fn_outer_;
  std::string fn_inner_;
  std::string fn_confout_;

  int isize_old, isize_new, isize_head, isize_output;
  int *index_old, *index_head, *index_output;;
  bool bGeo;

  mapping_rule map_r;

  int resnm_to_resnr(mapping_rule *map_r, char *res);
  void map_conf(mapping_rule *map_r, t_atoms *atoms, rvec x[], int ih, int index_head[], bool bGeo);
  void map_res(mapping_rule *map_r, t_atoms *atoms, rvec x[], int head, t_symtab *symtab, bool bGeo, bool bMapName);
  void map_traj(mapping_rule *map_r, t_atoms *atoms, rvec[], int ih, int index_head[], bool bGeo);
  void proc_res_name(mapping_rule *map_r, t_atoms *atoms, int isize, int index[]);
  void get_head(t_atoms *atoms, int isize, int index[], int *ih, int index_head[]);
  void make_output_index(mapping_rule *map_r, t_atoms *atoms,
                         int ih, int index_head[],
                         int *io, int index_output[]);
  int find_part_name_forward(mapping_rule *map_r, t_atoms *atoms, int from, char atomnm[]);

};


MappingAT2CG::MappingAT2CG()
    : cutoff_(0.0) {
}


void
MappingAT2CG::initOptions(IOptionsContainer *options,
                          TrajectoryAnalysisSettings *settings) {
  static const char *const desc[] = {"Mapping AT to CG"};
  settings->setHelpText(desc);
  options->addOption(FileNameOption("o")
                         .filetype(eftTopology).outputFile().required(true)
                         .store(&fn_topol_).defaultBasename("topol")
                         .description("Output CG topology"));
  options->addOption(FileNameOption("or")
                         .filetype(eftGenericData)
                         .inputFile().required(true)
                         .store(&fn_outer_)
                         .description("Outer file"));
  options->addOption(FileNameOption("ir")
                         .filetype(eftGenericData)
                         .inputFile().required(true)
                         .store(&fn_inner_)
                         .description("Outer file"));
  options->addOption(FileNameOption("c")
                         .filetype(eftTrajectory)
                         .outputFile().required(true)
                         .store(&fn_confout_)
                         .defaultBasename("trjout")
                         .description("Output trajectory"));
  options->addOption(BooleanOption("geo")
                         .store(&bGeo)
                         .defaultValue(false)
                         .description("mapping with geometric center (default is center of mass)"));

  settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
  settings->setFlag(TrajectoryAnalysisSettings::efUseTopX);
}


void
MappingAT2CG::initAnalysis(const TrajectoryAnalysisSettings &settings,
                           const TopologyInformation &top) {
  load_outer_rule(fn_outer_, &map_r);
  load_inner_rule(fn_inner_, &map_r);

  proc_res_name(&map_r, &top.topology()->atoms, isize_old, index_old);
  get_head(&top.topology()->atoms, isize_old, index_old, &isize_head, index_head);

  std::cout << "Loading rule files finished" << std::endl;

  f_confout = open_trx(fn_confout_.c_str(), "w");
}

void MappingAT2CG::analyzeFrame(int frnr, const t_trxframe &fr,
                                t_pbc *pbc, TrajectoryAnalysisModuleData *pdata) {
  rm_gropbc(fr.atoms, fr.x, fr.box);
  map_traj(&map_r, fr.atoms, fr.x, isize_head, index_head, bGeo);
  //put_atoms_in_box(int ePBC, const matrix box, int natoms, rvec x[])
  put_atoms_in_box(pbc->ePBC, fr.box, fr.atoms->nr, fr.x);
  //t_trxframe frame;

  write_trxframe(f_confout, (t_trxframe *) &fr, NULL);
}


void MappingAT2CG::finishAnalysis(int /*nframes*/) {
  std::cout << "finished!" << std::endl;
}


void MappingAT2CG::writeOutput() {

}


/** Helper methods. */
void MappingAT2CG::load_outer_rule(std::string fn, mapping_rule *map_r) {
  FILE *fp;
  int i, j;
  char temp[256];
  fp = fopen(fn.c_str(), "r");
  fgets(temp, 255, fp);
  printf("\n-loading outer_rule-: %s \n", temp);
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

void MappingAT2CG::load_inner_rule(std::string fn, mapping_rule *map_r) {
  FILE *fp;
  int i, j, k;
  char temp[256];
  fp = fopen(fn.c_str(), "r");

  fgets(temp, 100, fp);

  printf("\n-loading inner_rule-: %s \n", temp);

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

void MappingAT2CG::map_traj(mapping_rule *map_r,            //mapping rule
                            t_atoms *atoms,                //atoms struct
                            rvec
                            x[],                        //coordinate of original configuration
                            int ih,                        //head id currently
                            int index_head[],            //head index
                            bool bGeo
)                        //whether using geometric center
{
  int i;
  for (
      i = 0;
      i < ih;
      i++) {
    map_res(map_r, atoms, x, index_head[i],
            NULL, bGeo, FALSE);
  }
}

void MappingAT2CG::map_conf(mapping_rule *map_r, t_atoms *atoms, rvec
x[],
                            int ih,
                            int index_head[],
                            bool bGeo
) {
  static t_symtab *symtab = NULL;
  snew(symtab, 1);
  open_symtab(symtab);

  int i;
  for (
      i = 0;
      i < ih;
      i++) {
    map_res(map_r, atoms, x, index_head[i], symtab, bGeo,
            TRUE);
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
void MappingAT2CG::map_res(mapping_rule *map_r,
                           t_atoms *atoms,
                           rvec
                           x[],
                           int head,
                           t_symtab
                           *symtab,
                           bool bGeo,
                           bool bMapName
) {
  int i, j, k, l, n, m, nres, natom, flag_over;//-geo should be added

  n = -1;
//the cg particle number of current residue
  flag_over = 0;
  rvec xcom[CG_PART_TYPE];
  real mtot;

  nres =
      resnm_to_resnr(map_r, *atoms->resinfo[atoms->atom[head].resind].name);//paticle type number of specified residue

  for (
      j = 0;
      j < map_r->cg_particle_num[nres]; j++) //scan all particles in this residue
  {

    clear_rvec(xcom[j]);
    mtot = 0;
    for (
        k = 0;
        k < map_r->atom_num_each_particle[nres][j]; k++) //scan all atoms in the cg particle
    {
      n = find_part_name_forward(map_r,
                                 atoms,
                                 head,
                                 map_r->atom_name_each_particle[nres][j][k]);//the number of atoms in the specified configure within specified cg particle.

      if (bGeo) {
        m = 1;
      }
      else {
        m = atoms->atom[n].m;
      }
      for (
          l = 0;
          l < DIM; l++)
        xcom[j][l] +=
            m * x[n][l];
      mtot +=
          m;
    }

    svmul(1 / mtot, xcom[j], xcom[j]);
//	printf("%d---%lf %lf %lf----->\n",head+j,x[head+j+1][XX],x[head+j+1][YY],x[head+j+1][ZZ]);
    copy_rvec(xcom[j], x[head + j]
    );
//		printf("%d : %d--%s  :  %s\n",strlen(map_r->cg_particle_name[nres][j]),head+j,*atoms->atomname[head+j],map_r->cg_particle_name[nres][j]);
//		printf("%d\n",(*atoms->atomname[head+j])[0]);\


    if (bMapName) {
      atoms->atomname[head + j] =
          put_symtab(symtab, map_r
              ->cg_particle_name[nres][j]);
    }
//	strcpy(*atoms->atomname[head+j],map_r->cg_particle_name[nres][j]);
  }
}

void MappingAT2CG::proc_res_name(mapping_rule *map_r, t_atoms *atoms,
                                 int
                                 isize,
                                 int index[]
) {
  static t_symtab *symtab = NULL;
  int ai, i, j, resnum;
  printf("%d", isize);

  snew(symtab, 1);
  open_symtab(symtab);

  for (
      i = 0;
      i < map_r->
          cg_res_num;
      i++) {
    for (
        j = 0;
        j < map_r->cg_num[i]; j++) {

      atoms->resinfo[map_r->residue_num[i][j] - 1].
          name = put_symtab(symtab, map_r->cg_name[i]);

    }
  }
}

void MappingAT2CG::get_head(t_atoms *atoms,
                            int isize, int index[],
                            int *ih, int index_head[]) {
  int i, ia, current;
  index_head[0] = index[0];
  *ih = 1;//number of the head index
  current = 0;//current position
  for (i = 1; i < isize; i++) {
    ia = index[i];
    if (atoms->atom[ia].resind != current) {
      *ih = *ih + 1;
      index_head[*ih - 1] = ia;
      current = atoms->atom[ia].resind;
    }
  }
}
int MappingAT2CG::resnm_to_resnr(mapping_rule *map_r, char
res[]) {
  int i;
  for (
      i = 0;
      i < map_r->
          cg_res_num;
      i++) {
    if (
        strcmp(map_r
                   ->cg_name[i], res) == 0)
      return
          i;
  }
  printf("no proper particle nr, check your rule file!");
  exit(0);
}
void MappingAT2CG::make_output_index(mapping_rule *map_r, t_atoms *atoms,
                                     int
                                     ih,
                                     int index_head[],
                                     int *io,
                                     int index_output[]
) {
  int i, j, k, t, m, n;
  int head;
  *
      io = 0;
  for (
      i = 0;
      i < ih;
      i++) {
    head = index_head[i];
    t = resnm_to_resnr(map_r, *atoms->resinfo[atoms->atom[head].resind].name);
    m = map_r->cg_particle_num[t];
    for (
        j = 0;
        j < m;
        j++) {
      n = head + j;
      index_output[*io] =
          n;
      *
          io = *io + 1;
    }
  }
}
int MappingAT2CG::find_part_name_forward(mapping_rule *map_r, t_atoms *atoms, int
from,
                                         char atomnm[]
) {
  int i;
  for (
      i = from;
      i < from + atoms->
          nr;
      i++) {
    if (
        strcmp(*atoms
            ->atomname[i], atomnm) == 0) {
      return
          i;
    }
  }
  printf("no such atom name, plz check your rule file!");
  exit(0);
}


/*! \brief
 * The main function for the analysis template.
 */
int
main(int argc, char *argv[]) {
  return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<MappingAT2CG>(argc, argv);
}
