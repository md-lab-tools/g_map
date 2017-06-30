/*
  Copyright (C) 2015,2016
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


#ifndef G_MAP_G_MAP_H
#define G_MAP_G_MAP_H

#include <string>
#include <cstring>
#include <vector>

#include "gromacs/trajectoryanalysis.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/fileio/tpxio.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/index.h"
#include "gromacs/commandline/filenm.h"
#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/topology/mtop_util.h"
#include <iostream>
#include <gromacs/fileio/trxio.h>

#define CG_TYPE 20
#define CG_PART_TYPE 20
#define AA_RES_NUM 9999
#define AA_PART_NUM 999

typedef struct {
  char comment[2][255];       // comment [0]:outer [1]:inner
  int cg_res_num;             // total number of CG residues
  char cg_name[CG_TYPE][5];   //CG residue names

  int cg_num[CG_TYPE];        //corresponding AA residue numbers
  int residue_num[CG_TYPE][AA_RES_NUM];     // AA residue numbers in each CG residue

  int cg_particle_num[CG_TYPE];             // CG particle numbers for each CG residue
  char cg_particle_name[CG_TYPE][CG_PART_TYPE][5];      // CG particle names
  int atom_num_each_particle[CG_TYPE][CG_PART_TYPE];    // atom number in each CG particle
  char atom_name_each_particle[CG_TYPE][CG_PART_TYPE][AA_PART_NUM][5];    //atom name in each CG particle

} mapping_rule;

#define NFILE asize(fnm)

void load_outer_rule(std::string fn, mapping_rule *map_r);
void load_inner_rule(std::string fn, mapping_rule *map_r);

int resnm_to_resnr(mapping_rule *map_r, char *res);
void map_conf(mapping_rule *map_r, t_atoms *atoms, rvec x[], int ih, int index_head[], bool bGeo);
void map_res(mapping_rule *map_r, t_atoms *atoms, rvec x[], int head, t_symtab *symtab, bool bGeo, bool bMapName);
void map_traj(mapping_rule *map_r, t_atoms *atoms, rvec[], int ih, int index_head[], bool bGeo);
void proc_res_name(mapping_rule *map_r, t_atoms *atoms, int isize, int index[]);
void get_head(t_atoms *atoms, int isize, int index[], int *ih, int index_head[]);
void make_output_index(mapping_rule *map_r, t_atoms *atoms,
                       int ih, int index_head[],
                       int *io, int index_output[]);
int find_part_name_forward(t_atoms *atoms, int from, char atomnm[]);
#endif //G_MAP_G_MAP_H
