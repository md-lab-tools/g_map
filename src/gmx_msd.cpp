/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
//#include "gmxpre.h"

#include <cmath>
#include <cstring>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
//#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/gmxana/gstat.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/pbcutil/rmpbc.h"
//#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/index.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"


int main(int argc, char *argv[]) {
  const char *desc[] = {
      "this is a small test program meant to serve as a template ",
      "when writing your own analysis tools. The advantage of ",
      "using gromacs for this is that you have access to all ",
      "information in the topology, and your program will be ",
      "able to handle all types of coordinates and trajectory ",
      "files supported by gromacs. Go ahead and try it! ",
      "This test version just writes the coordinates of an ",
      "arbitrary atom to standard out for each frame. You can ",
      "select which atom you want to examine with the -n argument."
  };

  /*parameter input*/
  bool bGeo = FALSE;
  t_pargs pa[] = {
      {"-geo", FALSE, etBOOL, {&bGeo},
       "mapping with geometric center (default is center of mass)"}
  };

  t_filenm fnm[] = {
      {efTPS, "-s", "topol", ffREAD},         /* this is the input topology */
      {efTRX, "-f", "traj", ffREAD},         /* this is the input trajectory */
      {efNDX, "-n", "index", ffREAD},    /* this is the grp to be coarse grained*/
      {efTRX, "-o", "trajout", ffWRITE},        /* this is the output topology */
      {efSTO, "-c", "confout", ffWRITE},        /* this is the output trajectory */
      {efDAT, "-or", "outer", ffREAD},
      {efDAT, "-ir", "inner", ffREAD},
  };
#define NFILE asize(fnm)

  t_topology top;
  int ePBC;

  matrix box;
  const char *trx_file, *tps_file, *ndx_file;
  int isize_old, isize_new, isize_head, isize_output;
  int *index_old, *index_new, *index_head, *index_output;
  char *grpname_old;
  t_trxstatus *fp_confout;
  rvec *xdum;
  gmx_bool bTop;
  int axis, type;
  real dim_factor;
  gmx_output_env_t *oenv;

  if (!parse_common_args(&argc, argv,
                         PCA_CAN_VIEW | PCA_CAN_BEGIN | PCA_CAN_END | PCA_TIME_UNIT,
                         NFILE, fnm, asize(pa), pa, asize(desc), desc, 0, NULL, &oenv)) {
    return 0;
  }

  trx_file = ftp2fn_null(efTRX, NFILE, fnm);
  tps_file = ftp2fn_null(efTPS, NFILE, fnm);
  ndx_file = ftp2fn_null(efNDX, NFILE, fnm);

  bTop = read_tps_conf(tps_file, &top, &ePBC, &xdum, NULL, box, FALSE);
  if (!bTop) {
    gmx_fatal(FARGS, "Could not read a topology from %s. Try a tpr file instead.", tps_file);
  }
  get_index(&top.atoms, ndx_file, 1, &isize_old, &index_old, &grpname_old);
  snew(index_head, isize_old);
  snew(index_output, isize_old);

  /* The first time we read data is a little special */

  view_all(oenv, NFILE, fnm);

  return 0;
}
