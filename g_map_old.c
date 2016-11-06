/*
* $Id: template.c,v 1.4 2001/07/23 15:28:29 lindahl Exp $
* 
*                This source code is part of
* 
*                 G   R   O   M   A   C   S
* 
*          GROningen MAchine for Chemical Simulations
* 
*                        VERSION 3.0
* 
* Copyright (c) 1991-2001
* BIOSON Research Institute, Dept. of Biophysical Chemistry
* University of Groningen, The Netherlands
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
* 
* If you want to redistribute modifications, please consider that
* scientific software is very special. Version control is crucial -
* bugs must be traceable. We will be happy to consider code for
* inclusion in the official distribution, but derived work must not
* be called official GROMACS. Details are found in the README & COPYING
* files - if they are missing, get the official version at www.gromacs.org.
* 
* To help us fund GROMACS development, we humbly ask that you cite
* the papers on the package - you can find them in the top README file.
* 
* Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
* 
* And Hey:
* Gyas ROwers Mature At Cryogenic Speed
*/

/* This line is only for CVS version info */
static char *SRCID_template_c = "$Id: template.c,v 1.4 2001/07/23 15:28:29 lindahl Exp $";

#include "statutil.h"
#include "typedefs.h"
#include "smalloc.h"
#include "vec.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#define CG_TYPE 10
#define CG_PART_TYPE 20
#define AA_RES_NUM 9999
#define AA_PART_NUM 999

typedef struct{
	char comment[2][255];      						//comment [0]:outer [1]:inner
	int cg_res_num;								//total number of CG residues
	char cg_name[CG_TYPE][5];						//CG residue names
	//========================================================//
	int cg_num[CG_TYPE];							//corresponding AA residue numbers
	int residue_num[CG_TYPE][AA_RES_NUM];					//AA residue numbers in each CG residue
	//========================================================//
	int cg_particle_num[CG_TYPE];						//CG particle numbers for each CG residue
	char cg_particle_name[CG_TYPE][CG_PART_TYPE][5];			//CG particle names
	int atom_num_each_particle[CG_TYPE][CG_PART_TYPE];			//atom number in each CG particle
	char atom_name_each_particle[CG_TYPE][CG_PART_TYPE][AA_PART_NUM][5];	//atom name in each CG particle

}mapping_rule;

void proc_res_name(	mapping_rule *map_r,t_atoms *atoms,
			atom_id isize, atom_id index[])
{
	static t_symtab *symtab=NULL;
	int ai,i,j,resnum;
	printf("%d",isize);

	snew(symtab,1);
	open_symtab(symtab);

	for (i=0;i<map_r->cg_res_num;i++)
	{
		for (j=0;j<map_r->cg_num[i];j++)
		{

			atoms->resname[map_r->residue_num[i][j]-1]=put_symtab(symtab,map_r->cg_name[i]);

		}
	}
}

void get_head(	t_atoms *atoms,
		atom_id isize, atom_id index[],
		atom_id *ih,atom_id index_head[])
{
	int i,ia,current;
	index_head[0]=index[0];
	*ih=1;//number of the head index
	current=0;//current position
	for (i=1;i<isize;i++)
	{
		ia=index[i];
		if (atoms->atom[ia].resnr != current)
		{
			*ih=*ih+1;
			index_head[*ih-1]=ia;
			current=atoms->atom[ia].resnr;
		}
	}
}
int resnm_to_resnr(mapping_rule *map_r,char res[])
{
	int i;
	for (i=0;i<map_r->cg_res_num;i++)
	{
		if (strcmp(map_r->cg_name[i],res)==0)
		return i;
	}
	printf("no proper particle nr, check your rule file!");
	exit(0);
}
void make_output_index(	mapping_rule *map_r, t_atoms *atoms,
			atom_id ih, atom_id index_head[],
			atom_id *io, atom_id index_output[])
{
	int i,j,k,t,m,n;
	int head;
	*io=0;
	for (i=0;i<ih;i++)
	{
		head=index_head[i];
		t=resnm_to_resnr(map_r,*atoms->resname[atoms->atom[head].resnr]);
		m=map_r->cg_particle_num[t];
		for (j=0;j<m;j++)
		{
			n=head+j;
			index_output[*io]=n;
			*io=*io+1;
		}
	}
}
int find_part_name_forward(mapping_rule *map_r, t_atoms *atoms,int from,char atomnm[])
{
	int i;
	for (i=from;i<from+atoms->nr;i++)
	{
		if (strcmp(*atoms->atomname[i],atomnm)==0)
		{
			return i;
		}
	}
	printf("no such atom name, plz check your rule file!");
	exit(0);
}
void map_res(mapping_rule *map_r,t_atoms *atoms,rvec x[],atom_id head,t_symtab *symtab, bool bGeo,bool bMapName)
{
	int i,j,k,l,n,m,nres,natom,flag_over;//-geo should be added

	n=-1;//the cg particle number of current residue
	flag_over=0;
	rvec xcom[CG_PART_TYPE];
	real mtot;
	
	nres=resnm_to_resnr(map_r,*atoms->resname[atoms->atom[head].resnr]);//paticle type number of specified residue
	
	for (j=0;j<map_r->cg_particle_num[nres];j++) //scan all particles in this residue
	{
		
		clear_rvec(xcom[j]);
		mtot=0;
		for(k=0;k<map_r->atom_num_each_particle[nres][j];k++) //scan all atoms in the cg particle
		{
			n=find_part_name_forward(map_r,atoms,head,map_r->atom_name_each_particle[nres][j][k]);//the number of atoms in the specified configure within specified cg particle.	
			
			if (bGeo)
			{
				m=1;
			}
			else
			{
				m=atoms->atom[n].m;
			}
			for(l=0; l<DIM; l++)
				xcom[j][l] += m*x[n][l];
			mtot += m;	
		}
		
		svmul(1/mtot,xcom[j],xcom[j]);
	//	printf("%d---%lf %lf %lf----->\n",head+j,x[head+j+1][XX],x[head+j+1][YY],x[head+j+1][ZZ]);
		copy_rvec(xcom[j],x[head+j]);
//		printf("%d : %d--%s  :  %s\n",strlen(map_r->cg_particle_name[nres][j]),head+j,*atoms->atomname[head+j],map_r->cg_particle_name[nres][j]);
//		printf("%d\n",(*atoms->atomname[head+j])[0]);\


		if (bMapName) {atoms->atomname[head+j]=put_symtab(symtab,map_r->cg_particle_name[nres][j]);}
	//	strcpy(*atoms->atomname[head+j],map_r->cg_particle_name[nres][j]);
	}	
}
void map_traj(mapping_rule *map_r,t_atoms *atoms, rvec x[], 
		atom_id ih, atom_id index_head[], bool bGeo)
{
	static t_symtab *symtab=NULL;
//	snew(symtab,1);
  //      open_symtab(symtab);
	int i;
	for (i=0;i<ih;i++)
	{
		map_res(map_r,atoms,x,index_head[i],NULL,bGeo,FALSE);
	}
//	sfree(symtab);


	
}
void map_conf(mapping_rule *map_r,t_atoms *atoms, rvec x[],
		atom_id ih, atom_id index_head[], bool bGeo)
{
        static t_symtab *symtab=NULL;
        snew(symtab,1);
        open_symtab(symtab);

	int i;
	for(i=0;i<ih;i++)
	{
		map_res(map_r,atoms,x,index_head[i],symtab,bGeo,TRUE);	
	}
	sfree(symtab);

}
/*the loading function is not robost enough
so do please check your rule file carefully!*/
void load_outer_rule(char *fn,mapping_rule *map_r)
{
	FILE *fp;
	int i,j;
	char temp[256];
	fp = fopen(fn,"r");
	fgets(temp,255,fp);
	printf("\n-loading outer_rule-: %s \n",temp);
	strcpy(map_r->comment[0],temp);
	fgets(temp,255,fp);
	sscanf(temp,"%d",&map_r->cg_res_num);
	for (i=0;i<map_r->cg_res_num;i++)
	{
		fgets(temp,255,fp);
		sscanf(temp,"%s%d",map_r->cg_name[i],&map_r->cg_num[i]);

		for (j=0;j<map_r->cg_num[i];j++)
		{
			fgets(temp,255,fp);
			sscanf(temp,"%d",&map_r->residue_num[i][j]);
		}
	}
	fclose(fp);

}
void load_inner_rule(char *fn,mapping_rule *map_r)
{

	FILE *fp;
	int i,j,k;
	char temp[256];
	fp = fopen(fn,"r");

	fgets(temp,100,fp);

	printf("\n-loading inner_rule-: %s \n",temp);

	strcpy(map_r->comment[1],temp);

	fgets(temp,100,fp);

	int temp_num;
	sscanf(temp,"%d",&temp_num);

	if (temp_num!=map_r->cg_res_num)
	{
		printf("Check your rule file! The residue number is not equal\n");
		exit(0);
	}

	for (i=0;i<map_r->cg_res_num;i++)
	{
		fgets(temp,255,fp);
		sscanf(temp,"%s%d",map_r->cg_name[i],&map_r->cg_particle_num[i]);

		for (j=0;j<map_r->cg_particle_num[i];j++)
		{
			fgets(temp,255,fp);

			sscanf(temp,"%s%d",map_r->cg_particle_name[i][j],&map_r->atom_num_each_particle[i][j]);

			for (k=0;k<map_r->atom_num_each_particle[i][j];k++)
			{
				fgets(temp,100,fp);
				sscanf(temp,"%s",map_r->atom_name_each_particle[i][j][k]);
			}
		}
	}

}
int main(int argc,char *argv[])
{
	static char *desc[] = {
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

	static int n=1;
	static bool bGeo=FALSE;
	int ngrps=1;
	t_topology top;
	char       title[STRLEN];
	t_trxframe fr;
	rvec       *xtop;
	matrix     box;
	int        status;
	int        flags = TRX_READ_X;
	int        isize_old,isize_new,isize_head,isize_output;
	char       *grpname_old;
	atom_id    *index_old,*index_new,*index_head,*index_output;
	FILE       *fp_confout;
	bool	   bIndex,bOutRule,bInnerRule;

	/*parameter input*/
	t_pargs pa[] = {
		{ "-geo", FALSE, etBOOL, {&bGeo}, 
		"mapping with geometric center (default is center of mass)" }
	};
	/*file I/O part*/
	t_filenm fnm[] = {
		{ efTPS, "-s", "conf", ffREAD },         /* this is the input topology */
		{ efTRX, "-f", "traj", ffREAD },         /* this is the input trajectory */
		{ efNDX, "-n", "index", ffREAD },	/* this is the grp to be coarse grained*/
		{ efTRX, "-o", "trajout", ffWRITE },        /* this is the output topology */
		{ efSTO, "-c", "confout", ffWRITE },        /* this is the output trajectory */
		{ efDAT, "-or","outer", ffREAD },
		{ efDAT, "-ir","inner", ffREAD },
	};

#define NFILE asize(fnm)

	CopyRight(stderr,argv[0]);
	parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,
		NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);


	bIndex    = opt2bSet("-n",NFILE,fnm);



	read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xtop,NULL,box,TRUE);
	
//	sfree(xtop);
	get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&isize_old,&index_old,&grpname_old);
	
	snew(index_head,isize_old);
	snew(index_output,isize_old);

	/* The first time we read data is a little special */
	read_first_frame(&status,ftp2fn(efTRX,NFILE,fnm),&fr,flags);

	//pre proc
	fp_confout=open_trx(opt2fn("-o",NFILE,fnm),"w");

	mapping_rule map_r;

	load_outer_rule(opt2fn("-or",NFILE,fnm),&map_r);	/*load the outer mapping rule*/
	load_inner_rule(opt2fn("-ir",NFILE,fnm),&map_r);	/*load the inner mapping rule*/

	proc_res_name(&map_r,&top.atoms,isize_old, index_old);
	get_head(&top.atoms,isize_old, index_old,&isize_head,index_head);
	make_output_index(&map_r,&top.atoms,isize_head,index_head,&isize_output,index_output);
	
	printf("Loading rule file finished!\n\n");
	
	do {
		/* coordinates are available in the vector fr.x
		* you can find this and all other structures in
		* the types directory under the gromacs include dir.
		* Note how flags determines wheter to read x/v/f!
		*/
 		map_traj(&map_r,&top.atoms,fr.x,isize_head,index_head,bGeo);
		write_trxframe_indexed(fp_confout,&fr,isize_output,index_output);//need to be changed into new
		
	} while(read_next_frame(status,&fr));

	close_trj(status);

	map_conf(&map_r,&top.atoms,xtop,isize_head,index_head,bGeo);

	printf("%s\n",title);

	write_sto_conf_indexed(ftp2fn(efSTO,NFILE,fnm),title,&top.atoms,xtop,NULL,box,isize_output,index_output);//need to be changed into new

	thanx(stderr);

	return 0;
}

