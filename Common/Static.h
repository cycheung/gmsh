/* $Id: Static.h,v 1.4 2000-11-23 16:51:28 geuzaine Exp $ */
#ifndef _STATIC_H_
#define _STATIC_H_

/* This file defines the static structures for Gmsh. It should be
   included only once, in your 'main' file */

char TheFileName[NAME_STR_L], TheBaseFileName[NAME_STR_L];
char yyname[NAME_STR_L];
int  yyerrorstate;

Context_T   CTX ;
Mesh        M, *THEM, *LOCAL;
int         CurrentNodeNumber, CurrentSimplexNumber;
Tree_T     *EntitesVisibles = NULL;
List_T     *Post_ViewList = NULL;

double      LC, MiddleLC ;
int         LC_ORDER;
double      FACTEUR_MULTIPLICATIF=1.0, GLOBALSCALINGFACTOR=1.0;

/* Some garbage */

int  TYPBGMESH=WITHPOINTS;
int  FLAG_OLD_CIRCLE = 0 ; /* Pour David : cercles > Pi */


#endif
