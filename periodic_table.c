/* Useful periodic table data  and functions */

/* Copyright (c) 2007 MJ Rutter 
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3
 * of the Licence, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, see http://www.gnu.org/licenses/
 */ 

#define P_MAX_EL 103

static char *table[P_MAX_EL+2]={"X","H","He",
  "Li","Be", "B", "C", "N", "O", "F","Ne",
  "Na","Mg","Al","Si", "P", "S","Cl","Ar",
  "K","Ca",
  "Sc","Ti", "V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
  "Ga","Ge","As","Se","Br","Kr",
  "Rb","Sr",
   "Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
  "In","Sn","Sb","Te", "I","Xe",
  "Cs","Ba",
  "La",
  "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
  "Hf","Ta", "W","Re","Os","Ir","Pt","Au","Hg",
  "Tl","Pb","Bi","Po","At","Rn",
  "Fr","Ra",
  "Ac",
  "Th","Pa", "U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lw",
0};

char* atno2sym(unsigned no){
  if (no>P_MAX_EL) return (0);
  else return (table[no]);
}

unsigned int atsym2no(char* sym){
  char *p1,*p2;
  int i;

  for(i=1;i<P_MAX_EL;i++){
    p1=sym;
    p2=table[i];
    while(*p1==' ') p1++;
    while((*p1)&&(*p2)&&(!(((*p1)^(*p2))&95))){p1++;p2++;}
    if(((*p1==0)||(*p1==' ')||(*p1==':'))&&(*p2==0)) return i;
  }

  return(0);
}
