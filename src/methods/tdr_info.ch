/*****************************************************************************
 *                                                                           *
 *          UNURAN -- Universal Non-Uniform Random number generator          *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   FILE:      tdr_info.ch                                                  *
 *                                                                           *
 *   TYPE:      continuous univariate random variate                         *
 *   METHOD:    transformed density rejection                                *
 *                                                                           *
 *   DESCRIPTION:                                                            *
 *      Routines for creating info strings.                                  *
 *                                                                           *
 *****************************************************************************
 *                                                                           *
 *   Copyright (c) 2000-2006 Wolfgang Hoermann and Josef Leydold             *
 *   Department of Statistics and Mathematics, WU Wien, Austria              *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program; if not, write to the                           *
 *   Free Software Foundation, Inc.,                                         *
 *   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                  *
 *                                                                           *
 *****************************************************************************/

/*---------------------------------------------------------------------------*/
#ifdef UNUR_ENABLE_INFO
/*---------------------------------------------------------------------------*/

void
_unur_tdr_info( struct unur_gen *gen, int help )
     /*----------------------------------------------------------------------*/
     /* create character string that contains information about the          */
     /* given generator object.                                              */
     /*                                                                      */
     /* parameters:                                                          */
     /*   gen  ... pointer to generator object                               */
     /*   help ... whether to print additional comments                      */
     /*----------------------------------------------------------------------*/
{
  struct unur_string *info = gen->infostr;
  struct unur_distr *distr = gen->distr;

  /* generator ID */
  _unur_string_append(info,"generator ID: %s\n\n",gen->genid);

  /* distribution */
  _unur_string_append(info,"distribution: %s\n",distr->name);
  _unur_string_append(info,"   type      = continuous univariate distribution\n");
  _unur_string_append(info,"   functions = PDF dPDF\n");
  _unur_string_append(info,"   domain    = (%g, %g)\n",DISTR.domain[0],DISTR.domain[1]);
  _unur_string_append(info,"   center    = %g",unur_distr_cont_get_center(distr));
  if ( !(distr->set & UNUR_DISTR_SET_CENTER) ) {
    if ( distr->set & UNUR_DISTR_SET_MODE )
      _unur_string_append(info,"  (mode)");
    else {
      if (help)
	_unur_string_append(info,"\n\t[ Hint: %s\n\t\t%s ]",
			    "You should provide the location of the mode using \"mode\" or \"center\",",
			    "if this is not a typical point of the distribution");
    }
  }
  _unur_string_append(info,"\n\n");
      
  /* method */
  _unur_string_append(info,"method: TDR (Transformed Density Rejection)\n");
  _unur_string_append(info,"   variant   = ");
  switch (gen->variant & TDR_VARMASK_VARIANT) {
  case TDR_VARIANT_GW:
    _unur_string_append(info,"GW (original Gilks & Wild)\n"); break;
  case TDR_VARIANT_PS:
    _unur_string_append(info,"PS (proportional squeeze)\n"); break;
  case TDR_VARIANT_IA:
    _unur_string_append(info,"IA (immediate acceptance)\n"); break;
  }
  /* used transformation */
  _unur_string_append(info,"   T_c(x)    = ");
  switch( gen->variant & TDR_VARMASK_T ) {
  case TDR_VAR_T_LOG:
    _unur_string_append(info,"log(x)  ... c = 0\n"); break;
  case TDR_VAR_T_SQRT:
    _unur_string_append(info,"-1/sqrt(x)  ... c = -1/2\n"); break;
  case TDR_VAR_T_POW:
    _unur_string_append(info,"-x^(%g)  ... c = %g\n",GEN->c_T,GEN->c_T); break;
  }
  /* summary hat */
  _unur_string_append(info,"   area(hat) = %g\n",GEN->Atotal);

  _unur_string_append(info,"   rejection constant ");
  if (distr->set & UNUR_DISTR_SET_PDFAREA)
    _unur_string_append(info,"= %g\n",GEN->Atotal/DISTR.area);
  else
    _unur_string_append(info,"<= %g\n",GEN->Atotal/GEN->Asqueeze);

  _unur_string_append(info,"   area ratio squeeze/hat = %g\n",
		      GEN->Asqueeze/GEN->Atotal);

  _unur_string_append(info,"   # intervals = %d\n",GEN->n_ivs);

  if (help) {
    _unur_string_append(info,"\t[ Hint: %s\n\t\t%s ]\n",
			"You can decrease the rejection constant by setting \"max_sqhratio\"",
			"closer to 1" ); 
  }



/*   _unur_string_append(info,); */

/*

max_intervals
max_sqhratio
  
  help:
  set maximum number of intervals
  set bound for ratio  Asqueeze / Atotal
  set number of starting points

*/

} /* end of _unur_tdr_info() */

/*---------------------------------------------------------------------------*/
#endif
/*---------------------------------------------------------------------------*/
