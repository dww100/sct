C--------------------------------------------------------------------------
C
C Copyright 1981-2014 University College London
C
C Licensed under the Apache License, Version 2.0 (the "License");
C you may not use this file except in compliance with the License.
C You may obtain a copy of the License at
C
C    http://www.apache.org/licenses/LICENSE-2.0
C
C Unless required by applicable law or agreed to in writing, software
C distributed under the License is distributed on an "AS IS" BASIS,
C WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
C See the License for the specific language governing permissions and
C limitations under the License.
C
C--------------------------------------------------------------------------

      COMMON /PLOTS/ Q, I, EROR, CONSTANT, QVALP1, QVALPL, QVALF1,
     *               QVALFL, IVALP1, IVALPL, IVALF1, IVALFL, QTEMP,
     *               ITEMP, ERORU, ERORL
      DOUBLE PRECISION Q(10,512), I(10,512), EROR(10,512)
      DOUBLE PRECISION CONSTANT(10), QVALP1(10), QVALPL(10),QVALF1(10),
     *                 QVALFL(10), IVALP1(10), IVALPL(10), IVALF1(10),
     *                 IVALFL(10), ERORU(10, 512), ERORL(10,512)
      REAL QTEMP(512), ITEMP(512)
      COMMON /PLTONE/ RQMIN1, RQMAX1, RIMIN1, RIMAX1
      REAL RQMIN1, RQMAX1, RIMIN1, RIMAX1
      COMMON /PLOTAT2/ ISLINFIT, FSTFITPT, SYMBLIST, FSTPLTPT, LSTPLTPT,
     *                 SMEARD
      INTEGER ISLINFIT(10), FSTFITPT(10), SYMBLIST(10), FSTPLTPT(10),
     *        LSTPLTPT(10), SMEARD(10)
      COMMON /PLOTAT1/ PLOTCNT, PLTYPE, SYMBTYP, LINETYP, NUMPTS, PTYPE,
     *                 LSTFITPT, CURPLOT, STYLE, ISPLOT, REDRAW
      INTEGER PLOTCNT, PLTYPE(10), SYMBTYP(10), LINETYP(10), NUMPTS(10),
     *        PTYPE(10), LSTFITPT(10), CURPLOT, STYLE(10), ISPLOT(10),
     *        REDRAW(10), PLOTNUM
      COMMON /NAMES/ FILENAME
      CHARACTER FILENAME(10)*14
      INTEGER MAXCRV, MAXPLOT, MAXPTS
      PARAMETER (MAXCRV = 15)
      PARAMETER (MAXPLOT = 10)
      PARAMETER (MAXPTS = 512)
