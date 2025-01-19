%let pgm=utl-better-fit-with-two-paramater-beta-then-four-than-four-polynomial-on-unit-square;

%stop_submission;


Better fit with two paramater beta then four than four polynomial on unit square


HI RES GRAPHIC OUTPUT

github
https://tinyurl.com/m2xtstfj
https://github.com/rogerjdeangelis/utl-better-fit-with-two-paramater-beta-then-four-than-four-polynomial-on-unit-square/blob/main/pltcum.png


OBSERVATIONS
============
Note the beta CDF is only defined on the x~(0,1) ans y~[0,1] so you can't easily forecast, outside the square.
Also note the endpoints are forced to fit given they are mapped to 0 and 1
Fewer parameters is often better?
If you think the wobble in the graphh are real,  we could introduce trig functions (sin(c*x)**2?

As I sine note python sympy can be very useful


  PROCESS
     1 INPUT x,y pairs
     2 map to unit square
     3 non linear beta fit
     4 map back to original units amd compute NSE in original units
     5 stackoverflow multiple linear regression
     6 join all three results roger stackoverflow and original data
     7 plot all three original y and roger and stackoverflow predictions
     8 hi res graphic


Related repos
----------------------------------------------------------------------------------------------------------------------------
https://github.com/rogerjdeangelis/utl-using-beta-cdf-and-pdf-unitsquare-to-fit-non-linear-equations-and-simple-polynomials
https://github.com/rogerjdeangelis/utl-using-the-unitsquare-framework-and-cumulative-beta-to-fit-a-non-linear-equation
https://github.com/rogerjdeangelis/utl-general-framework-for-fitting-non-linear-equation-with-transcendental-functions-on-a-unit-square
https://github.com/rogerjdeangelis/utl-piecewise-linear-fit-compared-to-unit-square-beta-cdf-non-linear-fit-two-vs-four-parameters


/*               _     _
 _ __  _ __ ___ | |__ | | ___ _ __ ___
| `_ \| `__/ _ \| `_ \| |/ _ \ `_ ` _ \
| |_) | | | (_) | |_) | |  __/ | | | | |
| .__/|_|  \___/|_.__/|_|\___|_| |_| |_|
|_|
*/

/***********************************************************************************************************************************/
/*                                            |                                            |                                       */
/*                INPUT                       |        PROCESS                             |              OUTPUT                   */
/*                =====                       |        =======                             |              ======                   */
/*                                            |                                            |                                       */
/*                    X                       |                                            |                                       */
/*     200       450       700       950      | MAP TO UNIT SQUARE                         |    --+------+------+------+------+--  */
/* Y  --+---------+---------+---------+---- Y | ==================                         | Y  |                               |  */
/*    |                                   |   |                                            |    |  All Three Y values           |  */
/* 15 + Will fit the endpoints exactly    + 15| proc sql;                                  | 12 +                               +12*/
/*    |                                   |   |   create                                   |    |-- PERFECT FIT (950,11.7)-> X -|  */
/*    | Map to the unit square and fit    |   |      table mapped as                       |    |                               |  */
/*    | y = cdf('BETA', x , a, b)         |   |   select                                   |    |  r=roger                      |  */
/*    |                                   |   |     xs as xorg                             |    |  o=original                   |  */
/*    |  X  Y                (950,11.7) X |   |    ,ys as yorg                             |    |  s-stackoverflow              |  */
/*    | 200 0                             |   |    ,min(xs) as  minxs                      |    |                            s  |  */
/* 10 + 250 1.05                        * + 10|    ,max(xs) as  maxxs                      | 10 +                MSE            +10*/
/*    | 370 2.20                          |   |    ,min(ys) as  minys                      |    |                =====          |  */
/*    | 470 3.33                      *   |   |    ,max(ys) as  maxys                      |    |  ROGER         0.172      s   |  */
/*    | 520 4.31                     *    |   |    ,(-min(xs)+xs)/(max(xs)-min(xs)) as xs  |    |  STACKOVERFLOW 0.815      o   |  */
/*    | 590 4.70                   **     |   |    ,(-min(xs)+xs)/(max(xs)-min(xs)) as xs  |    |                               |  */
/*    | 660 5.41                 **       |   |    ,(-min(ys)+ys)/(max(ys)-min(ys)) as ys  |    |  ROGER EQUATION           r   |  */
/*    | 800 5.80              ***         |   |   from                                     |    |                               |  */
/*  5 + 880 6.70           ***            +  5|      sd1.have                              |  8 +  Y = cdf('BETA'          s    + 8*/
/*    | 925 9.05       ****               |   | ;quit;                                     |    |    ,(x-200)/750               |  */
/*    | 950 11.7   ****                   |   |                  UNIT SQUARE               |    |    ,0.4938,0.2961)*11.7  r    |  */
/*    |        ***                        |   |                ================            |    |                               |  */
/*    |     ***                           |   | XORG   YORG       XS       YS              |    |                               |  */
/*    |   **                              |   |                                            |    |  STACKOVERFLOW EQU       o    |  */
/*    |  *             SMOOTHED           |   |  200   0.00    0.00000  0.00000            |    |                       s       |  */
/*  0 + X (200,0)                         +  0|  250   1.05    0.06667  0.08974            |  6 +  Y= -10.0693                  + 6*/
/*    |                                   |   |  370   2.20    0.22667  0.18803            |    |   +7.4654E-8*x**3      r      |  */
/*    --+---------+---------+---------+----   |  470   3.33    0.36000  0.28462            |    |   -0.000123x**2     o         |  */
/*     200       450       700       950      |  520   4.31    0.42667  0.36838            |    |   + 0.07033  *x               |  */
/*                    X                       |  590   4.70    0.52000  0.40171            |    |                   r           |  */
/*                                            |  660   5.41    0.61333  0.46239            |    |                  os           |  */
/*   options validvarname=upcase;             |  800   5.80    0.80000  0.49573            |    |                o s            |  */
/*   libname sd1 "d:/sd1";                    |  880   6.70    0.90667  0.57265            |  4 +                s              + 4*/
/*   data sd1.have;                           |  925   9.05    0.96667  0.77350            |    |              s r              |  */
/*    input xs ys;                            |  950  11.70    1.00000  1.00000            |    |              o                |  */
/*   cards4;                                  |                                            |    |                               |  */
/*   200 0                                    |                                            |    |           s                   |  */
/*   250 1.05                                 |  FIT USING BETA CDF                        |    |           r                   |  */
/*   370 2.20                                 |  =================                         |    |           o                   |  */
/*   470 3.33                                 |                                            |  2 +                               + 2*/
/*   520 4.31                                 |  proc nlin data=mapped                     |    |                               |  */
/*   590 4.70                                 |    method=MARQUARDT CONVERGE=0.0001;       |    |        r                      |  */
/*   660 5.41                                 |    parms A = 1                             |    |        o                      |  */
/*   800 5.80                                 |          B = 1;                            |    |                               |  */
/*   880 6.70                                 |    model ys = cdf('BETA', xs, a, b);       |    |                               |  */
/*   925 9.05                                 |    output                                  |    |                               |  */
/*   950 11.70                                |       out=results                          |  0 +------ X --PERFECT FIT(0,200) -+ 0*/
/*   ;;;;                                     |       p=predicted                          |    |       s                       |  */
/*   run;quit;                                |       r=residual;                          |    |                               |  */
/*                                            |  run;quit;                                 |    |                               |  */
/*                                            |                                            |    |                               |  */
/*                                            |                                            |    |                               |  */
/*                                            |  MAP BACK ORGINAL UNITS AND COMPUTE MSE    |    |                               |  */
/*                                            |  =======================================   | -2 +                               +-2*/
/*                                            |                                            |    |                               |  */
/*                                            |  proc sql;                                 |    --+------+------+------+------+--  */
/*                                            |    create                                  |      0     250    500    750   1000   */
/*                                            |       table roger as                       |                                       */
/*                                            |    select                                  |                  X_ORG                */
/*                                            |      *                                     |                                       */
/*                                            |     ,avg(sse) as mse                       |                                       */
/*                                            |    from (                                  |       Y= -10.06925                    */
/*                                            |      select                                |        +7.465462E-8*x**3              */
/*                                            |        yorg                                |        +-0.00012226*x**2              */
/*                                            |       ,xs*(max(xorg)-min(xorg))            |        +    0.07033*x                 */
/*                                            |          +min(xorg)        as xorg         |                                       */
/*                                            |       ,predicted*(max(yorg)-min(yorg))     |                                       */
/*                                            |          +min(yorg) as ys                  |                                       */
/*                                            |       ,(calculated ys - yorg)**2 as sse    |                                       */
/*                                            |      from                                  |                                       */
/*                                            |        results )                           |                                       */
/*                                            |                                            |                                       */
/*                                            |  ;quit;                                    |                                       */
/*                                            |                                            |                                       */
/*                                            |  STACKOVERFLOW LINEAR REGRESSION OVERFIT   |                                       */
/*                                            |  ======================================    |                                       */
/*                                            |                                            |                                       */
/*                                            |  data fix;                                 |                                       */
/*                                            |   set sd1.have;                            |                                       */
/*                                            |    xs3=xs**3;                              |                                       */
/*                                            |    xs2=xs**2;                              |                                       */
/*                                            |  run;quit;                                 |                                       */
/*                                            |                                            |                                       */
/*                                            |  proc reg data=fix;                        |                                       */
/*                                            |    model ys = xs3 xs2 xs ;                 |                                       */
/*                                            |    output out=linreg p=predlin;            |                                       */
/*                                            |  run;quit;                                 |                                       */
/*                                            |                                            |                                       */
/*                                            |  JOIN ROGER, ORIGINAL AND STACKOVERFLOW    |                                       */
/*                                            |  ======================================    |                                       */
/*                                            |                                            |                                       */
/*                                            |  proc sql;                                 |                                       */
/*                                            |    create                                  |                                       */
/*                                            |      table allTre as                       |                                       */
/*                                            |    select                                  |                                       */
/*                                            |      org.xs  as x_org                      |                                       */
/*                                            |     ,org.ys  as Y_org                      |                                       */
/*                                            |     ,rog.ys  as y_rog                      |                                       */
/*                                            |     ,lin.predlin as y_stack                |                                       */
/*                                            |    from                                    |                                       */
/*                                            |      sd1.have as org                       |                                       */
/*                                            |    left join                               |                                       */
/*                                            |      roger as rog on org.xs=rog.xorg       |                                       */
/*                                            |    left join                               |                                       */
/*                                            |      linreg as lin on org.xs=lin.xs        |                                       */
/*                                            |  ;quit                                     |                                       */
/*                                            |                                            |                                       */
/*                                            |   PLOT ALL THREE                           |                                       */
/*                                            |   ==============                           |                                       */
/*                                            |                                            |                                       */
/*                                            |   options ls=64 ps=64;                     |                                       */
/*                                            |   proc plot data=allTre(                   |                                       */
/*                                            |    rename=y_org=                           |                                       */
/*                                            |    y12345678901234567890123456789);        |                                       */
/*                                            |    plot                                    |                                       */
/*                                            |     y12345678901234567890123456789         |                                       */
/*                                            |       *x_org='o'                           |                                       */
/*                                            |     y_stack*x_org='s'                      |                                       */
/*                                            |     y_rog*x_org='r'                        |                                       */
/*                                            |    /box overlay vref=0 11.7 ;              |                                       */
/*                                            |    run;quit;                               |                                       */
/*                                            |                                            |                                       */
/***********************************************************************************************************************************/

/*   _                   _
/ | (_)_ __  _ __  _   _| |_
| | | | `_ \| `_ \| | | | __|
| | | | | | | |_) | |_| | |_
|_| |_|_| |_| .__/ \__,_|\__|
            |_|
*/

libname sd1 "d:/sd1";

proc datasets lib=work kill nodetails nolist;
run;quit;

proc datasets lib=sd1 kill  nodetails nolist;
run;quit;

options validvarname=upcase;
data sd1.have;
 input xs ys;
cards4;
200 0
250 1.05
370 2.20
470 3.33
520 4.31
590 4.70
660 5.41
800 5.80
880 6.70
925 9.05
950 11.70
;;;;
run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*     XS      YS                                                                                                         */
/*                                                                                                                        */
/*    200     0.00                                                                                                        */
/*    250     1.05                                                                                                        */
/*    370     2.20                                                                                                        */
/*    470     3.33                                                                                                        */
/*    520     4.31                                                                                                        */
/*    590     4.70                                                                                                        */
/*    660     5.41                                                                                                        */
/*    800     5.80                                                                                                        */
/*    880     6.70                                                                                                        */
/*    925     9.05                                                                                                        */
/*    950    11.70                                                                                                        */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*___                                       _ _
|___ \  _ __ ___   __ _ _ __    _   _ _ __ (_) |_   ___  __ _ _   _  __ _ _ __ ___
  __) || `_ ` _ \ / _` | `_ \  | | | | `_ \| | __| / __|/ _` | | | |/ _` | `__/ _ \
 / __/ | | | | | | (_| | |_) | | |_| | | | | | |_  \__ \ (_| | |_| | (_| | | |  __/
|_____||_| |_| |_|\__,_| .__/   \__,_|_| |_|_|\__| |___/\__, |\__,_|\__,_|_|  \___|
                       |_|                                 |_|
*/
proc sql;
  create
     table mapped as
  select
    xs as xorg
   ,ys as yorg
   ,min(xs) as  minxs
   ,max(xs) as  maxxs
   ,min(ys) as  minys
   ,max(ys) as  maxys
   ,(-min(xs)+xs)/(max(xs)-min(xs)) as xs
   ,(-min(xs)+xs)/(max(xs)-min(xs)) as xs
   ,(-min(ys)+ys)/(max(ys)-min(ys)) as ys
  from
     sd1.have
;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/* Obs  XORG   YORG  MINXS  MAXXS  MINYS  MAXYS     XS       YS                                                           */
/*                                                                                                                        */
/*   1   200   0.00   200    950     0     11.7  0.00000  0.00000                                                         */
/*   2   250   1.05   200    950     0     11.7  0.06667  0.08974                                                         */
/*   3   370   2.20   200    950     0     11.7  0.22667  0.18803                                                         */
/*   4   470   3.33   200    950     0     11.7  0.36000  0.28462                                                         */
/*   5   520   4.31   200    950     0     11.7  0.42667  0.36838                                                         */
/*   6   590   4.70   200    950     0     11.7  0.52000  0.40171                                                         */
/*   7   660   5.41   200    950     0     11.7  0.61333  0.46239                                                         */
/*   8   800   5.80   200    950     0     11.7  0.80000  0.49573                                                         */
/*   9   880   6.70   200    950     0     11.7  0.90667  0.57265                                                         */
/*  10   925   9.05   200    950     0     11.7  0.96667  0.77350                                                         */
/*  11   950  11.70   200    950     0     11.7  1.00000  1.00000                                                         */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*____                       _ _                        _          _           __ _ _
|___ /   _ __   ___  _ __   | (_)_ __   ___  __ _ _ __ | |__   ___| |_ __ _   / _(_) |_
  |_ \  | `_ \ / _ \| `_ \  | | | `_ \ / _ \/ _` | `__|| `_ \ / _ \ __/ _` | | |_| | __|
 ___) | | | | | (_) | | | | | | | | | |  __/ (_| | |   | |_) |  __/ || (_| | |  _| | |_
|____/  |_| |_|\___/|_| |_| |_|_|_| |_|\___|\__,_|_|   |_.__/ \___|\__\__,_| |_| |_|\__|
*/

proc nlin data=mapped
  method=MARQUARDT CONVERGE=0.0001;
  parms A = 1
        B = 1;
  model ys = cdf('BETA', xs, a, b);
  output
     out=results
     p=predicted
     r=residual;
run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*                                                                                                                        */
/*                                    Sum of        Mean               Approx                                             */
/*  Source                    DF     Squares      Square    F Value    Pr > F                                             */
/*                                                                                                                        */
/*  Model                      2      2.7935      1.3967     911.36    <.0001                                             */
/*  Error                      9      0.0138     0.00153                                                                  */
/*  Uncorrected Total         11      2.8073                                                                              */
/*                                                                                                                        */
/*                                                                                                                        */
/*                                Approx                                                                                  */
/*  Parameter      Estimate    Std Error    Approximate 95% Confidence Limits                                             */
/*                                                                                                                        */
/*  A                0.4938       0.0635      0.3502      0.6373                                                          */
/*  B                0.2961       0.0346      0.2178      0.3744                                                          */
/*                                                                                                                        */
/*                                                                                                                        */
/* XORG     YORG    MINXS    MAXXS    MINYS    MAXYS       XS         YS      PREDICTED     RESIDUAL                      */
/*                                                                                                                        */
/*  200     0.00     200      950       0       11.7    0.00000    0.00000     0.00000      0.000000                      */
/*  250     1.05     200      950       0       11.7    0.06667    0.08974     0.11675     -0.027007                      */
/*  370     2.20     200      950       0       11.7    0.22667    0.18803     0.22286     -0.034823                      */
/*  470     3.33     200      950       0       11.7    0.36000    0.28462     0.29170     -0.007081                      */
/*  520     4.31     200      950       0       11.7    0.42667    0.36838     0.32457      0.043809                      */
/*  590     4.70     200      950       0       11.7    0.52000    0.40171     0.37089      0.030818                      */
/*  660     5.41     200      950       0       11.7    0.61333    0.46239     0.41942      0.042970                      */
/*  800     5.80     200      950       0       11.7    0.80000    0.49573     0.53551     -0.039782                      */
/*  880     6.70     200      950       0       11.7    0.90667    0.57265     0.63445     -0.061805                      */
/*  925     9.05     200      950       0       11.7    0.96667    0.77350     0.73247      0.041030                      */
/*  950    11.70     200      950       0       11.7    1.00000    1.00000     1.00000      0.000000                      */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*  _     _                _      _                     _       _             _              _ _
| || |   | |__   __ _  ___| | __ | |_ ___     ___  _ __(_) __ _(_)_ __   __ _| | _   _ _ __ (_) |_ ___
| || |_  | `_ \ / _` |/ __| |/ / | __/ _ \   / _ \| `__| |/ _` | | `_ \ / _` | || | | | `_ \| | __/ __|
|__   _| | |_) | (_| | (__|   <  | || (_) | | (_) | |  | | (_| | | | | | (_| | || |_| | | | | | |_\__ \
   |_|   |_.__/ \__,_|\___|_|\_\  \__\___/   \___/|_|  |_|\__, |_|_| |_|\__,_|_| \__,_|_| |_|_|\__|___/
                                                          |___/
*/
/*---- COMPUTE MSE ----*/

proc sql;
  create
     table roger as
  select
    *
   ,avg(sse) as mse
  from (
    select
      yorg
     ,xs*(max(xorg)-min(xorg))
        +min(xorg)        as xorg
     ,predicted*(max(yorg)-min(yorg))
        +min(yorg) as ys
     ,(calculated ys - yorg)**2 as sse
    from
      results )

;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*     YORG    XORG       YS        SSE        MSE                                                                        */
/*                                                                                                                        */
/*     0.00     200     0.0000    0.00000    0.17165                                                                      */
/*     1.05     250     1.3660    0.09984    0.17165                                                                      */
/*     2.20     370     2.6074    0.16599    0.17165                                                                      */
/*     3.33     470     3.4128    0.00686    0.17165                                                                      */
/*     4.31     520     3.7974    0.26273    0.17165                                                                      */
/*     4.70     590     4.3394    0.13001    0.17165                                                                      */
/*     5.41     660     4.9073    0.25276    0.17165                                                                      */
/*     5.80     800     6.2655    0.21665    0.17165                                                                      */
/*     6.70     880     7.4231    0.52289    0.17165                                                                      */
/*     9.05     925     8.5700    0.23044    0.17165                                                                      */
/*    11.70     950    11.7000    0.00000    0.17165                                                                      */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*___        _             _                        __ _                                                  _
| ___|   ___| |_ __ _  ___| | _______   _____ _ __ / _| | _____      __  _ __ ___  __ _ _ __ ___  ___ ___(_) ___  _ __
|___ \  / __| __/ _` |/ __| |/ / _ \ \ / / _ \ `__| |_| |/ _ \ \ /\ / / | `__/ _ \/ _` | `__/ _ \/ __/ __| |/ _ \| `_ \
 ___) | \__ \ || (_| | (__|   < (_) \ V /  __/ |  |  _| | (_) \ V  V /  | | |  __/ (_| | | |  __/\__ \__ \ | (_) | | | |
|____/  |___/\__\__,_|\___|_|\_\___/ \_/ \___|_|  |_| |_|\___/ \_/\_/   |_|  \___|\__, |_|  \___||___/___/_|\___/|_| |_|
                                                                                  |___/
*/

data fix;
 set sd1.have;
  xs3=xs**3;
  xs2=xs**2;
run;quit;

proc reg data=fix;
  model ys = xs3 xs2 xs ;
  output out=linreg p=predlin;
run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*                               Analysis of Variance                                                                     */
/*                                                                                                                        */
/*                                      Sum of           Mean                                                             */
/*  Source                   DF        Squares         Square    F Value    Pr > F                                        */
/*                                                                                                                        */
/*  Model                     3      111.03013       37.01004      45.40    <.0001                                        */
/*  Error                     7        5.70683        0.81526                                                             */
/*  Corrected Total          10      116.73696                                                                            */
/*                                                                                                                        */
/*                                                                                                                        */
/*  Root MSE              0.90292    R-Square     0.9511                                                                  */
/*  Dependent Mean        4.93182    Adj R-Sq     0.9302                                                                  */
/*  Coeff Var            18.30802                                                                                         */
/*                                                                                                                        */
/*                                                                                                                        */
/*                          Parameter Estimates                                                                           */
/*                                                                                                                        */
/*                       Parameter       Standard                                                                         */
/*  Variable     DF       Estimate          Error    t Value    Pr > |t|                                                  */
/*                                                                                                                        */
/*  Intercept     1      -10.06925        4.09274      -2.46      0.0434                                                  */
/*  XS3           1    7.465462E-8    3.011302E-8       2.48      0.0423                                                  */
/*  XS2           1    -0.00012226     0.00005193      -2.35      0.0508                                                  */
/*  XS            1        0.07033        0.02695       2.61      0.0349                                                  */
/*                                                                                                                        */
/*                                                                                                                        */
/*  Obs     XS      YS        XS3         XS2      PREDLIN                                                                */
/*                                                                                                                        */
/*    1    200     0.00      8000000     40000     -0.2958                                                                */
/*    2    250     1.05     15625000     62500      1.0393                                                                */
/*    3    370     2.20     50653000    136900      2.9981                                                                */
/*    4    470     3.33    103823000    220900      3.7310                                                                */
/*    5    520     4.31    140608000    270400      3.9420                                                                */
/*    6    590     4.70    205379000    348100      4.2012                                                                */
/*    7    660     5.41    287496000    435600      4.5573                                                                */
/*    8    800     5.80    512000000    640000      6.1743                                                                */
/*    9    880     6.70    681472000    774400      8.0212                                                                */
/*   10    925     9.05    791453125    855625      9.4663                                                                */
/*   11    950    11.70    857375000    902500     10.4151                                                                */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*__       _       _                                              _   __ _ _   _           _
 / /_     (_) ___ (_)_ __    _ __ __ ___      __   __ _ _ __   __| | / _(_) |_| |_ ___  __| |
| `_ \    | |/ _ \| | `_ \  | `__/ _` \ \ /\ / /  / _` | `_ \ / _` || |_| | __| __/ _ \/ _` |
| (_) |   | | (_) | | | | | | | | (_| |\ V  V /  | (_| | | | | (_| ||  _| | |_| ||  __/ (_| |
 \___/   _/ |\___/|_|_| |_| |_|  \__,_| \_/\_/    \__,_|_| |_|\__,_||_| |_|\__|\__\___|\__,_|
        |__/
*/

proc sql;
  create
    table allTre as
  select
    org.xs  as x_org
   ,org.ys  as Y_org
   ,rog.ys  as y_rog
   ,lin.predlin as y_stack
  from
    sd1.have as org
  left join
    roger as rog on org.xs=rog.xorg
  left join
    linreg as lin on org.xs=lin.xs
;quit

/**************************************************************************************************************************/
/*                                                                                                                        */
/*  Obs    X_ORG    Y_ORG     Y_ROG      Y_STACK                                                                          */
/*                                                                                                                        */
/*    1     200      0.00     0.0000     -0.2958                                                                          */
/*    2     250      1.05     1.3660      1.0393                                                                          */
/*    3     370      2.20     2.6074      2.9981                                                                          */
/*    4     470      3.33     3.4128      3.7310                                                                          */
/*    5     520      4.31     3.7974      3.9420                                                                          */
/*    6     590      4.70     4.3394      4.2012                                                                          */
/*    7     660      5.41     4.9073      4.5573                                                                          */
/*    8     800      5.80     6.2655      6.1743                                                                          */
/*    9     880      6.70     7.4231      8.0212                                                                          */
/*   10     925      9.05     8.5700      9.4663                                                                          */
/*   11     950     11.70    11.7000     10.4151                                                                          */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*____         _       _                                          _   __ _ _   _           _
|___  |  _ __ | | ___ | |_   _ __ __ ___      __   __ _ _ __   __| | / _(_) |_| |_ ___  __| |
   / /  | `_ \| |/ _ \| __| | `__/ _` \ \ /\ / /  / _` | `_ \ / _` || |_| | __| __/ _ \/ _` |
  / /   | |_) | | (_) | |_  | | | (_| |\ V  V /  | (_| | | | | (_| ||  _| | |_| ||  __/ (_| |
 /_/    | .__/|_|\___/ \__| |_|  \__,_| \_/\_/    \__,_|_| |_|\__,_||_| |_|\__|\__\___|\__,_|
        |_|
*/

options ls=64 ps=64;
proc plot data=allTre(
 rename=y_org=
 y12345678901234567890123456789);
 plot
  y12345678901234567890123456789
    *x_org='o'
  y_stack*x_org='s'
  y_rog*x_org='r'
 /box overlay vref=0 11.7 ;
 run;quit;

/*___    _     _                                         _
 / _ \  | |__ (_)  _ __ ___  ___   _ __  _ __ __ _ _ __ | |__
| (_) | | `_ \| | | `__/ _ \/ __| | `_ \| `__/ _` | `_ \| `_ \
 \__, | | | | | | | | |  __/\__ \ | |_) | | | (_| | |_) | | | |
   /_/  |_| |_|_| |_|  \___||___/ | .__/|_|  \__,_| .__/|_| |_|
                                  |_|             |_|
*/

proc sql;
  create
     table addorg as
  select
     roger.xs   as xs
    ,roger.ys   as roger
    ,roger.yorg as ys
    ,ops.predicted     as ops_y

  from
     roger left join ops
  on
     roger.xs = ops.xs
;quit;

%utlfkil(d:/png/pltcum.png);
ods listing gpath='d:/png/';
ods graphics / reset=all width=11in height=8.5in imagename="pltcum" imagefmt=png;

title "CDF on Unit Square";
proc sgplot data=addorg;
label  ys="Original"
       ops_y="Stackoverflow"
;
 series x=xs y=roger    / lineattrs=(thickness=2pt color=blue);
 series x=xs y=ys       / lineattrs=(thickness=3pt  color=green);
 series x=xs y=ops_y    / lineattrs=(thickness=1pt color=red);
 inset ( " " ="Stackoverflow overfit "
         " " ="a=A   -10.06925  "
         " " ="b=B 7.465462E-8  "
         " " ="c=C -0.00012226  "
         " " ="d=D     0.07033  "
         " " ="y =a*xs**3 + b*xs**2 +c*xs + d"
 ) / position=bottomright textattrs=( Size=12pt color=red) valuealign=left;

/*
 inset ( " " ="Rogers fit "
         " " ="a =  307.4 "
         " " ="b =  500.4 "
         " " ="c = 0.0422 "
         " " ="Ymapped =CDF('BETA',x,a,b)"
 ) / position=bottomright textattrs=( Size=12pt color=blue) valuealign=left;

*/
 inset ( " " ="Original Data "
 ) / position=right textattrs=( Size=12pt color=green) valuealign=left;

 inset ( " " ="Rogers fit mapped units"
         " " ="a =  0.4938 "
         " " ="b =  0.2961 "
         " " ="Ymapped =CDF('BETA',x,a,b)"
         " " ="Yorg    =CDF('BETA',(xorg-200)/750,0.4938,0.2961)*11.7"
 ) / position=topleft textattrs=( Size=12pt color=blue) valuealign=left;

keylegend / location = inside position=left;

run;quit;
ods graphics off;


/*              _
  ___ _ __   __| |
 / _ \ `_ \ / _` |
|  __/ | | | (_| |
 \___|_| |_|\__,_|

*/
