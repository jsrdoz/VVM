	  V     k820309    h          14.0        öÞBV                                                                                                           
       ../../DATA/t01/rad_driver.f RAD_DRIVER              RAD_DRIVER_RRTM INITIALIZE_RADIATION ISINITIALIZED_RADDRIVER TRACESINI P_FACTOR_XY P_COSZRS_XY                      @                              
       SHR_ORB_PARAMS                      @                              
       COMPUTERE_LIQUID COMPUTERE_ICE ALBEDO                                                     
       KIND_RB #         @                                                   #SHR_ORB_PARAMS%SIN    #SHR_ORB_PARAMS%COS    #SHR_ORB_PARAMS%FLOAT    #SHR_ORB_PARAMS%ABS    #SHR_ORB_PARAMS%SQRT 	   #SHR_ORB_PARAMS%ATAN 
   #IYEAR_AD    #ECCEN    #OBLIQ    #MVELP    #OBLIQR    #LAMBM0    #MVELPP    #LOG_PRINT                  @                                 SIN               @                                 COS               @                                 FLOAT               @                                 ABS               @                            	     SQRT               @                            
     ATAN           
   @                                                  
  @                                  	                 
  @                                  	                 
  @                                  	                   @                                  	                   @                                  	                   @                                  	                 
   @                                         %         H                                                 	       #COMPUTERE_LIQUID%MIN    #COMPUTERE_LIQUID%MAX    #COMPUTERE_LIQUID%PRESENT    #TEMPERATURE    #LANDFRAC    #ICEFRAC    #SNOWH                  @                                 MIN               @                                 MAX               @                                 PRESENT           
   @                                   	                
   @                                   	                
  @                                   	                
  @                                   	      %         H                                                 	       #COMPUTERE_ICE%MIN    #COMPUTERE_ICE%MAX    #COMPUTERE_ICE%INT    #TEMPERATURE                  @                                 MIN               @                                 MAX               @                                 INT           
   @                                   	      #         @                                                     #OCEAN !   #COSZRS "   #ASDIR #   #ALDIR $   #ASDIF %   #ALDIF &             
   @                              !                     
   @                              "                   	              &                                                       @                              #                   	               &                                                       @                              $                   	               &                                                       @                              %                   	               &                                                       @                              &                   	               &                                                          @                           '                          #MPI_BOTTOM (                @                           (                                   @                           )                          #MPI_IN_PLACE *                @                           *                                   @                           +                          #MPI_ARGV_NULL ,   -             @                           ,                                 p          p            p                                                 @                           -                          #MPI_ARGVS_NULL .   -             @                           .                                 p          p            p                                                 @                           /                          #MPI_ERRCODES_IGNORE 0                @                           0                                 p          p            p                                                 @                           1                          #MPI_STATUS_IGNORE 2                @                           2                                 p          p            p                                                 @                           3                          #MPI_STATUSES_IGNORE 4                @                           4                                 p          p          p            p          p                                                 @                           5                          #MPI_UNWEIGHTED 6                @                           6                                   @                           7                          #MPI_WEIGHTS_EMPTY 8                @                           8                    #         @                                   9                5   #RAD_DRIVER_RRTM%TINY :   #RAD_DRIVER_RRTM%FLOOR ;   #RAD_DRIVER_RRTM%FLOAT <   #RAD_DRIVER_RRTM%SUM =   #RAD_DRIVER_RRTM%COS >   #RAD_DRIVER_RRTM%MAX ?   #RAD_DRIVER_RRTM%MERGE @   #RAD_DRIVER_RRTM%MIN A   #RAD_DRIVER_RRTM%ALL B   #RAD_DRIVER_RRTM%REAL C   #RAD_DRIVER_RRTM%SPREAD D   #NX E   #NZM F   #LAT G   #PRES H   #PRESI I   #TABS J   #QV K   #QCL L   #QCI M   #TG N   #O3 O   #CO2 P   #CH4 Q   #N2O R   #O2 S   #CFC11 T   #CFC12 U   #CFC22 V   #CCL4 W   #DOLONGWAVE X   #DOSHORTWAVE Y   #DOPERPETUAL Z   #DOSEASONS [   #DOSOLARCONSTANT \   #SOLAR_CONSTANT ]   #ZENITH_ANGLE ^   #DAY _   #DAY0 `   #LATITUDE a   #LONGITUDE b   #P_FACTOR_SLICE c   #P_COSZRS_SLICE d   #OCEAN e   #GGR f   #CP g   #MASTERPROC h   #LWUP i   #LWDOWN j   #LWUPCLEARSKY k   #LWDOWNCLEARSKY l   #SWUP m   #SWDOWN n   #SWUPCLEARSKY o   #SWDOWNCLEARSKY p   #SWHEATINGRATE q   #SWHEATINGRATECLEARSKY r   #LWHEATINGRATE s   #LWHEATINGRATECLEARSKY t   #COSZRS u   #LWP v   #IWP w   #LIQUIDRE x   #ICERE y                                                                                                                                  @                            :     TINY               @                            ;     FLOOR               @                            <     FLOAT               @                            =     SUM               @                            >     COS               @                            ?     MAX               @                            @     MERGE               @                            A     MIN               @                            B     ALL               @             @              C     REAL               @                            D     SPREAD           
  @                               E                     
                                  F                     
                                  G                    
                                  H                    	    p          5  p        r F       5  p        r F                              
                                  I                    	    p           5  p        r F   n                                       1     5  p        r F   n                                      1                                    
                                  J                    	      p        5  p        r E   p          5  p        r E     5  p        r F       5  p        r E     5  p        r F                              
                                  K                    	      p        5  p        r E   p          5  p        r E     5  p        r F       5  p        r E     5  p        r F                              
                                  L                    	      p        5  p        r E   p          5  p        r E     5  p        r F       5  p        r E     5  p        r F                              
                                  M                    	      p        5  p        r E   p          5  p        r E     5  p        r F       5  p        r E     5  p        r F                              
                                  N                    	    p          5  p        r E       5  p        r E                              
                                  O                    	    p          5  p        r F       5  p        r F                              
                                  P                    	    p          5  p        r F       5  p        r F                              
                                  Q                    	    p          5  p        r F       5  p        r F                              
                                  R                    	    p          5  p        r F       5  p        r F                              
                                  S                    	     p          5  p        r F       5  p        r F                              
                                  T                    	 !   p          5  p        r F       5  p        r F                              
                                  U                    	 "   p          5  p        r F       5  p        r F                              
                                  V                    	 #   p          5  p        r F       5  p        r F                              
                                  W                    	 $   p          5  p        r F       5  p        r F                               
                                  X                     
                                  Y                     
                                  Z                     
                                  [                     
                                  \                     
                                  ]     	                
                                  ^     	                
  @                               _     	                
  @                               `     	               
                                  a                    	    p          5  p        r E       5  p        r E                              
                                  b                    	    p          5  p        r E       5  p        r E                              
                                  c                    	    p          5  p        r E       5  p        r E                              
                                  d                    	    p          5  p        r E       5  p        r E                               
  @                               e                     
                                  f     	                
                                  g     	                
                                  h                    D @                              i                    
 %      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       2    5  p        r E      5  p        r F   n                                      2                                    D @                              j                    
 &      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       2    5  p        r E      5  p        r F   n                                      2                                    D @                              k                    
 '      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       2    5  p        r E      5  p        r F   n                                      2                                    D @                              l                    
 (      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       2    5  p        r E      5  p        r F   n                                      2                                    D @                              m                    
 )      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       2    5  p        r E      5  p        r F   n                                      2                                    D @                              n                    
 *      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       2    5  p        r E      5  p        r F   n                                      2                                    D @                              o                    
 +      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       2    5  p        r E      5  p        r F   n                                      2                                    D @                              p                    
 ,      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       2    5  p        r E      5  p        r F   n                                      2                                    D @                              q                    
 -      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       1    5  p        r E      5  p        r F   n                                      1                                    D @                              r                    
 .      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       1    5  p        r E      5  p        r F   n                                      1                                    D @                              s                    
 /      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       1    5  p        r E      5  p        r F   n                                      1                                    D @                              t                    
 0      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       1    5  p        r E      5  p        r F   n                                      1                                     D                                 u     	                D @                              v                    
 1      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       1    5  p        r E      5  p        r F   n                                      1                                    D @                              w                    
 2      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       1    5  p        r E      5  p        r F   n                                      1                                    D @                              x                    
 3      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       1    5  p        r E      5  p        r F   n                                      1                                    D @                              y                    
 4      p        5  p        r E   p          5  p        r E      5  p        r F   n                                       1    5  p        r E      5  p        r F   n                                      1                           #         @                                   z                    #NX {   #NY |   #CP }   #IYEAR ~   #DAY0    #LATITUDE    #LONGITUDE    #DOPERPETUAL                                                
                                  {                     
                                  |                     
                                  }     	                
  @                               ~                     
  @                                    	                
  @                                                  	 R             &                   &                                                     
  @                                                  	 S             &                   &                                                     
                                                       @                                             #         @                                                      #TRACESINI%ABS    #TRACESINI%MAXVAL    #TRACESINI%TRIM    #TRACESINI%MAX    #TRACESINI%MIN    #NZM    #PRES    #PRESI    #GGR    #MASTERPROC                                                                                                                                                                                                                                                                                                                           @                                 ABS               @                                 MAXVAL               @                                 TRIM               @                                 MAX               @                                 MIN           
                                                      
                                                      	 T   p          5  p        r        5  p        r                               
                                                      	 U   p           5  p        r    n                                       1     5  p        r    n                                      1                                     
                                       	                
                                                      @                                                   	                &                   &                                                    @                                                   	                &                   &                                                  /      fn#fn     Ï   o   b   uapp(RAD_DRIVER    >  O   J  SHR_ORB_MOD *     f   J  CAM_RAD_PARAMETERIZATIONS    ó  H   J  PARKIND +   ;  >      SHR_ORB_PARAMS+SHR_ORB_MOD 3   y  <      SHR_ORB_PARAMS%SIN+SHR_ORB_MOD=SIN 3   µ  <      SHR_ORB_PARAMS%COS+SHR_ORB_MOD=COS 7   ñ  >      SHR_ORB_PARAMS%FLOAT+SHR_ORB_MOD=FLOAT 3   /  <      SHR_ORB_PARAMS%ABS+SHR_ORB_MOD=ABS 5   k  =      SHR_ORB_PARAMS%SQRT+SHR_ORB_MOD=SQRT 5   ¨  =      SHR_ORB_PARAMS%ATAN+SHR_ORB_MOD=ATAN 4   å  @   e   SHR_ORB_PARAMS%IYEAR_AD+SHR_ORB_MOD 1   %  @   e   SHR_ORB_PARAMS%ECCEN+SHR_ORB_MOD 1   e  @   e   SHR_ORB_PARAMS%OBLIQ+SHR_ORB_MOD 1   ¥  @   e   SHR_ORB_PARAMS%MVELP+SHR_ORB_MOD 2   å  @   e   SHR_ORB_PARAMS%OBLIQR+SHR_ORB_MOD 2   %  @   e   SHR_ORB_PARAMS%LAMBM0+SHR_ORB_MOD 2   e  @   e   SHR_ORB_PARAMS%MVELPP+SHR_ORB_MOD 5   ¥  @   e   SHR_ORB_PARAMS%LOG_PRINT+SHR_ORB_MOD ;   å  Ù       COMPUTERE_LIQUID+CAM_RAD_PARAMETERIZATIONS C   ¾  <      COMPUTERE_LIQUID%MIN+CAM_RAD_PARAMETERIZATIONS=MIN C   ú  <      COMPUTERE_LIQUID%MAX+CAM_RAD_PARAMETERIZATIONS=MAX K   6  @      COMPUTERE_LIQUID%PRESENT+CAM_RAD_PARAMETERIZATIONS=PRESENT G   v  @   e   COMPUTERE_LIQUID%TEMPERATURE+CAM_RAD_PARAMETERIZATIONS D   ¶  @   e   COMPUTERE_LIQUID%LANDFRAC+CAM_RAD_PARAMETERIZATIONS C   ö  @   e   COMPUTERE_LIQUID%ICEFRAC+CAM_RAD_PARAMETERIZATIONS A   6	  @   e   COMPUTERE_LIQUID%SNOWH+CAM_RAD_PARAMETERIZATIONS 8   v	  ¦       COMPUTERE_ICE+CAM_RAD_PARAMETERIZATIONS @   
  <      COMPUTERE_ICE%MIN+CAM_RAD_PARAMETERIZATIONS=MIN @   X
  <      COMPUTERE_ICE%MAX+CAM_RAD_PARAMETERIZATIONS=MAX @   
  <      COMPUTERE_ICE%INT+CAM_RAD_PARAMETERIZATIONS=INT D   Ð
  @   e   COMPUTERE_ICE%TEMPERATURE+CAM_RAD_PARAMETERIZATIONS 1            ALBEDO+CAM_RAD_PARAMETERIZATIONS 7     @   e   ALBEDO%OCEAN+CAM_RAD_PARAMETERIZATIONS 8   Û     e   ALBEDO%COSZRS+CAM_RAD_PARAMETERIZATIONS 7   g     e   ALBEDO%ASDIR+CAM_RAD_PARAMETERIZATIONS 7   ó     e   ALBEDO%ALDIR+CAM_RAD_PARAMETERIZATIONS 7        e   ALBEDO%ASDIF+CAM_RAD_PARAMETERIZATIONS 7        e   ALBEDO%ALDIF+CAM_RAD_PARAMETERIZATIONS .     `      RAD_DRIVER!MPI_FORTRAN_BOTTOM    ÷  H      MPI_BOTTOM 0   ?  b      RAD_DRIVER!MPI_FORTRAN_IN_PLACE    ¡  H      MPI_IN_PLACE 1   é  c      RAD_DRIVER!MPI_FORTRAN_ARGV_NULL    L  ¤      MPI_ARGV_NULL 2   ð  d      RAD_DRIVER!MPI_FORTRAN_ARGVS_NULL    T  ¤      MPI_ARGVS_NULL 7   ø  i      RAD_DRIVER!MPI_FORTRAN_ERRCODES_IGNORE $   a  ¤      MPI_ERRCODES_IGNORE 5     g      RAD_DRIVER!MPI_FORTRAN_STATUS_IGNORE "   l  ¤      MPI_STATUS_IGNORE 7     i      RAD_DRIVER!MPI_FORTRAN_STATUSES_IGNORE $   y  Ä      MPI_STATUSES_IGNORE 2   =  d      RAD_DRIVER!MPI_FORTRAN_UNWEIGHTED    ¡  H      MPI_UNWEIGHTED 5   é  g      RAD_DRIVER!MPI_FORTRAN_WEIGHTS_EMPTY "   P  H      MPI_WEIGHTS_EMPTY             RAD_DRIVER_RRTM %   +  =      RAD_DRIVER_RRTM%TINY &   h  >      RAD_DRIVER_RRTM%FLOOR &   ¦  >      RAD_DRIVER_RRTM%FLOAT $   ä  <      RAD_DRIVER_RRTM%SUM $      <      RAD_DRIVER_RRTM%COS $   \  <      RAD_DRIVER_RRTM%MAX &     >      RAD_DRIVER_RRTM%MERGE $   Ö  <      RAD_DRIVER_RRTM%MIN $     <      RAD_DRIVER_RRTM%ALL %   N  =      RAD_DRIVER_RRTM%REAL '     ?      RAD_DRIVER_RRTM%SPREAD #   Ê  @   a   RAD_DRIVER_RRTM%NX $   
  @   a   RAD_DRIVER_RRTM%NZM $   J  @   a   RAD_DRIVER_RRTM%LAT %     ´   a   RAD_DRIVER_RRTM%PRES &   >  &  a   RAD_DRIVER_RRTM%PRESI %   d   $  a   RAD_DRIVER_RRTM%TABS #   !  $  a   RAD_DRIVER_RRTM%QV $   ¬"  $  a   RAD_DRIVER_RRTM%QCL $   Ð#  $  a   RAD_DRIVER_RRTM%QCI #   ô$  ´   a   RAD_DRIVER_RRTM%TG #   ¨%  ´   a   RAD_DRIVER_RRTM%O3 $   \&  ´   a   RAD_DRIVER_RRTM%CO2 $   '  ´   a   RAD_DRIVER_RRTM%CH4 $   Ä'  ´   a   RAD_DRIVER_RRTM%N2O #   x(  ´   a   RAD_DRIVER_RRTM%O2 &   ,)  ´   a   RAD_DRIVER_RRTM%CFC11 &   à)  ´   a   RAD_DRIVER_RRTM%CFC12 &   *  ´   a   RAD_DRIVER_RRTM%CFC22 %   H+  ´   a   RAD_DRIVER_RRTM%CCL4 +   ü+  @   a   RAD_DRIVER_RRTM%DOLONGWAVE ,   <,  @   a   RAD_DRIVER_RRTM%DOSHORTWAVE ,   |,  @   a   RAD_DRIVER_RRTM%DOPERPETUAL *   ¼,  @   a   RAD_DRIVER_RRTM%DOSEASONS 0   ü,  @   a   RAD_DRIVER_RRTM%DOSOLARCONSTANT /   <-  @   a   RAD_DRIVER_RRTM%SOLAR_CONSTANT -   |-  @   a   RAD_DRIVER_RRTM%ZENITH_ANGLE $   ¼-  @   a   RAD_DRIVER_RRTM%DAY %   ü-  @   a   RAD_DRIVER_RRTM%DAY0 )   <.  ´   a   RAD_DRIVER_RRTM%LATITUDE *   ð.  ´   a   RAD_DRIVER_RRTM%LONGITUDE /   ¤/  ´   a   RAD_DRIVER_RRTM%P_FACTOR_SLICE /   X0  ´   a   RAD_DRIVER_RRTM%P_COSZRS_SLICE &   1  @   a   RAD_DRIVER_RRTM%OCEAN $   L1  @   a   RAD_DRIVER_RRTM%GGR #   1  @   a   RAD_DRIVER_RRTM%CP +   Ì1  @   a   RAD_DRIVER_RRTM%MASTERPROC %   2    a   RAD_DRIVER_RRTM%LWUP '   ¢3    a   RAD_DRIVER_RRTM%LWDOWN -   85    a   RAD_DRIVER_RRTM%LWUPCLEARSKY /   Î6    a   RAD_DRIVER_RRTM%LWDOWNCLEARSKY %   d8    a   RAD_DRIVER_RRTM%SWUP '   ú9    a   RAD_DRIVER_RRTM%SWDOWN -   ;    a   RAD_DRIVER_RRTM%SWUPCLEARSKY /   &=    a   RAD_DRIVER_RRTM%SWDOWNCLEARSKY .   ¼>    a   RAD_DRIVER_RRTM%SWHEATINGRATE 6   R@    a   RAD_DRIVER_RRTM%SWHEATINGRATECLEARSKY .   èA    a   RAD_DRIVER_RRTM%LWHEATINGRATE 6   ~C    a   RAD_DRIVER_RRTM%LWHEATINGRATECLEARSKY '   E  @   a   RAD_DRIVER_RRTM%COSZRS $   TE    a   RAD_DRIVER_RRTM%LWP $   êF    a   RAD_DRIVER_RRTM%IWP )   H    a   RAD_DRIVER_RRTM%LIQUIDRE &   J    a   RAD_DRIVER_RRTM%ICERE %   ¬K  Å       INITIALIZE_RADIATION (   qL  @   a   INITIALIZE_RADIATION%NX (   ±L  @   a   INITIALIZE_RADIATION%NY (   ñL  @   a   INITIALIZE_RADIATION%CP +   1M  @   a   INITIALIZE_RADIATION%IYEAR *   qM  @   a   INITIALIZE_RADIATION%DAY0 .   ±M  ¤   a   INITIALIZE_RADIATION%LATITUDE /   UN  ¤   a   INITIALIZE_RADIATION%LONGITUDE 1   ùN  @   a   INITIALIZE_RADIATION%DOPERPETUAL (   9O  @       ISINITIALIZED_RADDRIVER    yO        TRACESINI    Q  <      TRACESINI%ABS !   ÀQ  ?      TRACESINI%MAXVAL    ÿQ  =      TRACESINI%TRIM    <R  <      TRACESINI%MAX    xR  <      TRACESINI%MIN    ´R  @   a   TRACESINI%NZM    ôR  ´   a   TRACESINI%PRES     ¨S  &  a   TRACESINI%PRESI    ÎT  @   a   TRACESINI%GGR %   U  @   a   TRACESINI%MASTERPROC    NU  ¤       P_FACTOR_XY    òU  ¤       P_COSZRS_XY 