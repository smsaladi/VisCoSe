int locpenaltyd = -10000;
int locexpenaltyd = -3000;
/*
char locaminod[26] = "GASTPLIMVDNEQFYWKRHCXXX.-U";
*/
char locaminod[] = "ARNDCQEGHILKMFPSTWYVBZX.-U";
char locgrpd[] = 
{
	0, 3, 2, 2, 5, 2, 2, 0, 3, 1, 1, 3, 1, 4, 0, 0, 0, 4, 4, 1, 2, 2,
	6, 6, 6, 6,
};
int locn_disd[26][26] = 
        {
 2400, -600, -300, -300,  500, -200,    0,  500, -800, -800,
-1200, -400, -700,-2300,  300, 1100,  600,-3600,-2200,  100,
    0,    0,    0,    0,    0,    0,
 -600, 4700,  300, -300,-2200, 1500,  400,-1000,  600,-2400,
-2200, 2700,-1700,-3200, -900, -200, -200,-1600,-1800,-2000,
    0,    0,    0,    0,    0,    0,
 -300,  300, 3800, 2200,-1800,  700,  900,  400, 1200,-2800,
-3000,  800,-2200,-3100, -900,  900,  500,-3600,-1400,-2200,
    0,    0,    0,    0,    0,    0,
 -300, -300, 2200, 4700,-3200,  900, 2700,  100,  400,-3800,
-4000,  500,-3000,-4500, -700,  500,    0,-5200,-2800,-2900,
    0,    0,    0,    0,    0,    0,
  500,-2200,-1800,-3200,11500,-2400,-3000,-2000,-1300,-1100,
-1500,-2800, -900, -800,-3100,  100, -500,-1000, -500,    0,
    0,    0,    0,    0,    0,    0,
 -200, 1500,  700,  900,-2400, 2700, 1700,-1000, 1200,-1900,
-1600, 1500,-1000,-2600, -200,  200,    0,-2700,-1700,-1500,
    0,    0,    0,    0,    0,    0,
    0,  400,  900, 2700,-3000, 1700, 3600, -800,  400,-2700,
-2800, 1200,-2000,-3900, -500,  200, -100,-4300,-2700,-1900,
    0,    0,    0,    0,    0,    0,
  500,-1000,  400,  100,-2000,-1000, -800, 6600,-1400,-4500,
-4400,-1100,-3500,-5200,-1600,  400,-1100,-4000,-4000,-3300,
    0,    0,    0,    0,    0,    0,
 -800,  600, 1200,  400,-1300, 1200,  400,-1400, 6000,-2200,
-1900,  600,-1300, -100,-1100, -200, -300, -800, 2200,-2000,
    0,    0,    0,    0,    0,    0,
 -800,-2400,-2800,-3800,-1100,-1900,-2700,-4500,-2200, 4000,
 2800,-2100, 2500, 1000,-2600,-1800, -600,-1800, -700, 3100,
    0,    0,    0,    0,    0,    0,
-1200,-2200,-3000,-4000,-1500,-1600,-2800,-4400,-1900, 2800,
 4000,-2100, 2800, 2000,-2300,-2100,-1300, -700,    0, 1800,
    0,    0,    0,    0,    0,    0,
 -400, 2700,  800,  500,-2800, 1500, 1200,-1100,  600,-2100,
-2100, 3200,-1400,-3300, -600,  100,  100,-3500,-2100,-1700,
    0,    0,    0,    0,    0,    0,
 -700,-1700,-2200,-3000, -900,-1000,-2000,-3500,-1300, 2500,
 2800,-1400, 4300, 1600,-2400,-1400, -600,-1000, -200, 1600,
    0,    0,    0,    0,    0,    0,
-2300,-3200,-3100,-4500, -800,-2600,-3900,-5200, -100, 1000,
 2000,-3300, 1600, 7000,-3800,-2800,-2200, 3600, 5100,  100,
    0,    0,    0,    0,    0,    0,
  300, -900, -900, -700,-3100, -200, -500,-1600,-1100,-2600,
-2300, -600,-2400,-3800, 7600,  400,  100,-5000,-3100,-1800,
    0,    0,    0,    0,    0,    0,
 1100, -200,  900,  500,  100,  200,  200,  400, -200,-1800,
-2100,  100,-1400,-2800,  400, 2200, 1500,-3300,-1900,-1000,
    0,    0,    0,    0,    0,    0,
  600, -200,  500,    0, -500,    0, -100,-1100, -300, -600,
-1300,  100, -600,-2200,  100, 1500, 2500,-3500,-1900,    0,
    0,    0,    0,    0,    0,    0,
-3600,-1600,-3600,-5200,-1000,-2700,-4300,-4000, -800,-1800,
 -700,-3500,-1000, 3600,-5000,-3300,-3500,14200, 4100,-2600,
    0,    0,    0,    0,    0,    0,
-2200,-1800,-1400,-2800, -500,-1700,-2700,-4000, 2200, -700,
    0,-2100, -200, 5100,-3100,-1900,-1900, 4100, 7800,-1100,
    0,    0,    0,    0,    0,    0,
  100,-2000,-2200,-2900,    0,-1500,-1900,-3300,-2000, 3100,
 1800,-1700, 1600,  100,-1800,-1000,    0,-2600,-1100, 3400,
    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,
#if 0
   280,     1,    83,    94,    -6,    57,   112,   169,    -7,    52,
    41,    46,   -10,   -52,   142,   177,   169,  -210,   -57,   129,
     0,     0,     0,     0,     0,  -300,

     1,   314,    16,   -57,   -90,    90,   -35,   -49,    88,   -26,
   -46,   190,   -22,  -120,    45,    69,     5,    30,  -140,   -39,
     0,     0,     0,     0,     0,  -300,

    83,    16,   236,   168,  -118,    52,   103,    92,   122,   -22,
   -27,   130,   -99,   -91,    12,   135,    95,  -158,   -17,   -22,
     0,     0,     0,     0,     0,  -300,

    94,   -57,   168,   278,  -185,    95,   215,    99,    49,   -52,
   -97,    72,  -139,  -193,   -16,    82,    48,  -283,  -135,   -29,
     0,     0,     0,     0,     0,  -300,

    -6,   -90,  -118,  -185,   368,  -198,  -186,   -68,   -94,   -51,
  -182,  -163,  -237,  -165,   -69,    69,   -44,  -319,    11,   -12,
     0,     0,     0,     0,     0,  -300,

    57,    90,    52,    95,  -198,   280,   168,     4,   150,   -60,
    26,    94,   -32,  -164,    63,    30,     9,  -226,  -149,   -18,
     0,     0,     0,     0,     0,  -300,

   112,   -35,   103,   215,  -186,   168,   282,    80,    37,   -22,
   -50,    69,   -87,  -180,    21,    65,    30,  -305,  -111,     0,
     0,     0,     0,     0,     0,  -300,

   169,   -49,    92,    99,   -68,     4,    80,   335,   -42,   -59,
   -53,    28,   -91,   -83,    43,   153,    70,  -252,  -154,    49,
     0,     0,     0,     0,     0,  -300,

    -7,    88,   122,    49,   -94,   150,    37,   -42,   312,   -90,
    -4,    29,  -140,   -35,    29,    14,   -23,  -135,    21,   -29,
     0,     0,     0,     0,     0,  -300,

    52,   -26,   -22,   -52,   -51,   -60,   -22,   -59,   -90,   280,
   144,     0,    65,    70,   -57,    -3,    76,  -248,   -38,   210,
     0,     0,     0,     0,     0,  -300,

    41,   -46,   -27,   -97,  -182,    26,   -50,   -53,    -4,   144,
   345,    -3,   150,   115,    -6,   -12,    21,   -57,    -8,   148,
     0,     0,     0,     0,     0,  -300,

    46,   190,   130,    72,  -163,    94,    69,    28,    29,     0,
    -3,   327,    66,  -149,    19,    94,    99,  -153,   -98,   -19,
     0,     0,     0,     0,     0,  -300,

   -10,   -22,   -99,  -139,  -237,   -32,   -87,   -91,  -140,    65,
   150,    66,   275,   -18,   -97,   -37,    -5,  -264,  -166,    77,
     0,     0,     0,     0,     0,  -300,

   -52,  -120,   -91,  -193,  -165,  -164,  -180,   -83,   -35,    70,
   115,  -149,   -18,   344,  -117,   -29,   -66,   -42,   200,   -15,
     0,     0,     0,     0,     0,  -300,

   142,    45,    12,   -16,   -69,    63,    21,    43,    29,   -57,
    -6,    19,   -97,  -117,   324,   125,    69,  -238,  -175,    21,
     0,     0,     0,     0,     0,  -300,

   177,    69,   135,    82,    69,    30,    65,   153,    14,    -3,
   -12,    94,   -37,   -29,   125,   259,   177,   -57,   -40,    41,
     0,     0,     0,     0,     0,  -300,

   169,     5,    95,    48,   -44,     9,    30,    70,   -23,    76,
    21,    99,    -5,   -66,    69,   177,   278,  -208,   -43,   102,
     0,     0,     0,     0,     0,  -300,

  -210,    30,  -158,  -283,  -319,  -226,  -305,  -252,  -135,  -248,
   -57,  -153,  -264,   -42,  -238,   -57,  -208,   370,   -73,  -263,
     0,     0,     0,     0,     0,  -300,

   -57,  -140,   -17,  -135,    11,  -149,  -111,  -154,    21,   -38,
    -8,   -98,  -166,   200,  -175,   -40,   -43,   -73,   344,   -54,
     0,     0,     0,     0,     0,  -300,

   129,   -39,   -22,   -29,   -12,   -18,     0,    49,   -29,   210,
   148,   -19,    77,   -15,    21,    41,   102,  -263,   -54,   309,
     0,     0,     0,     0,     0,  -300,

    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0, -300,

    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0, -300,

    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0, -300,

    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0, -300,

    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,
    0,    0,    0,    0,    0,    0,

 -300, -300, -300, -300, -300, -300, -300, -300, -300, -300,
 -300, -300, -300, -300, -300, -300, -300, -300, -300, -300,
 -300, -300, -300, -300,    0,  700,
#endif
     };

