;; Mibcom: Metabolic Individual Based COMmunity model for cental-place forager mammals in heterogenous landscapes
;; copyright Leonna Szangolies, Florian Jeltsch, University of Potsdam

extensions [profiler csv]

globals [
  ;;species definition with mass, shelter need, foraging type and food preference
  s0mass s0massdev s1mass s1massdev  s2mass s2massdev s3mass s3massdev s4mass s4massdev s5mass s5massdev s6mass s6massdev s7mass s7massdev s8mass s8massdev s9mass s9massdev
  s0shelt s1shelt s2shelt s3shelt s4shelt s5shelt s6shelt s7shelt s8shelt s9shelt
  s0food s1food s2food s3food s4food s5food s6food s7food s8food s9food
  s0fortype s1fortype s2fortype s3fortype s4fortype s5fortype s6fortype s7fortype s8fortype s9fortype

  ;;landscape variables
  feed-matrix feed-struct season-var feed-var
  actcov cover_large cover_small
  drought_data

  ;;helper variables for the functions
  day r-search feed shelter-need maxrad minfeed spec-color movecost foodshare

  ;;model parameters
  spec_num mass_order n_immi max_nat_disp growthcost patch-move cover var

  ;; scenario definition
  density_dep small_habitat drought drought_day drought_start

  ;;output variables
  focal_spec                                 ;;mainly for plotting home ranges of only one species
  fail                                       ;;juveniles that did not find a home range (per time step)
  fail3                                      ;;individuals that dy because of no storage
  fail3_spec                                 ;;individuals that dy because of no storage for the different species
  failed-lact                                ;;offspring that gets lost during lactation due to shortage
  mort_order mort_stor mort_food mort_hr     ;;characteristics of individuals before death
  suc_juv                                    ;;successful juveniles
  repro                                      ;;weaned juveniles that can search their own home range
  rep_success_0 rep_success_1 rep_success_2 rep_success_3 rep_success_4 rep_success_5 rep_success_6 rep_success_7 rep_success_8 rep_success_9  ;;total weaned juveniles per female
  rep_success_hr_0 rep_success_hr_1 rep_success_hr_2 rep_success_hr_3 rep_success_hr_4 rep_success_hr_5 rep_success_hr_6 rep_success_hr_7 rep_success_hr_8 rep_success_hr_9  ;;total weaned juveniles per female that found a homerange
  patches0 patches1 patches2 patches3 patches4 patches5 patches6 patches7 patches8 patches9                                                    ;;number of foraging patches
  compet0 compet1 compet2 compet3 compet4 compet5 compet6 compet7 compet8 compet9                                                              ;;number of competitors
  order0 order1 order2 order3 order4 order5 order6 order7 order8 order9                                                                        ;;position in the foraging order
  vis_before0 vis_before1 vis_before2 vis_before3 vis_before4 vis_before5 vis_before6 vis_before7 vis_before8 vis_before9                      ;;number of individuals visiting a foraging patch before myself
  number number0 number1 number2 number3 number4 number5 number6 number7 number8 number9                                                       ;;number of individuals
  hr hr0 hr1 hr2 hr3 hr4 hr5 hr6 hr7 hr8 hr9                                                                                                   ;;mean home range size
  maxhrs maxhr0 maxhr1 maxhr2 maxhr3 maxhr4 maxhr5 maxhr6 maxhr7 maxhr8 maxhr9                                                                 ;;maximum home range size
  countmaxhr countmaxhr0 countmaxhr1 countmaxhr2 countmaxhr3 countmaxhr4 countmaxhr5 countmaxhr6 countmaxhr7 countmaxhr8 countmaxhr9           ;;how often maximum home range size is reached
  ages age0 age1 age2 age3 age4 age5 age6 age7 age8 age9                                                                                       ;;mean home range size
  stors stor0 stor1 stor2 stor3 stor4 stor5 stor6 stor7 stor8 stor9                                                                            ;;mean storage
  abs_stor0 abs_stor1 abs_stor2 abs_stor3 abs_stor4 abs_stor5 abs_stor6 abs_stor7 abs_stor8 abs_stor9                                          ;;mean absolute storage of non-pregnant
  rm rm0 rm1 rm2 rm3 rm4 rm5 rm6 rm7 rm8 rm9                                                                                                   ;;mean real mass
  hrpreg hrpreg0 hrpreg1 hrpreg2 hrpreg3 hrpreg4 hrpreg5 hrpreg6 hrpreg7 hrpreg8 hrpreg9                                                       ;;mean home range size during pregnancy and lactation
  storpreg storpreg0 storpreg1 storpreg2 storpreg3 storpreg4 storpreg5 storpreg6 storpreg7 storpreg8 storpreg9                                 ;;mean storage during pregnancy and lactation
  move0 move1 move2 move3 move4 move5 move6 move7 move8 move9                                                                                  ;;mean movement distance per day
  fmr0 fmr1 fmr2 fmr3 fmr4 fmr5 fmr6 fmr7 fmr8 fmr9                                                                                            ;;total field metabolic rate
  loco0 loco1 loco2 loco3 loco4 loco5 loco6 loco7 loco8 loco9                                                                                  ;;locomotion costs
  repro0 repro1 repro2 repro3 repro4 repro5 repro6 repro7 repro8 repro9                                                                        ;;reproduction costs
  grow0 grow1 grow2 grow3 grow4 grow5 grow6 grow7 grow8 grow9                                                                                  ;;growth costs
  basal0 basal1 basal2 basal3 basal4 basal5 basal6 basal7 basal8 basal9                                                                        ;;basal maintenance costs
  digest0 digest1 digest2 digest3 digest4 digest5 digest6 digest7 digest8 digest9                                                              ;;dirgestion costs
  prod0 prod1 prod2 prod3 prod4 prod5 prod6 prod7 prod8 prod9                                                                                  ;;investment in production: growth and reproduction
  in0 in1 in2 in3 in4 in5 in6 in7 in8 in9                                                                                                      ;;daily intake
  balance0 balance1 balance2 balance3 balance4 balance5 balance6 balance7 balance8 balance9                                                    ;;energy balance: intake versus field metabolic rate
  prodbalance0 prodbalance1 prodbalance2 prodbalance3 prodbalance4 prodbalance5 prodbalance6 prodbalance7 prodbalance8 prodbalance9 prod_sum   ;;production balance: intake versus energy for production
  inpatch0 inpatch1 inpatch2 inpatch3 inpatch4 inpatch5 inpatch6 inpatch7 inpatch8 inpatch9                                                    ;;intake per patch
  inhr0 inhr1 inhr2 inhr3 inhr4 inhr5 inhr6 inhr7 inhr8 inhr9                                                                                  ;;intake relative to home range size

  ;; drought output: characteristics during the drought period
  drought_hr0 drought_hr1 drought_hr2 drought_hr3 drought_hr4 drought_hr5 drought_hr6 drought_hr7 drought_hr8 drought_hr9
  drought_fmr0 drought_fmr1 drought_fmr2 drought_fmr3 drought_fmr4 drought_fmr5 drought_fmr6 drought_fmr7 drought_fmr8 drought_fmr9
  drought_in0 drought_in1 drought_in2 drought_in3 drought_in4 drought_in5 drought_in6 drought_in7 drought_in8 drought_in9
  drought_stor0 drought_stor1 drought_stor2 drought_stor3 drought_stor4 drought_stor5 drought_stor6 drought_stor7 drought_stor8 drought_stor9
]

patches-own [habitat patchfeed patchquali spec-id spec-list eaten save visited patchgrowth drought_patches]
turtles-own [mass species shelter food-pref for-type maxhr hrsize lococost feedrate age average_age preg
             pregcount age_first_repro gest_period lact_period sex order stime hunger core forage_large forage_small storage max_storage
             young mass-neonate mass-weaning mass-embryo real-mass real-young focal share competition forage-patches spec_ind_list
             fmr younggrowth daymoved fmr_basal fmr_growth fmr_repro fmr_loco juvs juvshr preg_next fmr_digest postural postural-time feedspeed vis_before move-factor-allometric dom just-immi
             production intake balance balanceprod]

to setup
  clear-all
  reset-ticks
  set-parameters                                       ;;setup parameters

  ifelse fragmentation_mode = "frag"                   ;;implemented landscape options: "frag" for a fragmentation algorithm and "slass" for single large and several small habitat patches
  [
  setup-patches-algo                                   ;;"frag": based on simple clumping algorithm
  ]
  [
  set cover_large cover * (1 - perc_small) set cover_small (cover * perc_small)      ;;"slass": place specific percentages of large and small habitat patches
  setup-patches-large
  setup-patches-small
  ]

  setup-turtles                                         ;;setup individuals
  find-hr                                               ;;let individuals search an initial home-range

  ask one-of turtles [set focal 1]                      ;;focal species for visualization

  init-out                                              ;;initialize output
  if drought_scenario = True [setup-drought]            ;;initialize drought scenario

end

;---------------------------------------------

to set-parameters                                      ;; characterize species und basic parameters
  let addmass 0
  ;set maxmass 0.1                                      ;; upper limit of mean body mass, now in interface
  set s0mass 0.1 * maxmass + addmass                             ;; mean body mass of species 0 in kg - normal distributed
  set s0massdev 0.2 * s0mass                           ;; std dev of body mass species 0 - normal distribution
  set s1mass 0.2 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s1massdev 0.2 * s1mass                           ;; std dev of body mass species 0 - normal distribution
  set s2mass 0.3 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s2massdev 0.2 * s2mass                           ;; std dev of body mass species 0 - normal distribution
  set s3mass 0.4 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s3massdev 0.2 * s3mass                           ;; std dev of body mass species 0 - normal distribution
  set s4mass 0.5 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s4massdev 0.2 * s4mass                           ;; std dev of body mass species 0 - normal distribution
  set s5mass 0.6 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s5massdev 0.2 * s5mass                           ;; std dev of body mass species 0 - normal distribution
  set s6mass 0.7 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s6massdev 0.2 * s6mass                           ;; std dev of body mass species 0 - normal distribution
  set s7mass 0.8 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s7massdev 0.2 * s7mass                           ;; std dev of body mass species 0 - normal distribution
  set s8mass 0.9 * maxmass + addmass                              ;; mean body mass of species 0 - normal distributed
  set s8massdev 0.2 * s8mass                           ;; std dev of body mass species 0 - normal distribution
  set s9mass maxmass + addmass                                    ;; mean body mass of species 0 - normal distributed
  set s9massdev 0.2 * s9mass                           ;; std dev of body mass species 0 - normal distribution

  set s0shelt 1;0.5                                      ;; index 0..1 of shelter use/need ; currently not  distinguished
  set s1shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s2shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s3shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s4shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s5shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s6shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s7shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s8shelt 1;0.5                                      ;; index 0..1 of shelter use/need
  set s9shelt 1;0.5                                      ;; index 0..1 of shelter use/need

  set s0food "omni"                                    ;; food preference; currently not distinguished
  set s1food "omni"                                    ;; food preference
  set s2food "omni"                                    ;; food preference
  set s3food "omni"                                    ;; food preference
  set s4food "omni"                                    ;; food preference
  set s5food "omni"                                    ;; food preference
  set s6food "omni"                                    ;; food preference
  set s7food "omni"                                    ;; food preference
  set s8food "omni"                                    ;; food preference
  set s9food "omni"                                    ;; food preference

  set s0fortype "central"                              ;; type of forage movement; currently not distinguished
  set s1fortype "central"                              ;; type of forage movement
  set s2fortype "central"                              ;; type of forage movement
  set s3fortype "central"                              ;; type of forage movement
  set s4fortype "central"                              ;; type of forage movement
  set s5fortype "central"                              ;; type of forage movement
  set s6fortype "central"                              ;; type of forage movement
  set s7fortype "central"                              ;; type of forage movement
  set s8fortype "central"                              ;; type of forage movement
  set s9fortype "central"                              ;; type of forage movement

  set actcov cover                                     ;; actual cover after changes
  set fail 0                                           ;; counting failing attempts by offspring to establish homerange
  set fail3 0                                          ;; counting death events caused by food shortage in attempts to re-establish homerange next day
  set failed-lact 0                                    ;; counting death of offspring during gestation and lactation
  set repro 0                                          ;; counting overall annual reproduction

;;parameters to be changed
  set cover cover_percentage * 101 * 101                ;; xx % of shrub cover * nr of patches; xx * 100 * 100
  set feed-matrix 0.0                                   ;; no food prod in matrix
  set feed-struct 1.7127 * feed-amount                  ;; 10% of 68.5g dry mass/ grid cell*day Buchmann et al. 2011
  set season-var 1.7127 * season-variability            ;; seasonal variability of the mean resource level
  set feed-var 0.7                                      ;; standard deviation of daily resource level
  ;set bold_prob 0.5                                    ;; percentage of bold individuals, now in the interface
  ;set mortbold 0.333 / 10000                           ;; increased mortality of bold individuals, now in the interface
  set mass_order prio_fact                              ;; percentage of individuals that are sorted for foraging depending on mass and possibly conspecifics density
  set density_dep false                                 ;; if there is density dependent mortality
  ;set hr_try_juv 10                                    ;; number of attempts that juveniles have to find a home range, now in the interface
  set focal_spec 1                                      ;; focal species of which daily homeranges are shown (cumulative)
  set patch-move 0                                      ;; extra movement within one patch
  ;set move-factor 3                                    ;; factor accounting for non-straight movement in cells, now in the interface
  set small_habitat True                                ;; whether small habitat patches can be home range core cells in the "slass" scneario
  set spec_num 10                                       ;; number of species
  ;set digestion 0.5                                    ;; costs for digestion, now in the interface
  set growthcost (7 + 6) * 1000 / 10                    ;; costs for synthezising 1g of flesh
  ;;set storage_addition 3                              ;; factor adapting maximum storage capacity, now in the interface
  set drought false                                     ;; there is no drought in the beginning of the simulation
  set drought_day 0                                     ;; there is no drought in the beginning of the simulation
  ifelse specs-included = "all" [set drought_start 1095 + 450] [set drought_start 450]  ;; a drought starts in summer time of year 5 in community simulations with all species and in year 2 in single-species simulations
end
;---------------------------------------------

to setup-patches-algo                                  ;; creates a clumped landscape with two types of patches
  let cov 0
  ask patches   [
      set habitat 0
      set patchquali 0
      set spec-id -1
      set spec-list [ ]
      set eaten 0
      set visited 0
   ]
  while [cov < cover]
  [
   ask patches
   [
    ifelse ((habitat > 0) and (any? neighbors with [habitat = 0]) and (cov < cover))
    [ask one-of neighbors with [habitat = 0] [set habitat 1 set cov cov + 1]]
    [if random-float 1 < (1 - clump) and habitat = 0 and (cov < cover) [set habitat 1 set cov cov + 1]]
   ]
  ]
  ask patches[set save (sum [habitat] of neighbors > 6)]

  resource-update
end
;..............................................

to setup-patches-large                                  ;; creates large habitat patches in the landscape
  let cov 0
  ask patches   [
      set patchquali 0
      set habitat 0
      set spec-id -1
      set spec-list [ ]
      set eaten 0
      set visited 0
   ]
  ask n-of 4 patches with [habitat = 0] [set habitat 1 set patchquali 1 set cov cov + 1]

  while [cov < cover_large]
  [
   ask patches with [habitat = 1]
   [
    if ((any? neighbors with [habitat = 0]) and (cov < cover_large))
      [ask one-of neighbors with [habitat = 0] [set habitat 1 set patchquali 1 set cov cov + 1]]
   ]
  ]
end

to setup-patches-small                                   ;; creates small habitat patches in the landscape

  ask n-of ceiling (cover_small / size_small) patches with [habitat = 0 and (count neighbors with [habitat = 0]) = 8]
  [
    ask (n-of size_small patches in-radius ceiling (size_small / 8 )) with [habitat = 0] [set habitat 1 set patchquali 2]
  ]
  ask patches[set save (sum [habitat] of neighbors > 6)]
  resource-update
end
;..............................................

to setup-drought                                          ;; setup drought scenarios with drought length and magnitude
  if drought_type = "data" [read-file]                    ;; read in file with observe drought occurrence
  if drought_scenario_combi > 0                           ;; either set drought_scenario_combi to 0 to manually define drought magnitude and length or use one of the predefined drought scenarios
  [
    if drought_scenario_combi = 1 [set drought_min 0.02 set drought_length 1]
    if drought_scenario_combi = 2 [set drought_min 0.02 set drought_length 7]
    if drought_scenario_combi = 3 [set drought_min 0.05 set drought_length 7]
    if drought_scenario_combi = 4 [set drought_min 0.05 set drought_length 20]
    if drought_scenario_combi = 5 [set drought_min 0.10 set drought_length 20]
    if drought_scenario_combi = 6 [set drought_min 0.10 set drought_length 55]
    if drought_scenario_combi = 7 [set drought_min 0.20 set drought_length 55]
  ]
end
 ;..............................................


to read-file                                              ;; read in drought file (for observed drought occurrence)
  file-open Drought_file
  set drought_data item Location csv:from-file Drought_file
  file-close
end

;--------------------------------------------------

to resource-update                                        ;;daily update of food resources, random variation normal distributed, without seasonality
  ask patches
   [
      ifelse habitat = 0 [set patchfeed feed-matrix] [set patchfeed (random-normal feed-struct feed-var)
      set spec-list [ ]]
   ]
  if drought_scenario = true
   [
      ifelse drought = true
        [
            ifelse drought_type = "data"
            [
              ask patches [set patchfeed (item (ticks - drought_start) drought_data) * patchfeed]
            ]
            [
          ifelse drought_type = "controlled"
          [
            ifelse drought_day < drought_length
            [
              ask patches [set patchfeed drought_min * patchfeed]
              set drought_day drought_day + 1
            ]
            [
              ifelse drought_day < (drought_length + drought_recover)
              [
                ask patches [set patchfeed (drought_min + (drought_day - drought_length) * ((1 - drought_min) / drought_recover)) * patchfeed]
                set drought_day drought_day + 1
              ]
              [
            set drought false
            set drought_day 0
              ]
            ]
          ]
          [
          ifelse drought_type = "aprubt"
          [
          ifelse drought_day < drought_length
            [
              ask patches [set patchfeed (0 + 1 / drought_length * drought_day) * patchfeed] ;0.5 *
              set drought_day drought_day + 1
            ]
            [
            set drought false
            set drought_day 0
            ]
          ]
          [
            ifelse drought_day < drought_length
            [
              set drought_day drought_day + 1
              ifelse drought_day < drought_length / 2
              [
                ask patches [set patchfeed (item drought_day [0.9999 0.975 0.958 0.9333 0.9 0.1 0.0666 0.042 0.025 0.011]) * patchfeed]
              ]
              [
                ifelse drought_day > (drought_length - (drought_length / 2))
                [
                  ask patches [set patchfeed (item (drought_length - drought_day) [0.9999 0.975 0.958 0.9333 0.9 0.1 0.0666 0.042 0.025 0.011]) * patchfeed]
                ]
                [
                  ask patches [set patchfeed 0]
                ]
              ]
            ]
            [
            set drought false
            set drought_day 0
            ]

          ]
          ]
      ]
      ]
         [
        if ticks = drought_start;
        [
          set drought true
          set drought_day 0
         ]
          ]
  ]
end
;---------------------------------------------

to resource-update-seasonal-sin                                    ;;daily update of food resources, random variation normal distributed, with seasonality
  ask patches
   [
      ifelse habitat = 0 [set patchfeed feed-matrix]
      [

        set patchfeed (random-normal ((sin (remainder ticks 360)) * season-var + feed-struct) feed-var)


        set spec-list [ ]
      ]
   ]
  if drought_scenario = true
   [
      ifelse drought = true
          [
            ifelse drought_type = "data"
            [
              ask patches [set patchfeed (item (ticks - drought_start) drought_data) * patchfeed]
            ]
            [
          ifelse drought_type = "controlled"
          [
            ifelse drought_day < drought_recover
            [
                set drought_day drought_day + 1
                ask patches [set patchfeed (1 - drought_day * ((1 - drought_min) / drought_recover)) * patchfeed]
            ]
            [
              ifelse drought_day < (drought_recover + drought_length)
              [
                ;print "Drought_day: " print drought_day
                ;print "Storage: " print mean [storage] of turtles
                ;print "HR: " print mean [hrsize] of turtles
                ;print "Intake: " print mean [intake] of turtles
                ;print "Preg: " print [preg] of turtles
                set drought_day drought_day + 1
                ask patches [set patchfeed drought_min * patchfeed]
              ]
              [
                ifelse drought_day < (drought_length + (2 * drought_recover))
                [
                  set drought_day drought_day + 1
                  ask patches [set patchfeed (drought_min + (drought_day - drought_length - drought_recover) * ((1 - drought_min) / drought_recover)) * patchfeed]
                ]
                [
                  ifelse drought_day = (drought_length + (2 * drought_recover))
                      [
                        set drought_day drought_day + 1
                      ]
                      [
                        set drought false
                        set drought_day 0
                      ]
                ]
              ]
            ]
          ]
          [
          ifelse drought_type = "aprubt"
          [
          ifelse drought_day < drought_length
            [
              ask patches [set patchfeed (0 + 1 / drought_length * drought_day) * patchfeed] ;0.5 *
              set drought_day drought_day + 1
            ]
            [
            set drought false
            set drought_day 0
            ]
          ]
          [
            ifelse drought_day < drought_length
            [
              set drought_day drought_day + 1
              ifelse drought_day < 10
              [
                ask patches [set patchfeed (item drought_day [0.9999 0.975 0.958 0.9333 0.9 0.1 0.0666 0.042 0.025 0.011]) * patchfeed]
              ]
              [
                ifelse drought_day > (drought_length - 10)
                [
                  ask patches [set patchfeed (item (drought_length - drought_day) [0.9999 0.975 0.958 0.9333 0.9 0.1 0.0666 0.042 0.025 0.011]) * patchfeed]
                ]
                [
                  ask patches [set patchfeed 0]
                ]
              ]
            ]
            [
            set drought false
            set drought_day 0
            ]

          ]
          ]
          ]
      ]
         [
        if ticks = drought_start
        [
          ;print "Storage: " print mean [storage] of turtles
          ;print "Age: " print mean [age] of turtles
          ;print "Number: " print count turtles
          ;print "Mass: " print mean [mass] of turtles
          set drought true
          set drought_day 0
         ]
      ]
  ]
end
;---------------------------------------------


to setup-turtles                                       ;; create xx initial individuals, location at 0,0; identifies species, mass,
  crt starting_indivs
                                                       ;; scenario 10 species, species identity 0 .. 9
  let slope -1.5                                       ;; allometric frequency distribution with mass^-1.5, see Buchmann et al. 2012 Ecography
  let s0m s0mass ^ (slope)
  let s1m s1mass ^ (slope)
  let s2m s2mass ^ (slope)
  let s3m s3mass ^ (slope)
  let s4m s4mass ^ (slope)
  let s5m s5mass ^ (slope)
  let s6m s6mass ^ (slope)
  let s7m s7mass ^ (slope)
  let s8m s8mass ^ (slope)
  let s9m s9mass ^ (slope)
  let mass-sum s0m + s1m + s2m + s3m + s4m + s5m + s6m + s7m + s8m + s9m

  ask turtles
  [
   let zuf random-float 1.0

   ifelse zuf < s0m / mass-sum
   [ set species 0 ]
    [ ifelse zuf < (s0m + s1m) / mass-sum [set species 1 ]
      [ ifelse zuf < (s0m + s1m + s2m) / mass-sum [set species 2 ]
        [ ifelse zuf < (s0m + s1m + s2m + s3m) / mass-sum [set species 3 ]
          [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m) / mass-sum [set species 4 ]
            [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m) / mass-sum [set species 5 ]
              [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m + s6m) / mass-sum [set species 6 ]
                [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m + s6m + s7m) / mass-sum [set species 7 ]
                  [ ifelse zuf < (s0m + s1m + s2m + s3m + s4m + s5m + s6m + s7m + s8m) / mass-sum [set species 8 ]
                    [set species 9] ;;last else
                  ];;ifelse8
                ];;ifelse7
              ];;ifelse6
            ];;ifelse5
          ];; ifelse4
        ] ;; ifelse3
      ] ;; ifelse2
    ] ;;ifelse1
    if specs-included != "all" and specs-included != "two" [set species specs-included]
    specify_turtle
    if starting_indivs = 1 [set sex "female"]
  ] ;;ask turtles
end
;.............................................

to specify_turtle                                               ;;characterize and parameterize individuals
  set hrsize 0
  let shy 0
  let zuf1 random-float 1.0
  set color (species * 10 + 5)
  ifelse (zuf1 < bold_prob) [set shy 0][set shy 1]              ;; shy=0 means bold individuals

    if (species = 0) [set mass random-normal s0mass s0massdev   ;; body mass normal distributed
      set shelter s0shelt * shy
      set food-pref s0food
      set for-type s0fortype]
    if (species = 1) [set mass random-normal s1mass s1massdev
      set shelter s1shelt * shy
      set food-pref s1food
      set for-type s1fortype]
    if (species = 2) [set mass random-normal s2mass s2massdev
      set shelter s2shelt * shy
      set food-pref s2food
      set for-type s2fortype]
    if (species = 3) [set mass random-normal s3mass s3massdev
      set shelter s3shelt * shy
      set food-pref s3food
      set for-type s3fortype]
    if (species = 4) [set mass random-normal s4mass s4massdev
      set shelter s4shelt * shy
      set food-pref s4food
      set for-type s4fortype]
    if (species = 5) [set mass random-normal s5mass s5massdev
      set shelter s5shelt * shy
      set food-pref s5food
      set for-type s5fortype]
    if (species = 6) [set mass random-normal s6mass s6massdev
      set shelter s6shelt * shy
      set food-pref s6food
      set for-type s6fortype]
    if (species = 7) [set mass random-normal s7mass s7massdev
      set shelter s7shelt * shy
      set food-pref s7food
      set for-type s7fortype]
    if (species = 8) [set mass random-normal s8mass s8massdev
      set shelter s8shelt * shy
      set food-pref s8food
      set for-type s8fortype]
    if (species = 9) [set mass random-normal s9mass s9massdev
      set shelter s9shelt * shy
      set food-pref s9food
      set for-type s9fortype]

    if mass < 0 [set mass 0 die]

    set preg 0
    set pregcount 0

    set average_age 1766.53 * (mass ^ 0.21)                        ;; average lifespan: allometric formula of mammals after Hamilton et al 2011 [g, days]
    set age_first_repro 293.17 * (mass ^ 0.27)                     ;; age of first reproduction: allometric fromula after Hamilton et al 2011
    set gest_period 64.14 * (mass ^ 0.24)                          ;; gestation period allometric after Hamilton et al 2011
    set lact_period  57.16 * (mass ^ 0.22)                         ;; lactation period allometric after Hamilton et al 2011
    set hunger 0
    set forage_small 0
    set forage_large 0

    set move-factor-allometric move-factor

    set age random (average_age - 180) + 180                       ;;initially random age between 180 days and the average lifespan
    let zuf random-float 1.0
    ifelse zuf < 0.5 [set sex "male"] [set sex "female"]           ;;random sex

    life-history                                                  ;;calculate life history allometries
    set real-mass mass
    set storage 0                                                 ;;initially no storage
    set mass-embryo 0
    set spec_ind_list []
    set juvs 0
    set juvshr 0

    calc-maxhr                                                   ;;calculate allometries
    calc-lococost
    calc-feedrate real-mass
    calc-energetics

  if age > age_first_repro                                       ;; initialization with some already pregnant individuals
  [
  if sex = "female"
  [
  set preg random 2
  if preg = 1
  [
    set real-young young
    set pregcount random gest_period + lact_period
    repeat min (list gest_period pregcount)
    [
      embryo
    ]
    if pregcount > gest_period
    [
      repeat pregcount - gest_period
      [
        juvenile
      ]
    ]
  ]
    let young-mass real-young * mass-embryo
    calc-feedrate (young-mass + real-mass)
    set max_storage storage_addition * (294.8) * ((real-mass + young-mass) ^ 1.19) + feedrate
  ]
  ]

end
;---------------------------------------------
to calc-energetics
    set max_storage storage_addition * (294.8) * (real-mass ^ 1.19) + feedrate ;; allometric storage after Lindstedt and Boyce 1985 ;;alternative: 1572 * mass + feedrate
    set share (real-mass / 0.001) ^ (-0.25)                           ;; allometric share of available food see Buchmann 2011
end
;----------------------------------------------

to calc-maxhr                                                      ;; calculate max home range after Kelt & Van Vuren 2001 in ha, see Buchmann et al 2011
  let maxhrx1 56.23 * (real-mass ^ 0.91)                           ;; max hr for herbivores and omnivores, larger one is used
  let maxhrx2 47.863 * (real-mass ^ 1.18)
  ifelse maxhrx2 > maxhrx1 [set maxhr maxhrx2] [set maxhr maxhrx1]
  set maxhr maxhr * 10000                                          ;; in m2
  set maxhr sqrt (maxhr / pi)                                      ;; radius of maxhr in m
  set maxhr maxhr / 10                                             ;; radius in patch length (=10m)
end
;---------------------------------------------

to calc-lococost                                                   ;; calculate movement costs, specific for food type
  calc-postural
  set lococost 10.7 * (real-mass ^ 0.68) + postural                ;; costs for mammals in J/m; mass in kg after Calder 1996
  set lococost lococost / 10000                                    ;; costs in g dry biomass/m; after Nagy'99 p.263 Buchmann 2011
  set lococost lococost * 10                                       ;; lococost in patchlength (= 10m)
end
;---------------------------------------------

to calc-postural
  calc-postural-time
  set postural postural-time / (0.01398 * ((real-mass * 1000) ^ 0.217) * (e ^ (-0.002 * ((log (real-mass * 1000) e) ^ 2 )))) ;J/m
end
;----------------------------------------------

to calc-postural-time
  set postural-time (6.03 * (real-mass ^ 0.697) - 2.963 * (real-mass ^ 0.737)) ;J/s
end
;-----------------------------------------------

to calc-feed-postural
  set feedspeed postural-time * (0.71 * (real-mass ^ 0.7) / 60) ;J/g (g/s Shipley 1994)
  set feedspeed feedspeed * 100 ; g/g
end
;----------------------------------------------

to calc-feedrate [m]                                                   ; basal metabolic rate from Savage et al. 2004
 set feedrate 25.6 * (m ^ 0.737);
end
;----------------------------------------------

to life-history
  set young round (2.24 * (mass ^ (-0.13)))                       ;; allomatric number of offspring after Hamilton et al. 2011
  set mass-neonate 47.86 * (mass ^ 0.93)                          ;; allometric neonate mass after Hamilton 2011
  set mass-weaning 295.12 * (mass ^ 0.91)                         ;; allometric weaning mass after Hamilton 2011
end
;---------------------------------------------

to embryo                                                         ;; embryo growth after Rickleffs 2010
  let A exp(0.865 + 1.006 * log  mass-neonate e)
  let M0 0.0001
  let b log (A / M0) e
  let k exp(0.627 - 0.905 * log gest_period e)
  set mass-embryo A * exp(- b * exp(- k * pregcount)) / 1000
end
;-----------------------------------------------

to juvenile                                                        ;; lactation growth after Sibly 2013
  let m mass * 1000
  let jm mass-embryo * 1000
  let factor ((m / jm) ^ (1 / 3) - 1)
  let prod factor * 3 / lact_period * log ((1 - ((mass-neonate / m) ^ (1 / 3))) / (1 - ((mass-weaning / m) ^ (1 / 3)))) e
  set jm prod * jm
  set mass-embryo mass-embryo + jm / 1000
end

;-----------------------------------------------

to growth                                                          ;; ongoing growth after Sibly 2013
  let m mass * 1000
  let realm real-mass * 1000
  let factor ((m / realm) ^ (1 / 3) - 1)
  let prod factor * 3 / lact_period * log ((1 - ((mass-neonate / m) ^ (1 / 3))) / (1 - ((mass-weaning / m) ^ (1 / 3)))) e
  set realm prod * realm
  set real-mass min list mass (real-mass + realm / 1000)

  calc-lococost
  calc-feedrate real-mass
  calc-energetics
end

;-----------------------------------------------

to change_order                                                    ;; order x% of turtles according to body mass, remaining turtles are ordered randomly
  ask turtles
  [
     let zuf random-float 1.0
     let sp species
     ifelse zuf < mass_order [set order mass ] [set order random-float maxmass]
  ]
end
;-----------------------------------------------

to change_order_density_radius                                     ;; order x% of turtles according to body mass and conspecifics in the surroundings, remaining turtles are ordered randomly
  ask turtles
  [
     let zuf random-float 1.0
     let sp species
     let mhr maxhr
     ifelse zuf < mass_order [set order mass * (1 - count turtles in-radius mhr with [species = sp] / count turtles in-radius mhr)] [set order random-float maxmass]
  ]
end
;-----------------------------------------------

to one-step                                                      ;; one time step

  if ( debug = 1 ) and ( ticks > 360 ) [
    profiler:start
  ]

  set day ticks mod 360                                          ;;day of the year
  set fail 0                                                     ;;initializations
  set fail3 0
  set fail3_spec [0 0 0 0 0 0 0 0 0 0]
  set repro 0

  ask patches
  [
    set pcolor scale-color green habitat 1 0
    set spec-list [ ]
    set eaten 0
    set visited 0
  ]

  ask turtles
  [
   set age age + 1                                               ;;aging
   set fmr 0
   set preg_next false
  ]

  ifelse seasonal-resources [resource-update-seasonal-sin] [resource-update]  ;; function: update resources

  change_order_density_radius                                      ;; function: order mass_order% of turtles according to body mass
  check-hr                                                         ;; function: check if old homerange is still O.K., adapt if possible
  patch-use                                                        ;; function: save the use of different patch-types
  maintenance                                                      ;; function: energy allocation to maintenance, growth, juveniles
  foreach sort-on [ (- order)] turtles                             ;; offspring check hr dependent on mass_oder% of turtles
  [ ?1 -> ask ?1                                                   ;; ask turtles
   [
   if sex = "female"
    [
     ifelse preg = 1
      [set pregcount pregcount + 1]
        [let spec species ;
          if (real-mass >= maturity * mass) [set preg 1 set real-young young set mass-embryo 0]] ;; getting pregnant when maturity mass is reached
     if pregcount >  gest_period + lact_period                     ;; days of gestation/pregnancy PLUS days of lactation==> allometric after Hamilton et al 2011
      [ offspring set pregcount 0 set preg_next true ]
    ] ;; end if sex = female
   mort                                                            ;; mortality function
   ] ;;end ?1 -> ask ?1
   ] ;; end foreach sort

  output                                                           ;; saves output variables
  ask turtles with [preg_next = true] [set preg 0]

  if ( debug = 1 ) and ( ticks > 360 ) [
    profiler:stop
    print profiler:report
    profiler:reset
  ]
end
;-----------------------------------------------

to go                                                              ;; to go function, multiple time steps
  one-step
  tick
end
;---------------------------------------------

to lifetime_success                                               ;; calculate and save lifetime reproductive success of females (before dying)
  if sex = "female"
  [
  if species = 0 [set rep_success_0 lput (juvs) rep_success_0 set rep_success_hr_0 lput (juvshr) rep_success_hr_0]
  if species = 1 [set rep_success_1 lput (juvs) rep_success_1 set rep_success_hr_1 lput (juvshr) rep_success_hr_1]
  if species = 2 [set rep_success_2 lput (juvs) rep_success_2 set rep_success_hr_2 lput (juvshr) rep_success_hr_2]
  if species = 3 [set rep_success_3 lput (juvs) rep_success_3 set rep_success_hr_3 lput (juvshr) rep_success_hr_3]
  if species = 4 [set rep_success_4 lput (juvs) rep_success_4 set rep_success_hr_4 lput (juvshr) rep_success_hr_4]
  if species = 5 [set rep_success_5 lput (juvs) rep_success_5 set rep_success_hr_5 lput (juvshr) rep_success_hr_5]
  if species = 6 [set rep_success_6 lput (juvs) rep_success_6 set rep_success_hr_6 lput (juvshr) rep_success_hr_6]
  if species = 7 [set rep_success_7 lput (juvs) rep_success_7 set rep_success_hr_7 lput (juvshr) rep_success_hr_7]
  if species = 8 [set rep_success_8 lput (juvs) rep_success_8 set rep_success_hr_8 lput (juvshr) rep_success_hr_8]
  if species = 9 [set rep_success_9 lput (juvs) rep_success_9 set rep_success_hr_9 lput (juvshr) rep_success_hr_9]
  ]
end
;------------------------------------------------

to offspring                                                         ;; offspring, inherits everything from mother
  set repro repro + real-young                                       ;; to count overall annual reproduction
  let f focal
  let stor storage / (real-young + 1)
  set storage stor                                                   ;; mothers storage equally divided to mother and offspring
  set juvs juvs + real-young
  let juvs-found-hr 0
  hatch real-young
    [
     specify_turtle
     set age 0
     set preg 0
     set pregcount 0
     set real-mass (mass-weaning / 1000)
     calc-lococost                                                   ;; note: has to be calculated again according to real mass
     calc-feedrate real-mass
     calc-energetics
     set storage stor
     set focal f
     find-hr-offspring                                               ;; note: to offspring takes place after gestation AND lactation!
     set juvs-found-hr juvs-found-hr + 1
    ]
  set juvshr juvshr + juvs-found-hr
  calc-feedrate real-mass                                            ;; adapt mothers feedrate and max storage to without pregnancy
  set max_storage storage_addition * (294.8) * (real-mass ^ 1.19) + feedrate

end
;.............................................

to mort                                                                              ;; mortality function related to life span, normally distributed and boldness
  if shelter = 0 and random-float 1 < mortbold [ lifetime_success die ]              ;; additional mortality for bold indivdiuals (bold=> shelter= 0)
  if age > random-normal average_age ( 0.1 * average_age ) [lifetime_success die ]   ;; average life span is allometric, see above
  if density_dep [                                                                   ;; optional mortality due to density dependence
    let spec_id species
    let spec_ind count turtles with [species = spec_id]
    let capac 1000 * ( mass ^ (-0.75))
    if random-float 1.0 < (spec_ind / capac) [die]
  ]
end
;.............................................


to find-hr                                                         ;; find suitable homerange for initial distribution
  foreach sort-on [ (- order)] turtles
  [ ?1 -> ask ?1
  ;ask turtles
  [
  set maxrad maxhr
  set minfeed feedrate
  set foodshare share
  set shelter-need shelter
  set spec-color (species * 10 + 5)
  set movecost lococost
  set move-factor move-factor-allometric
  ;set feedtime feedspeed
  let spec species
  let success 0
  let moved 0
  let try 0
  let stor storage
  while [success = 0 and try < hr_try_init]                           ;; by default 100 attempts for each mammal of inital community to find suitable hr
    [
      set try try + 1
      set r-search 0
      set feed 0
      set moved 0
      ifelse small_habitat [move-to one-of patches with [habitat > 0]][move-to one-of patches with [patchquali = 1]]   ;; random search for potential hr core cell
      set core patchquali
      set feed digestion * patchfeed * foodshare - patch-move * movecost                                ;; calculate food of core cell
      set eaten 1
      if (storage + feed) >= max_storage                                           ;; if enough food in core cell - end search
       [
         set success 1
         set eaten 1                                               ;; 1 indicates patch as hr-patch forlater food reduction
         set spec-id spec                                          ;; identify patch as part of hr of species
         set spec-list lput spec-id spec-list                      ;; add species to patch-specific species list
       ]
      while [r-search < maxhr and (storage + feed) < max_storage]  ;; if not enough food in core cell - search in neighborhood
      [
       set r-search r-search + 1
       ask patches in-radius r-search with [patchfeed > 0 and eaten = 0]
        [
          if (stor + feed) >  (movecost * distance myself) and eaten = 0
          [
           set eaten 2                                             ;; 2 indicates potential use as hr-patch for later food reduction
           ifelse save                                             ;; for shy individuals there are edge effects: patches 'in the open' are less frequently visited or shorter time => reduced food intake
            [set moved moved + 2 * distance myself * move-factor + patch-move                 ;; calculate moved distance
              set feed feed + digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost]   ;;calculate food intake minus digestion and locomotion costs
            [set moved moved + (2 * distance myself * move-factor + patch-move) * (1 - shelter-need)
              set feed feed + (digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost) * (1 - shelter-need)]
        ]  ;; end ask-patches in r
       ]  ;; end while r-search
      ]

     ifelse storage + feed >= minfeed                          ;; if there was sufficient food intake, define this cell as home range core cell
       [
          set success 1
          set storage min (list max_storage (feed - minfeed))
          set patchfeed patchfeed * (1 - foodshare)
          set eaten 0
          ask patches in-radius (r-search + 1) with [eaten = 2]
           [
              ifelse save
                [set patchfeed patchfeed * (1 - foodshare)]        ;; reduce remaining food in patch
                [set patchfeed patchfeed * (1 - foodshare * (1 - shelter-need))]
              set spec-id spec                                     ;; identify patch as part of hr of species
              set spec-list lput spec-id spec-list                 ;; add species to patch-specific species list
              if spec = focal_spec [set pcolor spec-color]         ;; show only hr of focal species
              set eaten 0
           ] ;; end ask parches in-radius
          let patchafter sum [patchfeed] of patches
       ] ;; end ifelse feed >=minfeed cond1
       [
          ask patches in-radius (r-search) [set eaten 0]            ;; set back to not-eaten
       ] ;; end ifelse feed >=minfeed cond2

    ] ;; end while success = 0

    if success = 1 [                                                ;; when successfully found a home range, calculate home range size, moved distance, digestion costs, locomotion costs and storage
      set hrsize 2 * r-search
      set daymoved moved
      set fmr_digest 0.1 * (feed + moved * movecost)
      set fmr_loco moved * movecost
      set storage min (list max_storage (storage + feed))]
    if success = 0 [die]                                            ;; when unsuccessful after all triels: die

  ] ;; end ask turtles
  ]
end
;----------------------------------------------

to mother-feedrate
ask turtles
  [
    if age > 0                                                         ;; only for non-offspring
     [
     let young-mass 0
     if preg = 1                                                       ;; for gestating or lactating females calculate the summed mass of growen offspring (growth different in gestating and lactating phase)
      [
       let ym mass-embryo
       ifelse (pregcount < gest_period)
        [
         embryo
         set young-mass real-young * mass-embryo
        ]
      [
        if (pregcount < gest_period + lact_period)
         [
          juvenile
          set young-mass real-young * mass-embryo
         ]
      ] ;;end ifelse pregcount
          set younggrowth (mass-embryo - ym)                         ;; calculate the growth difference
          calc-feedrate (young-mass + real-mass)                     ;; calculate the amount of needed food for the mother
          set max_storage storage_addition * (294.8) * ((real-mass + young-mass) ^ 1.19) + feedrate    ;; calculate the storage potential of the mother
     ]
    ]
  ]
end

;----------------------------------------------
to check-hr                                                              ;; check if existing homerange is still sufficiant and reduce food
 mother-feedrate                                                         ;; function: mother individuals calculate how much food they need
 foreach sort-on [ (- order)] turtles                                    ;; ask turtles to forage, potentially ordered
  [ ?1 -> ask ?1
   [
    if age > 0                                                           ;; only for non-offspring
    [
     let m real-mass
     set real-mass real-mass + storage / 3930                            ;; adapt locomotion costs for storage
     calc-lococost
     set move-factor move-factor-allometric                              ;; initialize helper variables
     set real-mass m
     let spec species
     set maxrad maxhr
    set foodshare share
    set shelter-need shelter
    set movecost lococost
    let success 0
    set r-search 0
    set feed 0
    let stor storage
    set feed digestion * patchfeed * foodshare - patch-move * movecost   ;; calculate food of core cell
    set patchfeed patchfeed * (1 - foodshare)                            ;; reduce remaining food in core patch
    let in patchfeed * foodshare
    set eaten 1
    let moved 0
    let before visited
    if patchfeed < 0 [set patchfeed 0]
    set visited visited + 1
    if (storage + feed) >= max_storage                                                   ;; if enough food in core cell - forage and end
       [
         set success 1
         set hunger 0
         set pcolor spec-color
         set spec-id spec
         set spec-list lput spec-id spec-list
       ]
        while [r-search < maxhr and (storage + feed) < max_storage]                       ;; if not enough food in core cell - search in neighborhood
      [
       set r-search r-search + 1
         ask patches in-radius r-search with [eaten = 0 and patchfeed > 0]
        [
         if (stor + feed) >  (movecost * distance myself) and eaten = 0
         [
         if spec = focal_spec [set pcolor spec-color]                    ;; show only hr of focal species
         set eaten 2
         set spec-id spec
         set spec-list lput spec-id spec-list
         ifelse save                                                     ;; edge effect for shy individuals: patches 'in the open' are less frequently or for shorter time visited
           [
             set moved moved + 2 * distance myself * move-factor + patch-move       ;; calculate moved distance
             set feed feed + digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost   ;;calculate food intake minus digestion and locomotion costs
             set in in + patchfeed * foodshare                                      ;; calculate absolut food intake
             set patchfeed patchfeed * (1 - foodshare)                              ;; reduce resources
             set before before + visited                                            ;; increase visitation rates
             set visited visited + 1
           ]
           [
             set moved moved + (2 * distance myself * move-factor + patch-move) * (1 - shelter-need)
             set feed feed + (digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost) * (1 - shelter-need)
             set in in + (patchfeed * foodshare) * (1 - shelter-need)
             set patchfeed patchfeed * (1 - foodshare * (1 - shelter-need))  ;;reduce food in patch
             set before before + visited * (1 - shelter-need)
             set visited visited + 1 * (1 - shelter-need)
           ]
         if patchfeed < 0 [set patchfeed 0]
        ] ;; end ask patches in-radius
        ]
       ] ;; end while r-search
        set fmr_digest 0.1 * (feed / digestion + moved * movecost)           ;; calculate individual metabolism etc.
        set fmr_loco moved * movecost
        set fmr fmr + fmr_loco + fmr_digest
        set daymoved moved
        set storage min (list max_storage (storage + feed))
        set hrsize 2 * r-search
        set vis_before before
        set intake in
        ask patches in-radius r-search [set eaten 0]

      ]
    ]
  ]
end
;----------------------------------------------------------
to maintenance                                                        ;; distribute energy to maintenance, growth, and reproduction
ask turtles
  [
    set production 0
    set fmr_repro 0
    set fmr_growth 0
    if age > 0                                                         ;; only for non-offspring
     [
     let young-mass 0
     let synth 0
     if preg = 1                                                       ;; gestating or lactating females calculate costs for growing offspring
      [
       set young-mass real-young * mass-embryo
       set synth (real-young * younggrowth * growthcost)
      ]
      let old feedrate
      calc-feedrate (young-mass + real-mass)
      let act_feedrate feedrate
      set feedrate old
      ifelse storage > (act_feedrate + synth)                          ;; if there is enough energy for maintenenace and possibly reproduction
        [
          let m real-mass
          growth                                                       ;; calculate possible growth of adult
          let diff real-mass - m
          let diffcost diff * growthcost
          set storage storage - act_feedrate - synth                   ;; reduce energy by maintenance and possibly reproduction
          set production production + (synth / growthcost)
          ifelse storage > diffcost                                    ;; if there is enough energy for growth: grow, reduce energy, and recalculate allometries
          [set storage storage - diffcost
           set production production + diff
           calc-feedrate (real-mass + young-mass)
           set act_feedrate feedrate
           calc-feedrate (real-mass)
           set fmr fmr + act_feedrate + diffcost + synth
           set fmr_basal feedrate
           set fmr_repro synth + (act_feedrate - feedrate)
           set fmr_growth diffcost]
          [set real-mass m                                              ;; if there is not wnough energy for growth, there is no growth, update metabolism variables
              calc-lococost
              calc-feedrate (young-mass + real-mass)
              calc-energetics
              set act_feedrate feedrate
              calc-feedrate (real-mass)
              set fmr_basal feedrate
              set fmr_repro synth + (act_feedrate - feedrate)
              set fmr fmr + act_feedrate + synth
            ]
        ]
      [                                                                 ;; if there is not enough energy for maintenance and possibly reproduction
        ifelse preg = 1                                                 ;; for reproducing individuals first reduce reproduction investment if there is not enough energy for maintenance and reproduction
        [
          ifelse pregcount > gest_period                                ;; during lactation period loose one juvenile after the other until there is enough energy for the remaining juveniles
            [
                while [storage < act_feedrate + (real-young * younggrowth * growthcost) and real-young > 0]
              [
              set real-young real-young - 1
              set failed-lact failed-lact + 1
              set young-mass real-young * mass-embryo
                set old feedrate
                calc-feedrate (young-mass + real-mass)
                set act_feedrate feedrate
                set feedrate old
              ]
              ifelse storage > act_feedrate + (real-young * younggrowth * growthcost)
              [set storage storage - act_feedrate - (real-young * younggrowth * growthcost)
               set production production + real-young * younggrowth
               set fmr fmr + act_feedrate + (real-young * younggrowth * growthcost)
               set act_feedrate feedrate
               calc-feedrate (real-mass)
               set fmr_basal feedrate
               set fmr_repro (real-young * younggrowth * growthcost) + (act_feedrate - feedrate)
               if real-young = 0 [set preg 0 set pregcount 0 ]
              ]
              [                                                            ;; if there is not enough energy for maintenance even after loosing all juveniles: die
                set mort_order replace-item species mort_order (sentence order (item species mort_order))
                set mort_stor replace-item species mort_stor (sentence (storage / max_storage) (item species mort_stor))
                set mort_food replace-item species mort_food (sentence (competition / forage-patches) (item species mort_food))
                set mort_hr replace-item species mort_hr (sentence hrsize (item species mort_hr))
                set fail3_spec replace-item species fail3_spec (item species fail3_spec + 1)
                set preg 0
                set pregcount 0
                set storage 0
                set fail3 fail3 + 1
                lifetime_success
                die
              ]
            ]
        [                                                                   ;; during gestation period loose the entire pregnancy if there is not enough energy for all
          set preg 0
          set pregcount 0
          calc-feedrate real-mass
          ;set feedrate 55.1 * (real-mass ^ 0.74)
          ifelse storage > feedrate
          [set storage storage - feedrate
           set fmr fmr + feedrate
           set fmr_basal feedrate]
          [                                                                 ;; if there is not enough energy for maintenance even after loosing the pregnancy: die
           set mort_order replace-item species mort_order (sentence order (item species mort_order))
           set mort_stor replace-item species mort_stor (sentence (storage / max_storage) (item species mort_stor))
           set mort_food replace-item species mort_food (sentence (competition / forage-patches) (item species mort_food))
           set mort_hr replace-item species mort_hr (sentence hrsize (item species mort_hr))
           set fail3_spec replace-item species fail3_spec (item species fail3_spec + 1)
           set storage 0
           set fail3 fail3 + 1
           lifetime_success
           die
          ]
        ]
        ]
        [                                                                     ;; non-reproducing individuals die if they have not enough energy for maintenance
           set mort_order replace-item species mort_order (sentence order (item species mort_order))
           set mort_stor replace-item species mort_stor (sentence (storage / max_storage) (item species mort_stor))
           set mort_food replace-item species mort_food (sentence (competition / forage-patches) (item species mort_food))
           set mort_hr replace-item species mort_hr (sentence hrsize (item species mort_hr))
           set fail3_spec replace-item species fail3_spec (item species fail3_spec + 1)
           set storage 0
           set fail3 fail3 + 1
           lifetime_success
           die
        ]

      ]
    if intake > 0
        [
        set balance (intake * digestion) / fmr
        set balanceprod production / intake
        ]
    ]
  ]
  ;]
end
;----------------------------------------------

to find-hr-offspring                                                         ;; offspring searches for new homerange similar to home range search in initialization
  set maxrad maxhr
  set minfeed feedrate
  set foodshare share
  set shelter-need shelter
  set spec-color (species * 10 + 5)
  set movecost lococost
  set move-factor move-factor-allometric
  set max_nat_disp 3.31 * (mass ^ 0.65)                                      ;; allometric maximum natal dispersal after Sutherland et al (2000) in km
  let spec species
  let success 0
  let try 0
  let xx xcor
  let yy ycor
  let moved 0
  let stor storage
  while [success = 0 and try < hr_try_juv]                                   ;; hr_try_juv attempts for each mammal to find suitable hr
    [
      setxy xx yy
      set try try + 1
      let patch-before patch-here
      ifelse small_habitat
          [move-to one-of patches in-radius (100 * max_nat_disp) with [habitat > 0]]     ;; random search for hr core cell in maximum natal dispersal distance in grid cell length = 10m
          [move-to one-of patches in-radius (100 * max_nat_disp) with [patchquali = 1]]
      set core patchquali
      set r-search 0
      set moved 0
      set feed digestion * patchfeed * foodshare - lococost * distance patch-before - patch-move * movecost   ;; food intake in core cell minus digestion and locomotion costs
      set eaten 1
      if (storage + feed) >= max_storage                                     ;; if enough food in core cell - end
       [
         set success 1
         set pcolor spec-color
         set eaten 1                                                         ;; 1 indicates patch as hr-patch for later food reduction
         set spec-id spec                                                    ;; identify patch as part of hr of species
         set spec-list lput spec-id spec-list                                ;; add species to patch-specific species list
       ]
      while [r-search < maxhr and (storage + feed) < max_storage]
     [
       set r-search r-search + 1
       ask patches in-radius (r-search) with [patchfeed > 0 and eaten = 0]
        [
          if (stor + feed) >  (movecost * distance myself) and eaten = 0
          [
          set eaten 2                                                        ;;  2 indicates potential use as hr-patch for later food reduction
          ifelse save                                                        ;; edge effect for shy individuals: patches 'in the open' are less frequently visited or fot shorter time
            [set moved moved + 2 * distance myself * move-factor + patch-move                                                              ;;calculate moved distance
              set feed feed + digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost]    ;;calculate food intake minus digestion and locomotion costs
            [set moved moved + (2 * distance myself * move-factor + patch-move) * (1 - shelter-need)
              set feed feed + (digestion * patchfeed * foodshare - 2 * movecost * distance myself * move-factor - patch-move * movecost) * (1 - shelter-need)]
        ] ;; end ask patches in radius
       ] ;; end while r-search
      ]

      ifelse (storage + feed) >= minfeed                                   ;; if enough food found
       [
         set success 1
         set storage  min ( list max_storage (storage + feed - minfeed))
         set patchfeed patchfeed * (1 - foodshare)
         set eaten 0
         set daymoved moved
         ask patches in-radius (r-search + 1) with [eaten = 2]
           [
              ifelse save
                [set patchfeed patchfeed * (1 - foodshare)]                  ;; reduce remaining food in patch
                [set patchfeed patchfeed * (1 - foodshare * (1 - shelter-need))]
              set spec-id spec                                               ;; identify patch as part of hr of species
              set spec-list lput spec-id spec-list                           ;; add species to patch-specific species list
              if spec = focal_spec [set pcolor spec-color]                   ;; show only hr of focal species
              set eaten 0
             ;; end if eaten = 1
           ] ;; end ask patches in-radius
         ] ;; end ifelse feed >=minfeed cond 1
         [
          ask patches in-radius (r-search) [set eaten 0]                     ;; set back to not-eaten
         ] ;;end ifelse feed >=minfeed cond 2
    ] ;; end while success = 0
   if success = 1 [set hrsize 2 * r-search set suc_juv replace-item species suc_juv (item species suc_juv + 1)]     ;;successful individuals are counted
   if success = 0 [set fail (fail + 1) die]                                                    ;; unsuccessful individuals die or emigrate from area
end
;----------------------------------------------

to patch-use                                                                        ;; calculate the use of small and large patches, currently not in use
  ask turtles
  [
    set forage_large count (patches in-radius (hrsize / 2) with [patchquali = 1])
    set forage_small count (patches in-radius (hrsize / 2) with [patchquali = 2])

    set competition sum [length spec-list] of patches in-radius (hrsize / 2)
    set forage-patches count (patches in-radius (hrsize / 2) with [habitat = 1])
  ]
end
;----------------------------------------------


to init-out        ;;initialize output variables

  set mort_order [0 0 0 0 0 0 0 0 0 0]
  set mort_stor [0 0 0 0 0 0 0 0 0 0]
  set mort_food [0 0 0 0 0 0 0 0 0 0]
  set mort_hr [0 0 0 0 0 0 0 0 0 0]
  set suc_juv [0 0 0 0 0 0 0 0 0 0]

  set rep_success_0 []
  set rep_success_1 []
  set rep_success_2 []
  set rep_success_3 []
  set rep_success_4 []
  set rep_success_5 []
  set rep_success_6 []
  set rep_success_7 []
  set rep_success_8 []
  set rep_success_9 []

  set rep_success_hr_0 []
  set rep_success_hr_1 []
  set rep_success_hr_2 []
  set rep_success_hr_3 []
  set rep_success_hr_4 []
  set rep_success_hr_5 []
  set rep_success_hr_6 []
  set rep_success_hr_7 []
  set rep_success_hr_8 []
  set rep_success_hr_9 []

  set number []
  set number0 []
  set number1 []
  set number2 []
  set number3 []
  set number4 []
  set number5 []
  set number6 []
  set number7 []
  set number8 []
  set number9 []

  set patches0 []
  set patches1 []
  set patches2 []
  set patches3 []
  set patches4 []
  set patches5 []
  set patches6 []
  set patches7 []
  set patches8 []
  set patches9 []

  set compet0 []
  set compet1 []
  set compet2 []
  set compet3 []
  set compet4 []
  set compet5 []
  set compet6 []
  set compet7 []
  set compet8 []
  set compet9 []

  set order0 []
  set order1 []
  set order2 []
  set order3 []
  set order4 []
  set order5 []
  set order6 []
  set order7 []
  set order8 []
  set order9 []

  set vis_before0 []
  set vis_before1 []
  set vis_before2 []
  set vis_before3 []
  set vis_before4 []
  set vis_before5 []
  set vis_before6 []
  set vis_before7 []
  set vis_before8 []
  set vis_before9 []

  set hr []
  set hr0 []
  set hr1 []
  set hr2 []
  set hr3 []
  set hr4 []
  set hr5 []
  set hr6 []
  set hr7 []
  set hr8 []
  set hr9 []

  set maxhrs []
  set maxhr0 []
  set maxhr1 []
  set maxhr2 []
  set maxhr3 []
  set maxhr4 []
  set maxhr5 []
  set maxhr6 []
  set maxhr7 []
  set maxhr8 []
  set maxhr9 []

  set countmaxhr []
  set countmaxhr0 []
  set countmaxhr1 []
  set countmaxhr2 []
  set countmaxhr3 []
  set countmaxhr4 []
  set countmaxhr5 []
  set countmaxhr6 []
  set countmaxhr7 []
  set countmaxhr8 []
  set countmaxhr9 []

  set ages []
  set age0 []
  set age1 []
  set age2 []
  set age3 []
  set age4 []
  set age5 []
  set age6 []
  set age7 []
  set age8 []
  set age9 []

  set stors []
  set stor0 []
  set stor1 []
  set stor2 []
  set stor3 []
  set stor4 []
  set stor5 []
  set stor6 []
  set stor7 []
  set stor8 []
  set stor9 []

  set rm []
  set rm0 []
  set rm1 []
  set rm2 []
  set rm3 []
  set rm4 []
  set rm5 []
  set rm6 []
  set rm7 []
  set rm8 []
  set rm9 []

  set fmr0 []
  set fmr1 []
  set fmr2 []
  set fmr3 []
  set fmr4 []
  set fmr5 []
  set fmr6 []
  set fmr7 []
  set fmr8 []
  set fmr9 []

  set loco0 []
  set loco1 []
  set loco2 []
  set loco3 []
  set loco4 []
  set loco5 []
  set loco6 []
  set loco7 []
  set loco8 []
  set loco9 []

  set repro0 []
  set repro1 []
  set repro2 []
  set repro3 []
  set repro4 []
  set repro5 []
  set repro6 []
  set repro7 []
  set repro8 []
  set repro9 []

  set grow0 []
  set grow1 []
  set grow2 []
  set grow3 []
  set grow4 []
  set grow5 []
  set grow6 []
  set grow7 []
  set grow8 []
  set grow9 []

  set basal0 []
  set basal1 []
  set basal2 []
  set basal3 []
  set basal4 []
  set basal5 []
  set basal6 []
  set basal7 []
  set basal8 []
  set basal9 []

  set digest0 []
  set digest1 []
  set digest2 []
  set digest3 []
  set digest4 []
  set digest5 []
  set digest6 []
  set digest7 []
  set digest8 []
  set digest9 []

  set prod0 []
  set prod1 []
  set prod2 []
  set prod3 []
  set prod4 []
  set prod5 []
  set prod6 []
  set prod7 []
  set prod8 []
  set prod9 []

  set in0 []
  set in1 []
  set in2 []
  set in3 []
  set in4 []
  set in5 []
  set in6 []
  set in7 []
  set in8 []
  set in9 []

  set inpatch0 []
  set inpatch1 []
  set inpatch2 []
  set inpatch3 []
  set inpatch4 []
  set inpatch5 []
  set inpatch6 []
  set inpatch7 []
  set inpatch8 []
  set inpatch9 []

  set inhr0 []
  set inhr1 []
  set inhr2 []
  set inhr3 []
  set inhr4 []
  set inhr5 []
  set inhr6 []
  set inhr7 []
  set inhr8 []
  set inhr9 []

  set balance0 []
  set balance1 []
  set balance2 []
  set balance3 []
  set balance4 []
  set balance5 []
  set balance6 []
  set balance7 []
  set balance8 []
  set balance9 []

  set prodbalance0 []
  set prodbalance1 []
  set prodbalance2 []
  set prodbalance3 []
  set prodbalance4 []
  set prodbalance5 []
  set prodbalance6 []
  set prodbalance7 []
  set prodbalance8 []
  set prodbalance9 []

set move0 []
set move1 []
set move2 []
set move3 []
set move4 []
set move5 []
set move6 []
set move7 []
set move8 []
set move9 []

  if drought_scenario
    [
      set drought_hr0 []
      set drought_hr1 []
      set drought_hr2 []
      set drought_hr3 []
      set drought_hr4 []
      set drought_hr5 []
      set drought_hr6 []
      set drought_hr7 []
      set drought_hr8 []
      set drought_hr9 []

      set drought_fmr0 []
      set drought_fmr1 []
      set drought_fmr2 []
      set drought_fmr3 []
      set drought_fmr4 []
      set drought_fmr5 []
      set drought_fmr6 []
      set drought_fmr7 []
      set drought_fmr8 []
      set drought_fmr9 []

      set drought_in0 []
      set drought_in1 []
      set drought_in2 []
      set drought_in3 []
      set drought_in4 []
      set drought_in5 []
      set drought_in6 []
      set drought_in7 []
      set drought_in8 []
      set drought_in9 []

      set drought_stor0 []
      set drought_stor1 []
      set drought_stor2 []
      set drought_stor3 []
      set drought_stor4 []
      set drought_stor5 []
      set drought_stor6 []
      set drought_stor7 []
      set drought_stor8 []
      set drought_stor9 []
    ]

end

to output                        ;;calculate output variables for all species

  let specs [species] of turtles
  let specs_unique remove-duplicates specs
  set spec_num length specs_unique

set number lput count turtles number
set number0 lput count turtles with [species = 0] number0
set number1 lput count turtles with [species = 1] number1
set number2 lput count turtles with [species = 2] number2
set number3 lput count turtles with [species = 3] number3
set number4 lput count turtles with [species = 4] number4
set number5 lput count turtles with [species = 5] number5
set number6 lput count turtles with [species = 6] number6
set number7 lput count turtles with [species = 7] number7
set number8 lput count turtles with [species = 8] number8
set number9 lput count turtles with [species = 9] number9


  if any? turtles with [age > 1]
  [
set hr lput mean [hrsize] of turtles with [age > 1] hr
set maxhrs lput mean [hrsize / maxhr] of turtles with [age > 1] maxhrs
set countmaxhr lput count turtles with [age > 1 and hrsize / maxhr > 0.9] countmaxhr
set stors lput mean [storage / max_storage] of turtles with [age > 1] stors
set ages lput mean [age] of turtles with [age > 1] ages
  ]

  if any? turtles with [species = 0 and age > 1]
  [
set hr0 lput mean [hrsize] of turtles with [species = 0 and age > 1] hr0
set maxhr0 lput mean [hrsize / maxhr] of turtles with [species = 0 and age > 1] maxhr0
set countmaxhr0 lput count turtles with [species = 0 and age > 1 and hrsize / maxhr > 0.9] countmaxhr0
set age0 lput mean [age] of turtles with [species = 0 and age > 1] age0
set stor0 lput mean [storage / max_storage] of turtles with [species = 0 and age > 1] stor0
set move0 lput mean [daymoved] of turtles with [species = 0 and age > 1] move0
set patches0 lput mean [forage-patches] of turtles with [species = 0 and age > 1] patches0
set compet0 lput mean [competition / forage-patches] of turtles with [species = 0 and age > 1] compet0
set order0 lput mean [order] of turtles with [species = 0 and age > 1] order0
set vis_before0 lput mean [vis_before / forage-patches] of turtles with [species = 0 and age > 1] vis_before0
set fmr0 lput mean [fmr] of turtles with [species = 0 and age > 1] fmr0
set loco0 lput mean [fmr_loco / fmr] of turtles with [species = 0 and age > 1] loco0
set repro0 lput mean [fmr_repro / fmr] of turtles with [species = 0 and age > 1] repro0
set grow0 lput mean [fmr_growth / fmr] of turtles with [species = 0 and age > 1] grow0
set basal0 lput mean [fmr_basal / fmr] of turtles with [species = 0 and age > 1] basal0
set digest0 lput mean [fmr_digest / fmr] of turtles with [species = 0 and age > 1] digest0
set prod0 lput mean [production] of turtles with [species = 0 and age > 1] prod0
set in0 lput mean [intake] of turtles with [species = 0 and age > 1] in0
set inpatch0 lput mean [intake / forage-patches] of turtles with [species = 0 and age > 1] inpatch0
set inhr0 lput mean [intake / (hrsize + 1)] of turtles with [species = 0 and age > 1] inhr0
set balance0 lput median [balance] of turtles with [species = 0 and age > 1] balance0
set prodbalance0 lput median [balanceprod] of turtles with [species = 0 and age > 1] prodbalance0
  if drought_scenario
    [
     if drought
     [
        set drought_hr0 lput mean [hrsize] of turtles with [species = 0 and age > 1] drought_hr0
        set drought_fmr0 lput mean [fmr] of turtles with [species = 0 and age > 1] drought_fmr0
        set drought_in0 lput mean [intake] of turtles with [species = 0 and age > 1] drought_in0
        set drought_stor0 lput mean [storage / max_storage] of turtles with [species = 0 and age > 1] drought_stor0
     ]
    ]
  ]


  if any? turtles with [species = 1 and age > 1]
  [
set hr1 lput mean [hrsize] of turtles with [species = 1 and age > 1] hr1
set maxhr1 lput mean [hrsize / maxhr] of turtles with [species = 1 and age > 1] maxhr1
set countmaxhr1 lput count turtles with [species = 1 and age > 1 and hrsize / maxhr > 0.9] countmaxhr1
set age1 lput mean [age] of turtles with [species = 1 and age > 1] age1
set stor1 lput mean [storage / max_storage] of turtles with [species = 1 and age > 1] stor1
set move1 lput mean [daymoved] of turtles with [species = 1 and age > 1] move1
set patches1 lput mean [forage-patches] of turtles with [species = 1 and age > 1] patches1
set compet1 lput mean [competition / forage-patches] of turtles with [species = 1 and age > 1] compet1
set order1 lput mean [order] of turtles with [species = 1 and age > 1] order1
set vis_before1 lput mean [vis_before / forage-patches] of turtles with [species = 1 and age > 1] vis_before1
set fmr1 lput mean [fmr] of turtles with [species = 1 and age > 1] fmr1
set loco1 lput mean [fmr_loco / fmr] of turtles with [species = 1 and age > 1] loco1
set repro1 lput mean [fmr_repro / fmr] of turtles with [species = 1 and age > 1] repro1
set grow1 lput mean [fmr_growth / fmr] of turtles with [species = 1 and age > 1] grow1
set basal1 lput mean [fmr_basal / fmr] of turtles with [species = 1 and age > 1] basal1
set digest1 lput mean [fmr_digest / fmr] of turtles with [species = 1 and age > 1] digest1
set prod1 lput mean [production] of turtles with [species = 1 and age > 1] prod1
set in1 lput mean [intake] of turtles with [species = 1 and age > 1] in1
set inpatch1 lput mean [intake / forage-patches] of turtles with [species = 1 and age > 1] inpatch1
set inhr1 lput mean [intake / (hrsize + 1)] of turtles with [species = 1 and age > 1] inhr1
set balance1 lput median [balance] of turtles with [species = 1 and age > 1] balance1
set prodbalance1 lput median [balanceprod] of turtles with [species = 1 and age > 1] prodbalance1
    if drought_scenario
    [
     if drought
     [
        set drought_hr1 lput mean [hrsize] of turtles with [species = 1 and age > 1] drought_hr1
        set drought_fmr1 lput mean [fmr] of turtles with [species = 1 and age > 1] drought_fmr1
        set drought_in1 lput mean [intake] of turtles with [species = 1 and age > 1] drought_in1
        set drought_stor1 lput mean [storage / max_storage] of turtles with [species = 1 and age > 1] drought_stor1
     ]
    ]
  ]


  if any? turtles with [species = 2 and age > 1]
  [
set hr2 lput mean [hrsize] of turtles with [species = 2 and age > 1] hr2
set maxhr2 lput mean [hrsize / maxhr] of turtles with [species = 2 and age > 1] maxhr2
set countmaxhr2 lput count turtles with [species = 2 and age > 1 and hrsize / maxhr > 0.9] countmaxhr2
set age2 lput mean [age] of turtles with [species = 2 and age > 1] age2
set stor2 lput mean [storage / max_storage] of turtles with [species = 2 and age > 1] stor2
set move2 lput mean [daymoved] of turtles with [species = 2 and age > 1] move2
set patches2 lput mean [forage-patches] of turtles with [species = 2 and age > 1] patches2
set compet2 lput mean [competition / forage-patches] of turtles with [species = 2 and age > 1] compet2
set order2 lput mean [order] of turtles with [species = 2 and age > 1] order2
set vis_before2 lput mean [vis_before / forage-patches] of turtles with [species = 2 and age > 1] vis_before2
set fmr2 lput mean [fmr] of turtles with [species = 2 and age > 1] fmr2
set loco2 lput mean [fmr_loco / fmr] of turtles with [species = 2 and age > 1] loco2
set repro2 lput mean [fmr_repro / fmr] of turtles with [species = 2 and age > 1] repro2
set grow2 lput mean [fmr_growth / fmr] of turtles with [species = 2 and age > 1] grow2
set basal2 lput mean [fmr_basal / fmr] of turtles with [species = 2 and age > 1] basal2
set digest2 lput mean [fmr_digest / fmr] of turtles with [species = 2 and age > 1] digest2
set prod2 lput mean [production] of turtles with [species = 2 and age > 1] prod2
set in2 lput mean [intake] of turtles with [species = 2 and age > 1] in2
set inpatch2 lput mean [intake / forage-patches] of turtles with [species = 2 and age > 1] inpatch2
set inhr2 lput mean [intake / (hrsize + 1)] of turtles with [species = 2 and age > 1] inhr2
set balance2 lput median [balance] of turtles with [species = 2 and age > 1] balance2
set prodbalance2 lput median [balanceprod] of turtles with [species = 2 and age > 1] prodbalance2
    if drought_scenario
    [
     if drought
     [
        set drought_hr2 lput mean [hrsize] of turtles with [species = 2 and age > 1] drought_hr2
        set drought_fmr2 lput mean [fmr] of turtles with [species = 2 and age > 1] drought_fmr2
        set drought_in2 lput mean [intake] of turtles with [species = 2 and age > 1] drought_in2
        set drought_stor2 lput mean [storage / max_storage] of turtles with [species = 2 and age > 1] drought_stor2
     ]
    ]
  ]


  if any? turtles with [species = 3 and age > 1]
  [
set hr3 lput mean [hrsize] of turtles with [species = 3 and age > 1] hr3
set maxhr3 lput mean [hrsize / maxhr] of turtles with [species = 3 and age > 1] maxhr3
set countmaxhr3 lput count turtles with [species = 3 and age > 1 and hrsize / maxhr > 0.9] countmaxhr3
set age3 lput mean [age] of turtles with [species = 3 and age > 1] age3
set stor3 lput mean [storage / max_storage] of turtles with [species = 3 and age > 1] stor3
set move3 lput mean [daymoved] of turtles with [species = 3 and age > 1] move3
set patches3 lput mean [forage-patches] of turtles with [species = 3 and age > 1] patches3
set compet3 lput mean [competition / forage-patches] of turtles with [species = 3 and age > 1] compet3
set order3 lput mean [order] of turtles with [species = 3 and age > 1] order3
set vis_before3 lput mean [vis_before / forage-patches] of turtles with [species = 3 and age > 1] vis_before3
set fmr3 lput mean [fmr] of turtles with [species = 3 and age > 1] fmr3
set loco3 lput mean [fmr_loco / fmr] of turtles with [species = 3 and age > 1] loco3
set repro3 lput mean [fmr_repro / fmr] of turtles with [species = 3 and age > 1] repro3
set grow3 lput mean [fmr_growth / fmr] of turtles with [species = 3 and age > 1] grow3
set basal3 lput mean [fmr_basal / fmr] of turtles with [species = 3 and age > 1] basal3
set digest3 lput mean [fmr_digest / fmr] of turtles with [species = 3 and age > 1] digest3
set prod3 lput mean [production] of turtles with [species = 3 and age > 1] prod3
set in3 lput mean [intake] of turtles with [species = 3 and age > 1] in3
set inpatch3 lput mean [intake / forage-patches] of turtles with [species = 3 and age > 1] inpatch3
set inhr3 lput mean [intake / (hrsize + 1)] of turtles with [species = 3 and age > 1] inhr3
set balance3 lput median [balance] of turtles with [species = 3 and age > 1] balance3
set prodbalance3 lput median [balanceprod] of turtles with [species = 3 and age > 1] prodbalance3
    if drought_scenario
    [
     if drought
     [
        set drought_hr3 lput mean [hrsize] of turtles with [species = 3 and age > 1] drought_hr3
        set drought_fmr3 lput mean [fmr] of turtles with [species = 3 and age > 1] drought_fmr3
        set drought_in3 lput mean [intake] of turtles with [species = 3 and age > 1] drought_in3
        set drought_stor3 lput mean [storage / max_storage] of turtles with [species = 3 and age > 1] drought_stor3
     ]
    ]
  ]


  if any? turtles with [species = 4 and age > 1]
  [
set hr4 lput mean [hrsize] of turtles with [species = 4 and age > 1] hr4
set maxhr4 lput mean [hrsize / maxhr] of turtles with [species = 4 and age > 1] maxhr4
set countmaxhr4 lput count turtles with [species = 4 and age > 1 and hrsize / maxhr > 0.9] countmaxhr4
set age4 lput mean [age] of turtles with [species = 4 and age > 1] age4
set stor4 lput mean [storage / max_storage] of turtles with [species = 4 and age > 1] stor4
set move4 lput mean [daymoved] of turtles with [species = 4 and age > 1] move4
set patches4 lput mean [forage-patches] of turtles with [species = 4 and age > 1] patches4
set compet4 lput mean [competition / forage-patches] of turtles with [species = 4 and age > 1] compet4
set order4 lput mean [order] of turtles with [species = 4 and age > 1] order4
set vis_before4 lput mean [vis_before / forage-patches] of turtles with [species = 4 and age > 1] vis_before4
set fmr4 lput mean [fmr] of turtles with [species = 4 and age > 1] fmr4
set loco4 lput mean [fmr_loco / fmr] of turtles with [species = 4 and age > 1] loco4
set repro4 lput mean [fmr_repro / fmr] of turtles with [species = 4 and age > 1] repro4
set grow4 lput mean [fmr_growth / fmr] of turtles with [species = 4 and age > 1] grow4
set basal4 lput mean [fmr_basal / fmr] of turtles with [species = 4 and age > 1] basal4
set digest4 lput mean [fmr_digest / fmr] of turtles with [species = 4 and age > 1] digest4
set prod4 lput mean [production] of turtles with [species = 4 and age > 1] prod4
set in4 lput mean [intake] of turtles with [species = 4 and age > 1] in4
set inpatch4 lput mean [intake / forage-patches] of turtles with [species = 4 and age > 1] inpatch4
set inhr4 lput mean [intake / (hrsize + 1)] of turtles with [species = 4 and age > 1] inhr4
set balance4 lput median [balance] of turtles with [species = 4 and age > 1] balance4
set prodbalance4 lput median [balanceprod] of turtles with [species = 4 and age > 1] prodbalance4
    if drought_scenario
    [
     if drought
     [
        set drought_hr4 lput mean [hrsize] of turtles with [species = 4 and age > 1] drought_hr4
        set drought_fmr4 lput mean [fmr] of turtles with [species = 4 and age > 1] drought_fmr4
        set drought_in4 lput mean [intake] of turtles with [species = 4 and age > 1] drought_in4
        set drought_stor4 lput mean [storage / max_storage] of turtles with [species = 4 and age > 1] drought_stor4
     ]
    ]
  ]


  if any? turtles with [species = 5 and age > 1]
  [
set hr5 lput mean [hrsize] of turtles with [species = 5 and age > 1] hr5
set maxhr5 lput mean [hrsize / maxhr] of turtles with [species = 5 and age > 1] maxhr5
set countmaxhr5 lput count turtles with [species = 5 and age > 1 and hrsize / maxhr > 0.9] countmaxhr5
set age5 lput mean [age] of turtles with [species = 5 and age > 1] age5
set stor5 lput mean [storage / max_storage] of turtles with [species = 5 and age > 1] stor5
set move5 lput mean [daymoved] of turtles with [species = 5 and age > 1] move5
set patches5 lput mean [forage-patches] of turtles with [species = 5 and age > 1] patches5
set compet5 lput mean [competition / forage-patches] of turtles with [species = 5 and age > 1] compet5
set order5 lput mean [order] of turtles with [species = 5 and age > 1] order5
set vis_before5 lput mean [vis_before / forage-patches] of turtles with [species = 5 and age > 1] vis_before5
set fmr5 lput mean [fmr] of turtles with [species = 5 and age > 1] fmr5
set loco5 lput mean [fmr_loco / fmr] of turtles with [species = 5 and age > 1] loco5
set repro5 lput mean [fmr_repro / fmr] of turtles with [species = 5 and age > 1] repro5
set grow5 lput mean [fmr_growth / fmr] of turtles with [species = 5 and age > 1] grow5
set basal5 lput mean [fmr_basal / fmr] of turtles with [species = 5 and age > 1] basal5
set digest5 lput mean [fmr_digest / fmr] of turtles with [species = 5 and age > 1] digest5
set prod5 lput mean [production] of turtles with [species = 5 and age > 1] prod5
set in5 lput mean [intake] of turtles with [species = 5 and age > 1] in5
set inpatch5 lput mean [intake / forage-patches] of turtles with [species = 5 and age > 1] inpatch5
set inhr5 lput mean [intake / (hrsize + 1)] of turtles with [species = 5 and age > 1] inhr5
set balance5 lput median [balance] of turtles with [species = 5 and age > 1] balance5
set prodbalance5 lput median [balanceprod] of turtles with [species = 5 and age > 1] prodbalance5
    if drought_scenario
    [
     if drought
     [
        set drought_hr5 lput mean [hrsize] of turtles with [species = 5 and age > 1] drought_hr5
        set drought_fmr5 lput mean [fmr] of turtles with [species = 5 and age > 1] drought_fmr5
        set drought_in5 lput mean [intake] of turtles with [species = 5 and age > 1] drought_in5
        set drought_stor5 lput mean [storage / max_storage] of turtles with [species = 5 and age > 1] drought_stor5
     ]
    ]
  ]

  if any? turtles with [species = 6 and age > 1]
  [
set hr6 lput mean [hrsize] of turtles with [species = 6 and age > 1] hr6
set maxhr6 lput mean [hrsize / maxhr] of turtles with [species = 6 and age > 1] maxhr6
set countmaxhr6 lput count turtles with [species = 6 and age > 1 and hrsize / maxhr > 0.9] countmaxhr6
set age6 lput mean [age] of turtles with [species = 6 and age > 1] age6
set stor6 lput mean [storage / max_storage] of turtles with [species = 6 and age > 1] stor6
set move6 lput mean [daymoved] of turtles with [species = 6 and age > 1] move6
set patches6 lput mean [forage-patches] of turtles with [species = 6 and age > 1] patches6
set compet6 lput mean [competition / forage-patches] of turtles with [species = 6 and age > 1] compet6
set order6 lput mean [order] of turtles with [species = 6 and age > 1] order6
set vis_before6 lput mean [vis_before / forage-patches] of turtles with [species = 6 and age > 1] vis_before6
set fmr6 lput mean [fmr] of turtles with [species = 6 and age > 1] fmr6
set loco6 lput mean [fmr_loco / fmr] of turtles with [species = 6 and age > 1] loco6
set repro6 lput mean [fmr_repro / fmr] of turtles with [species = 6 and age > 1] repro6
set grow6 lput mean [fmr_growth / fmr] of turtles with [species = 6 and age > 1] grow6
set basal6 lput mean [fmr_basal / fmr] of turtles with [species = 6 and age > 1] basal6
set digest6 lput mean [fmr_digest / fmr] of turtles with [species = 6 and age > 1] digest6
set prod6 lput mean [production] of turtles with [species = 6 and age > 1] prod6
set in6 lput mean [intake] of turtles with [species = 6 and age > 1] in6
set inpatch6 lput mean [intake / forage-patches] of turtles with [species = 6 and age > 1] inpatch6
set inhr6 lput mean [intake / (hrsize + 1)] of turtles with [species = 6 and age > 1] inhr6
set balance6 lput median [balance] of turtles with [species = 6 and age > 1] balance6
set prodbalance6 lput median [balanceprod] of turtles with [species = 6 and age > 1] prodbalance6
    if drought_scenario
    [
     if drought
     [
        set drought_hr6 lput mean [hrsize] of turtles with [species = 6 and age > 1] drought_hr6
        set drought_fmr6 lput mean [fmr] of turtles with [species = 6 and age > 1] drought_fmr6
        set drought_in6 lput mean [intake] of turtles with [species = 6 and age > 1] drought_in6
        set drought_stor6 lput mean [storage / max_storage] of turtles with [species = 6 and age > 1] drought_stor6
     ]
    ]
  ]


  if any? turtles with [species = 7 and age > 1]
  [
set hr7 lput mean [hrsize] of turtles with [species = 7 and age > 1] hr7
set maxhr7 lput mean [hrsize / maxhr] of turtles with [species = 7 and age > 1] maxhr7
set countmaxhr7 lput count turtles with [species = 7 and age > 1 and hrsize / maxhr > 0.9] countmaxhr7
set age7 lput mean [age] of turtles with [species = 7 and age > 1] age7
set move7 lput mean [daymoved] of turtles with [species = 7 and age > 1] move7
set stor7 lput mean [storage / max_storage] of turtles with [species = 7 and age > 1] stor7
set patches7 lput mean [forage-patches] of turtles with [species = 7 and age > 1] patches7
set compet7 lput mean [competition / forage-patches] of turtles with [species = 7 and age > 1] compet7
set order7 lput mean [order] of turtles with [species = 7 and age > 1] order7
set vis_before7 lput mean [vis_before / forage-patches] of turtles with [species = 7 and age > 1] vis_before7
set fmr7 lput mean [fmr] of turtles with [species = 7 and age > 1] fmr7
set loco7 lput mean [fmr_loco / fmr] of turtles with [species = 7 and age > 1] loco7
set repro7 lput mean [fmr_repro / fmr] of turtles with [species = 7 and age > 1] repro7
set grow7 lput mean [fmr_growth / fmr] of turtles with [species = 7 and age > 1] grow7
set basal7 lput mean [fmr_basal / fmr] of turtles with [species = 7 and age > 1] basal7
set digest7 lput mean [fmr_digest / fmr] of turtles with [species = 7 and age > 1] digest7
set prod7 lput mean [production] of turtles with [species = 7 and age > 1] prod7
set in7 lput mean [intake] of turtles with [species = 7 and age > 1] in7
set inpatch7 lput mean [intake / forage-patches] of turtles with [species = 7 and age > 1] inpatch7
set inhr7 lput mean [intake / (hrsize + 1)] of turtles with [species = 7 and age > 1] inhr7
set balance7 lput median [balance] of turtles with [species = 7 and age > 1] balance7
set prodbalance7 lput median [balanceprod] of turtles with [species = 7 and age > 1] prodbalance7
    if drought_scenario
    [
     if drought
     [
        set drought_hr7 lput mean [hrsize] of turtles with [species = 7 and age > 1] drought_hr7
        set drought_fmr7 lput mean [fmr] of turtles with [species = 7 and age > 1] drought_fmr7
        set drought_in7 lput mean [intake] of turtles with [species = 7 and age > 1] drought_in7
        set drought_stor7 lput mean [storage / max_storage] of turtles with [species = 7 and age > 1] drought_stor7
     ]
    ]
  ]


  if any? turtles with [species = 8 and age > 1]
  [
set hr8 lput mean [hrsize] of turtles with [species = 8 and age > 1] hr8
set maxhr8 lput mean [hrsize / maxhr] of turtles with [species = 8 and age > 1] maxhr8
set countmaxhr8 lput count turtles with [species = 8 and age > 1 and hrsize / maxhr > 0.9] countmaxhr8
set age8 lput mean [age] of turtles with [species = 8 and age > 1] age8
set stor8 lput mean [storage / max_storage] of turtles with [species = 8 and age > 1] stor8
set move8 lput mean [daymoved] of turtles with [species = 8 and age > 1] move8
set patches8 lput mean [forage-patches] of turtles with [species = 8 and age > 1] patches8
set compet8 lput mean [competition / forage-patches] of turtles with [species = 8 and age > 1] compet8
set order8 lput mean [order] of turtles with [species = 8 and age > 1] order8
set vis_before8 lput mean [vis_before / forage-patches] of turtles with [species = 8 and age > 1] vis_before8
set fmr8 lput mean [fmr] of turtles with [species = 8 and age > 1] fmr8
set loco8 lput mean [fmr_loco / fmr] of turtles with [species = 8 and age > 1] loco8
set repro8 lput mean [fmr_repro / fmr] of turtles with [species = 8 and age > 1] repro8
set grow8 lput mean [fmr_growth / fmr] of turtles with [species = 8 and age > 1] grow8
set basal8 lput mean [fmr_basal / fmr] of turtles with [species = 8 and age > 1] basal8
set digest8 lput mean [fmr_digest / fmr] of turtles with [species = 8 and age > 1] digest8
set prod8 lput mean [production] of turtles with [species = 8 and age > 1] prod8
set in8 lput mean [intake] of turtles with [species = 8 and age > 1] in8
set inpatch8 lput mean [intake / forage-patches] of turtles with [species = 8 and age > 1] inpatch8
set inhr8 lput mean [intake / (hrsize + 1)] of turtles with [species = 8 and age > 1] inhr8
set balance8 lput median [balance] of turtles with [species = 8 and age > 1] balance8
set prodbalance8 lput median [balanceprod] of turtles with [species = 8 and age > 1] prodbalance8
    if drought_scenario
    [
     if drought
     [
        set drought_hr8 lput mean [hrsize] of turtles with [species = 8 and age > 1] drought_hr8
        set drought_fmr8 lput mean [fmr] of turtles with [species = 8 and age > 1] drought_fmr8
        set drought_in8 lput mean [intake] of turtles with [species = 8 and age > 1] drought_in8
        set drought_stor8 lput mean [storage / max_storage] of turtles with [species = 8 and age > 1] drought_stor8
     ]
    ]
  ]

  if any? turtles with [species = 9 and age > 1]
  [
set hr9 lput mean [hrsize] of turtles with [species = 9 and age > 1] hr9
set maxhr9 lput mean [hrsize / maxhr] of turtles with [species = 9 and age > 1] maxhr9
set countmaxhr9 lput count turtles with [species = 9 and age > 1 and hrsize / maxhr > 0.9] countmaxhr9
set age9 lput mean [age] of turtles with [species = 9 and age > 1] age9
set stor9 lput mean [storage / max_storage] of turtles with [species = 9 and age > 1] stor9
set move9 lput mean [daymoved] of turtles with [species = 9 and age > 1] move9
set patches9 lput mean [forage-patches] of turtles with [species = 9 and age > 1] patches9
set compet9 lput mean [competition / forage-patches] of turtles with [species = 9 and age > 1] compet9
set order9 lput mean [order] of turtles with [species = 9 and age > 1] order9
set vis_before9 lput mean [vis_before / forage-patches] of turtles with [species = 9 and age > 1] vis_before9
set fmr9 lput mean [fmr] of turtles with [species = 9 and age > 1] fmr9
set loco9 lput mean [fmr_loco / fmr] of turtles with [species = 9 and age > 1] loco9
set repro9 lput mean [fmr_repro / fmr] of turtles with [species = 9 and age > 1] repro9
set grow9 lput mean [fmr_growth / fmr] of turtles with [species = 9 and age > 1] grow9
set basal9 lput mean [fmr_basal / fmr] of turtles with [species = 9 and age > 1] basal9
set digest9 lput mean [fmr_digest / fmr] of turtles with [species = 9 and age > 1] digest9
set prod9 lput mean [production] of turtles with [species = 9 and age > 1] prod9
set in9 lput mean [intake] of turtles with [species = 9 and age > 1] in9
set inpatch9 lput mean [intake / forage-patches] of turtles with [species = 9 and age > 1] inpatch9
set inhr9 lput mean [intake / (hrsize + 1)] of turtles with [species = 9 and age > 1] inhr9
set balance9 lput median [balance] of turtles with [species = 9 and age > 1] balance9
set prodbalance9 lput median [balanceprod] of turtles with [species = 9 and age > 1] prodbalance9
    if drought_scenario
    [
     if drought
     [
        set drought_hr9 lput mean [hrsize] of turtles with [species = 9 and age > 1] drought_hr9
        set drought_fmr9 lput mean [fmr] of turtles with [species = 9 and age > 1] drought_fmr9
        set drought_in9 lput mean [intake] of turtles with [species = 9 and age > 1] drought_in9
        set drought_stor9 lput mean [storage / max_storage] of turtles with [species = 9 and age > 1] drought_stor9
     ]
    ]
  ]
end
;-------------------------------------------------------------------------------------------------------------------

to-report nanmean [l]   ;;function calculate the mean if possible, otherwise reprot 999 without error
  ifelse length l > 0
  [report mean l]
  [report 999]
end
@#$#@#$#@
GRAPHICS-WINDOW
212
13
688
490
-1
-1
4.634
1
6
1
1
1
0
1
1
1
0
100
0
100
1
1
1
ticks
30.0

BUTTON
29
17
109
51
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
29
55
111
89
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

PLOT
1356
15
1609
219
Species richness
species ID
nr individuals
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 1 -12895429 true "histogram [species] of turtles ;with [hrsize > 0]" "histogram [species] of turtles ;with [hrsize > 0]"

TEXTBOX
426
483
593
503
1km
11
0.0
1

PLOT
698
15
1347
220
Species numbers
time
number of ind.
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"sp0" 1.0 0 -16777216 true "" "plot count turtles with [species = 0]"
"sp1" 1.0 0 -7500403 true "" "plot count turtles with [species = 1]"
"sp2" 1.0 0 -2674135 true "" "plot count turtles with [species = 2]"
"sp3" 1.0 0 -955883 true "" "plot count turtles with [species = 3]"
"sp4" 1.0 0 -6459832 true "" "plot count turtles with [species = 4]"
"sp5" 1.0 0 -1184463 true "" "plot count turtles with [species = 5]"
"sp6" 1.0 0 -10899396 true "" "plot count turtles with [species = 6]"
"sp7" 1.0 0 -13840069 true "" "plot count turtles with [species = 7]"
"sp8" 1.0 0 -14835848 true "" "plot count turtles with [species = 8]"
"sp9" 1.0 0 -11221820 true "" "plot count turtles with [species = 9]"

PLOT
1360
438
1609
643
Leftover patchfeed
time
patchfeed
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"patchfeed" 1.0 0 -16777216 true "plot sum [patchfeed] of patches" "plot sum [patchfeed] of patches"

SLIDER
30
278
202
311
size_small
size_small
1
100
1.0
1
1
NIL
HORIZONTAL

PLOT
699
227
1348
432
Mean energy storage
time
energy storage
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"sp0" 1.0 0 -16777216 true "" "plot mean [storage] of turtles with [species = 0]"
"sp1" 1.0 0 -7500403 true "" "plot mean [storage] of turtles with [species = 1]"
"sp2" 1.0 0 -2674135 true "" "plot mean [storage] of turtles with [species = 2]"
"sp3" 1.0 0 -955883 true "" "plot mean [storage] of turtles with [species = 3]"
"sp4" 1.0 0 -6459832 true "" "plot mean [storage] of turtles with [species = 4]"
"sp5" 1.0 0 -1184463 true "" "plot mean [storage] of turtles with [species = 5]"
"sp6" 1.0 0 -10899396 true "" "plot mean [storage] of turtles with [species = 6]"
"sp7" 1.0 0 -13840069 true "" "plot mean [storage] of turtles with [species = 7]"
"sp8" 1.0 0 -14835848 true "" "plot mean [storage] of turtles with [species = 8]"
"sp9" 1.0 0 -11221820 true "" "plot mean [storage] of turtles with [species = 9]"

CHOOSER
29
186
204
231
clump
clump
0.9 0.99 0.999 0.9999
1

CHOOSER
33
484
171
529
specs-included
specs-included
"all" 0 1 2 3 4 5 6 7 8 9
1

PLOT
700
438
1350
644
Mean field metabolic rate
time
FMR [g food]
0.0
12.0
0.0
12.0
true
true
"" ""
PENS
"sp0" 1.0 0 -16777216 true "" "plot mean [fmr] of turtles with [species = 0]"
"sp1" 1.0 0 -7500403 true "" "plot mean [fmr] of turtles with [species = 1]"
"sp2" 1.0 0 -2674135 true "" "plot mean [fmr] of turtles with [species = 2]"
"sp3" 1.0 0 -955883 true "" "plot mean [fmr] of turtles with [species = 3]"
"sp4" 1.0 0 -6459832 true "" "plot mean [fmr] of turtles with [species = 4]"
"sp5" 1.0 0 -1184463 true "" "plot mean [fmr] of turtles with [species = 5]"
"sp6" 1.0 0 -10899396 true "" "plot mean [fmr] of turtles with [species = 6]"
"sp7" 1.0 0 -13840069 true "" "plot mean [fmr] of turtles with [species = 7]"
"sp8" 1.0 0 -14835848 true "" "plot mean [fmr] of turtles with [species = 8]"
"sp9" 1.0 0 -11221820 true "" "plot mean [fmr] of turtles with [species = 9]"

SLIDER
31
317
203
350
feed-amount
feed-amount
0
50
16.0
1
1
NIL
HORIZONTAL

CHOOSER
29
133
203
178
fragmentation_mode
fragmentation_mode
"frag" "slass"
0

PLOT
1358
229
1612
431
Species number
time
# species
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot spec_num"

CHOOSER
522
700
667
745
debug
debug
0 1
0

INPUTBOX
233
528
301
588
digestion
0.5
1
0
Number

INPUTBOX
308
528
373
588
mortbold
3.33E-5
1
0
Number

INPUTBOX
382
528
444
588
hr_try_juv
10.0
1
0
Number

INPUTBOX
450
594
515
654
maxmass
0.1
1
0
Number

INPUTBOX
234
593
301
653
prio_fact
0.2
1
0
Number

INPUTBOX
308
594
373
654
move-factor
3.0
1
0
Number

INPUTBOX
382
594
444
654
bold_prob
1.0
1
0
Number

INPUTBOX
452
528
515
588
hr_try_init
100.0
1
0
Number

INPUTBOX
523
529
583
589
cover_percentage
0.05
1
0
Number

INPUTBOX
524
595
583
655
starting_indivs
1000.0
1
0
Number

TEXTBOX
29
113
179
131
Landscape Setup:
12
0.0
1

PLOT
700
653
1351
803
Individuals failing to find enough food
individuals
time
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"juveniles" 1.0 0 -16777216 true "" "plot fail"
"adults" 1.0 0 -7500403 true "" "plot fail3"

SWITCH
30
360
202
393
seasonal-resources
seasonal-resources
0
1
-1000

SLIDER
30
402
202
435
season-variability
season-variability
0
10
4.0
1
1
NIL
HORIZONTAL

INPUTBOX
599
529
666
589
maturity
0.95
1
0
Number

INPUTBOX
599
595
666
655
storage_addition
3.0
1
0
Number

SLIDER
31
757
173
790
drought_length
drought_length
0
100
1.0
1
1
NIL
HORIZONTAL

CHOOSER
32
613
170
658
drought_type
drought_type
"aprubt" "continous" "controlled" "data"
2

SLIDER
31
796
173
829
drought_min
drought_min
0
1
0.02
0.01
1
NIL
HORIZONTAL

SLIDER
32
717
174
750
drought_recover
drought_recover
0
100
1.0
1
1
NIL
HORIZONTAL

INPUTBOX
195
700
488
760
Drought_file
/home/leonna/Documents/Netlogo/Test1_drought2012.csv
1
0
String

CHOOSER
32
666
174
711
drought_scenario_combi
drought_scenario_combi
0 1 2 3 4 5 6 7 8
1

INPUTBOX
195
770
258
830
Location
0.0
1
0
Number

SLIDER
30
238
202
271
perc_small
perc_small
0
1
0.0
0.05
1
NIL
HORIZONTAL

TEXTBOX
33
463
221
486
Species Setup:
12
0.0
1

TEXTBOX
235
505
423
528
Parameters;
12
0.0
1

TEXTBOX
33
551
221
574
Drought scenario:
12
0.0
1

TEXTBOX
197
679
385
702
For droughts based on data:
12
0.0
1

SWITCH
32
572
197
605
drought_scenario
drought_scenario
0
1
-1000

TEXTBOX
524
677
674
695
Debugging:
12
0.0
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment_drought_data" repetitions="10" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="4745"/>
    <exitCondition>count turtles = 0</exitCondition>
    <metric>spec_num</metric>
    <metric>number0</metric>
    <metric>number1</metric>
    <metric>number2</metric>
    <metric>number3</metric>
    <metric>number4</metric>
    <metric>number5</metric>
    <metric>number6</metric>
    <metric>number7</metric>
    <metric>number8</metric>
    <metric>number9</metric>
    <enumeratedValueSet variable="perc_small">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="size_small">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bold_prob">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fragmentation_mode">
      <value value="&quot;frag&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="clump">
      <value value="0.9999"/>
      <value value="0.999"/>
      <value value="0.99"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="loss">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_scenario">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starting_indivs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="feed-amount">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prio_fact">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immi_prob">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seasonal-resources">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="season-variability">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_length">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_min">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_recover">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="specs-included">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_type">
      <value value="&quot;data&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Drought_file">
      <value value="&quot;Data1_drought1952_002.csv&quot;"/>
      <value value="&quot;Data1_drought1952_005.csv&quot;"/>
      <value value="&quot;Data1_drought1952_01.csv&quot;"/>
      <value value="&quot;Data1_drought1952_02.csv&quot;"/>
      <value value="&quot;Data1_drought2009_002.csv&quot;"/>
      <value value="&quot;Data1_drought2009_005.csv&quot;"/>
      <value value="&quot;Data1_drought2009_01.csv&quot;"/>
      <value value="&quot;Data1_drought2009_02.csv&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Location">
      <value value="0"/>
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
      <value value="6"/>
      <value value="7"/>
      <value value="8"/>
      <value value="9"/>
      <value value="10"/>
      <value value="11"/>
      <value value="12"/>
      <value value="13"/>
      <value value="14"/>
      <value value="15"/>
      <value value="16"/>
      <value value="17"/>
      <value value="18"/>
      <value value="19"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_drought_community" repetitions="20" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="3650"/>
    <exitCondition>count turtles = 0</exitCondition>
    <metric>spec_num</metric>
    <metric>number0</metric>
    <metric>number1</metric>
    <metric>number2</metric>
    <metric>number3</metric>
    <metric>number4</metric>
    <metric>number5</metric>
    <metric>number6</metric>
    <metric>number7</metric>
    <metric>number8</metric>
    <metric>number9</metric>
    <metric>fmr0</metric>
    <metric>fmr1</metric>
    <metric>fmr2</metric>
    <metric>fmr3</metric>
    <metric>fmr4</metric>
    <metric>fmr5</metric>
    <metric>fmr6</metric>
    <metric>fmr7</metric>
    <metric>fmr8</metric>
    <metric>fmr9</metric>
    <metric>in0</metric>
    <metric>in1</metric>
    <metric>in2</metric>
    <metric>in3</metric>
    <metric>in4</metric>
    <metric>in5</metric>
    <metric>in6</metric>
    <metric>in7</metric>
    <metric>in8</metric>
    <metric>in9</metric>
    <metric>balance0</metric>
    <metric>balance1</metric>
    <metric>balance2</metric>
    <metric>balance3</metric>
    <metric>balance4</metric>
    <metric>balance5</metric>
    <metric>balance6</metric>
    <metric>balance7</metric>
    <metric>balance8</metric>
    <metric>balance9</metric>
    <metric>repro0</metric>
    <metric>repro1</metric>
    <metric>repro2</metric>
    <metric>repro3</metric>
    <metric>repro4</metric>
    <metric>repro5</metric>
    <metric>repro6</metric>
    <metric>repro7</metric>
    <metric>repro8</metric>
    <metric>repro9</metric>
    <metric>drought_hr0</metric>
    <metric>drought_hr1</metric>
    <metric>drought_hr2</metric>
    <metric>drought_hr3</metric>
    <metric>drought_hr4</metric>
    <metric>drought_hr5</metric>
    <metric>drought_hr6</metric>
    <metric>drought_hr7</metric>
    <metric>drought_hr8</metric>
    <metric>drought_hr9</metric>
    <metric>drought_stor0</metric>
    <metric>drought_stor1</metric>
    <metric>drought_stor2</metric>
    <metric>drought_stor3</metric>
    <metric>drought_stor4</metric>
    <metric>drought_stor5</metric>
    <metric>drought_stor6</metric>
    <metric>drought_stor7</metric>
    <metric>drought_stor8</metric>
    <metric>drought_stor9</metric>
    <metric>sublist loco0 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco1 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco2 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco3 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco4 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco5 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco6 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco7 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco8 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco9 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet0 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet1 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet2 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet3 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet4 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet5 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet6 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet7 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet8 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet9 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>nanmean rep_success_0</metric>
    <metric>nanmean rep_success_1</metric>
    <metric>nanmean rep_success_2</metric>
    <metric>nanmean rep_success_3</metric>
    <metric>nanmean rep_success_4</metric>
    <metric>nanmean rep_success_5</metric>
    <metric>nanmean rep_success_6</metric>
    <metric>nanmean rep_success_7</metric>
    <metric>nanmean rep_success_8</metric>
    <metric>nanmean rep_success_9</metric>
    <metric>nanmean rep_success_hr_0</metric>
    <metric>nanmean rep_success_hr_1</metric>
    <metric>nanmean rep_success_hr_2</metric>
    <metric>nanmean rep_success_hr_3</metric>
    <metric>nanmean rep_success_hr_4</metric>
    <metric>nanmean rep_success_hr_5</metric>
    <metric>nanmean rep_success_hr_6</metric>
    <metric>nanmean rep_success_hr_7</metric>
    <metric>nanmean rep_success_hr_8</metric>
    <metric>nanmean rep_success_hr_9</metric>
    <enumeratedValueSet variable="bold_prob">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fragmentation_mode">
      <value value="&quot;frag&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="clump">
      <value value="0.9999"/>
      <value value="0.999"/>
      <value value="0.99"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starting_indivs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="feed-amount">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prio_fact">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seasonal-resources">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_scenario">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="season-variability">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_scenario_combi">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
      <value value="6"/>
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_recover">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="specs-included">
      <value value="&quot;all&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_type">
      <value value="&quot;controlled&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment_drought_singlespecies" repetitions="20" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="1095"/>
    <exitCondition>count turtles = 0</exitCondition>
    <metric>number0</metric>
    <metric>number1</metric>
    <metric>number2</metric>
    <metric>number3</metric>
    <metric>number4</metric>
    <metric>number5</metric>
    <metric>number6</metric>
    <metric>number7</metric>
    <metric>number8</metric>
    <metric>number9</metric>
    <metric>fmr0</metric>
    <metric>fmr1</metric>
    <metric>fmr2</metric>
    <metric>fmr3</metric>
    <metric>fmr4</metric>
    <metric>fmr5</metric>
    <metric>fmr6</metric>
    <metric>fmr7</metric>
    <metric>fmr8</metric>
    <metric>fmr9</metric>
    <metric>in0</metric>
    <metric>in1</metric>
    <metric>in2</metric>
    <metric>in3</metric>
    <metric>in4</metric>
    <metric>in5</metric>
    <metric>in6</metric>
    <metric>in7</metric>
    <metric>in8</metric>
    <metric>in9</metric>
    <metric>balance0</metric>
    <metric>balance1</metric>
    <metric>balance2</metric>
    <metric>balance3</metric>
    <metric>balance4</metric>
    <metric>balance5</metric>
    <metric>balance6</metric>
    <metric>balance7</metric>
    <metric>balance8</metric>
    <metric>balance9</metric>
    <metric>repro0</metric>
    <metric>repro1</metric>
    <metric>repro2</metric>
    <metric>repro3</metric>
    <metric>repro4</metric>
    <metric>repro5</metric>
    <metric>repro6</metric>
    <metric>repro7</metric>
    <metric>repro8</metric>
    <metric>repro9</metric>
    <metric>drought_hr0</metric>
    <metric>drought_hr1</metric>
    <metric>drought_hr2</metric>
    <metric>drought_hr3</metric>
    <metric>drought_hr4</metric>
    <metric>drought_hr5</metric>
    <metric>drought_hr6</metric>
    <metric>drought_hr7</metric>
    <metric>drought_hr8</metric>
    <metric>drought_hr9</metric>
    <metric>drought_stor0</metric>
    <metric>drought_stor1</metric>
    <metric>drought_stor2</metric>
    <metric>drought_stor3</metric>
    <metric>drought_stor4</metric>
    <metric>drought_stor5</metric>
    <metric>drought_stor6</metric>
    <metric>drought_stor7</metric>
    <metric>drought_stor8</metric>
    <metric>drought_stor9</metric>
    <metric>sublist loco0 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco1 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco2 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco3 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco4 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco5 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco6 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco7 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco8 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist loco9 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet0 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet1 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet2 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet3 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet4 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet5 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet6 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet7 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet8 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>sublist compet9 drought_start (drought_start + drought_length + drought_recover + drought_recover + 2)</metric>
    <metric>nanmean rep_success_0</metric>
    <metric>nanmean rep_success_1</metric>
    <metric>nanmean rep_success_2</metric>
    <metric>nanmean rep_success_3</metric>
    <metric>nanmean rep_success_4</metric>
    <metric>nanmean rep_success_5</metric>
    <metric>nanmean rep_success_6</metric>
    <metric>nanmean rep_success_7</metric>
    <metric>nanmean rep_success_8</metric>
    <metric>nanmean rep_success_9</metric>
    <metric>nanmean rep_success_hr_0</metric>
    <metric>nanmean rep_success_hr_1</metric>
    <metric>nanmean rep_success_hr_2</metric>
    <metric>nanmean rep_success_hr_3</metric>
    <metric>nanmean rep_success_hr_4</metric>
    <metric>nanmean rep_success_hr_5</metric>
    <metric>nanmean rep_success_hr_6</metric>
    <metric>nanmean rep_success_hr_7</metric>
    <metric>nanmean rep_success_hr_8</metric>
    <metric>nanmean rep_success_hr_9</metric>
    <enumeratedValueSet variable="bold_prob">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="fragmentation_mode">
      <value value="&quot;frag&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="clump">
      <value value="0.9999"/>
      <value value="0.999"/>
      <value value="0.99"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="starting_indivs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="feed-amount">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="prio_fact">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seasonal-resources">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="season-variability">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_scenario">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_scenario_combi">
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
      <value value="6"/>
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="specs-included">
      <value value="0"/>
      <value value="1"/>
      <value value="2"/>
      <value value="3"/>
      <value value="4"/>
      <value value="5"/>
      <value value="6"/>
      <value value="7"/>
      <value value="8"/>
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="drought_type">
      <value value="&quot;controlled&quot;"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
