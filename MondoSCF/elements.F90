!===============================================================================
!
!       Droits de reproduction et de diffusion réservés. © 2000 CEA/CNRS.
!              (Laboratoire de Dynamique Moléculaire/IBS/DSV) 2000
!
!===============================================================================
!
!                Copyright © 2000 CEA/CNRS. All Rights Reserved.
!              (Laboratoire de Dynamique Moléculaire/IBS/DSV) 2000
!
!===============================================================================
!                              The Elements Module
!===============================================================================
!
! . Scalar parameter data:
!
!   NELEMENTS                      The number of elements.
!
! . Array parameter data:
!
!   MASS                           The element masses (amu).
!   RADII                          The element radii (Angstroms).
!   SYMBOL                         The element symbols.
!
! . Element names:
!
!   Hydrogen      Helium        Lithium       Beryllium      Boron
!   Carbon        Nitrogen      Oxygen        Fluorine       Neon
!   Sodium        Magnesium     Aluminium     Silicon        Phosphorus
!   Sulphur       Chlorine      Argon         Potassium      Calcium
!   Scandium      Titanium      Vanadium      Chromium       Manganese
!   Iron          Cobalt        Nickel        Copper         Zinc
!   Gallium       Germanium     Arsenic       Selenium       Bromine
!   Krypton       Rubidium      Strontium     Yttrium        Zirconium
!   Niobium       Molybdenum    Technetium    Ruthenium      Rhodium
!   Palladium     Silver        Cadmium       Indium         Tin
!   Antimony      Tellurium     Iodine        Xenon          Caesium
!   Barium        Lanthanum     Cerium        Praseodymium   Neodymium
!   Promethium    Samarium      Europium      Gadolinium     Terbium
!   Dysprosium    Holmium       Erbium        Thulium        Ytterbium
!   Lutetium      Hafnium       Tantalum      Tungsten       Rhenium
!   Osmium        Iridium       Platinum      Gold           Mercury
!   Thallium      Lead          Bismuth       Polonium       Astatine
!   Radon         Francium      Radium        Actinium       Thorium
!   Protactinium  Uranium       Neptunium     Plutonium      Americium
!   Curium        Berkelium     Californium   Einsteinium    Fermium
!
! . Notes:
!
!   All data values are parameters and so cannot be changed.
!
!   The element masses were obtained from the CRC handbook (see "constants.f90"
!   for the full reference) and they are the values for the elements taking into
!   account the observed isotopic abundancies. These values might not be the
!   best for certain types of calculation (e.g. molecular dynamics simulations
!   in which it might be preferable to use masses appropriate for a particular
!   isotope).
!
!===============================================================================
MODULE ELEMENTS

! . Module declarations.
USE DEFINITIONS, ONLY : DP

IMPLICIT NONE
PUBLIC

! . The number of elements.
INTEGER, PARAMETER :: NELEMENTS = 100

! . The element symbols.
CHARACTER ( LEN = 2 ), DIMENSION(1:NELEMENTS), PARAMETER :: SYMBOL = (/ &
           'H ', 'He', 'Li', 'Be', 'B ', 'C ', 'N ', 'O ', 'F ', 'Ne',  &
           'Na', 'Mg', 'Al', 'Si', 'P ', 'S ', 'Cl', 'Ar', 'K ', 'Ca',  &
           'Sc', 'Ti', 'V ', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',  &
           'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y ', 'Zr',  &
           'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn',  &
           'Sb', 'Te', 'I ', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',  &
           'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',  &
           'Lu', 'Hf', 'Ta', 'W ', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',  &
           'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th',  &
           'Pa', 'U ', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm'  /)

! . The element masses.
REAL ( KIND = DP ), DIMENSION(1:NELEMENTS), PARAMETER :: MASS = (/         &
            1.00794_DP,  4.0026_DP,   6.941_DP,    9.0122_DP,  10.811_DP,  &
           12.011_DP,   14.0067_DP,  15.9994_DP,  18.9984_DP,  20.1797_DP, &
           22.9898_DP,  24.305_DP,   26.9815_DP,  28.0855_DP,  30.9738_DP, &
           32.066_DP,   35.453_DP,   39.948_DP,   39.098_DP,   40.078_DP,  &
           44.9559_DP,  47.88_DP,    50.9415_DP,  51.9961_DP,  54.9381_DP, &
           55.847_DP,   58.9332_DP,  58.6934_DP,  63.546_DP,   65.39_DP,   &
           69.723_DP,   72.61_DP,    74.9216_DP,  78.96_DP,    79.904_DP,  &
           83.80_DP,    85.4678_DP,  87.62_DP,    88.9059_DP,  91.224_DP,  &
           92.9064_DP,  95.94_DP,    98.0_DP,    101.07_DP,   102.9055_DP, &
          106.42_DP,   107.8682_DP, 112.411_DP,  114.82_DP,   118.710_DP,  &
          121.757_DP,  127.60_DP,   126.9045_DP, 131.29_DP,   132.9054_DP, &
          137.327_DP,  138.9055_DP, 140.115_DP,  140.9077_DP, 144.24_DP,   &
          147.0_DP,    150.36_DP,   151.965_DP,  157.25_DP,   158.9253_DP, &
          162.50_DP,   164.9303_DP, 167.26_DP,   168.9342_DP, 173.04_DP,   &
          174.967_DP,  178.49_DP,   180.9479_DP, 183.85_DP,   186.207_DP,  &
          190.2_DP,    192.22_DP,   195.08_DP,   196.9665_DP, 200.59_DP,   &
          204.3833_DP, 207.2_DP,    208.9804_DP, 210.0_DP,    210.0_DP,    &
          222.0_DP,    223.0_DP,    226.0_DP,    227.0_DP,    232.0381_DP, &
          231.0_DP,    238.0289_DP, 237.0_DP,    242.0_DP,    243.0_DP,    &
          247.0_DP,    247.0_DP,    249.0_DP,    254.0_DP,    253.0_DP    /)

! . The element radii.
REAL ( KIND = DP ), DIMENSION(1:NELEMENTS), PARAMETER :: RADII = (/        &
                              0.23_DP, 0.20_DP, 0.68_DP, 0.35_DP, 0.83_DP, &
                              0.68_DP, 0.68_DP, 0.68_DP, 0.64_DP, 0.60_DP, &
                              0.97_DP, 1.10_DP, 1.35_DP, 1.20_DP, 1.05_DP, &
                              1.02_DP, 0.99_DP, 0.90_DP, 1.33_DP, 0.99_DP, &
                              1.44_DP, 1.47_DP, 1.33_DP, 1.35_DP, 1.35_DP, &
                              1.34_DP, 1.33_DP, 1.50_DP, 1.52_DP, 1.45_DP, &
                              1.22_DP, 1.17_DP, 1.21_DP, 1.22_DP, 1.21_DP, &
                              1.10_DP, 1.47_DP, 1.12_DP, 1.78_DP, 1.56_DP, &
                              1.48_DP, 1.47_DP, 1.35_DP, 1.40_DP, 1.45_DP, &
                              1.50_DP, 1.59_DP, 1.69_DP, 1.63_DP, 1.46_DP, &
                              1.46_DP, 1.47_DP, 1.40_DP, 1.30_DP, 1.67_DP, &
                              1.34_DP, 1.87_DP, 1.83_DP, 1.82_DP, 1.81_DP, &
                              1.80_DP, 1.80_DP, 1.99_DP, 1.79_DP, 1.76_DP, &
                              1.75_DP, 1.74_DP, 1.73_DP, 1.72_DP, 1.94_DP, &
                              1.72_DP, 1.57_DP, 1.43_DP, 1.37_DP, 1.35_DP, &
                              1.37_DP, 1.32_DP, 1.50_DP, 1.50_DP, 1.70_DP, &
                              1.55_DP, 1.54_DP, 1.54_DP, 1.68_DP, 1.60_DP, &
                              2.00_DP, 2.00_DP, 1.90_DP, 1.88_DP, 1.79_DP, &
                              1.61_DP, 1.58_DP, 1.55_DP, 1.53_DP, 1.51_DP, &
                              1.40_DP, 1.40_DP, 1.40_DP, 1.40_DP, 1.40_DP /)

END MODULE ELEMENTS
