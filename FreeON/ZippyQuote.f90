!------------------------------------------------------------------------------
!    This code is part of the MondoSCF suite of programs for linear scaling
!    electronic structure theory and ab initio molecular dynamics.
!
!    Copyright (2004). The Regents of the University of California. This
!    material was produced under U.S. Government contract W-7405-ENG-36
!    for Los Alamos National Laboratory, which is operated by the University
!    of California for the U.S. Department of Energy. The U.S. Government has
!    rights to use, reproduce, and distribute this software.  NEITHER THE
!    GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
!    OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
!
!    This program is free software; you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by the
!    Free Software Foundation; either version 2 of the License, or (at your
!    option) any later version. Accordingly, this program is distributed in
!    the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
!    PURPOSE. See the GNU General Public License at www.gnu.org for details.
!
!    While you may do as you like with this software, the GNU license requires
!    that you clearly mark derivative software.  In addition, you are encouraged
!    to return derivative works to the MondoSCF group for review, and possible
!    disemination in future releases.
!------------------------------------------------------------------------------
MODULE ZippyQuote
  USE Order
  USE ControlStructures

  IMPLICIT NONE

CONTAINS
  !===========================================================================================================================
  ! LETS NOT TAKE OURSELVES TOO SERIOUSLY ...
  !===========================================================================================================================
  SUBROUTINE ZippySez(C)
    TYPE(Controls)  :: C
    INTEGER         :: N,PU
    !------------------------------------------------------------------------------------------------------------------------!
    PU=6
    WRITE(PU,*)'Zippy sez:'
    N=Random((/1,45/))
    SELECT CASE(N)
    CASE (1)
      WRITE(PU,*)"Are we having fun yet?"
    CASE (2)
      WRITE(PU,*)"I am demographically correct."
    CASE (3)
      WRITE(PU,*)"I just became one with my browser software."
    CASE (4)
      WRITE(PU,*)"Virtual reality isn't what it used to be."
    CASE (5)
      WRITE(PU,*)"I want a mega-meal in a mega-mall."
    CASE (6)
      WRITE(PU,*)"Adopt my lifestyle or I'll have to press charges."
    CASE (7)
      WRITE(PU,*)"If you can't say something nice say something surrealistic."
    CASE (8)
      WRITE(PU,*)"I'm afraid! I need something in a heavy cream sauce."
    CASE (9)
      WRITE(PU,*)"I can silence Joan Rivers with a single slice of Kraft cheese."
    CASE (10)
      WRITE(PU,*)"I am protected by the power of stain-reistant Scotchguard"
    CASE (11)
      WRITE(PU,*)"I just accepted provolone into my life."
    CASE (12)
      WRITE(PU,*)"Frivolity is a stern taskmaster."
    CASE (13)
      WRITE(PU,*)"All life is a blur of Republicans and meat."
    CASE (14)
      WRITE(PU,*)"I'm Zippy the Pinhead and I'm totally committed to the festive mode."
    CASE (15)
      WRITE(PU,*)"I just felt a paradigm shift."
    CASE (16)
      WRITE(PU,*)"My boxer shorts just went on a rampage through a Long Island bowling alley."
    CASE (17)
      WRITE(PU,*)"Glazed donuts are the building blocks of the universe."
    CASE (18)
      WRITE(PU,*)"Nobody brings small problems into a laundromat."
    CASE (19)
      WRITE(PU,*)"Consciousness is vastly overrated."
    CASE (20)
      WRITE(PU,*)"I hope my sensitive female side is wearing sensible leather pumps."
    CASE (21)
      WRITE(PU,*)"Reality distorts my sense of television."
    CASE (22)
      WRITE(PU,*)"AIEEEEE!  I am having an UNDULATING EXPERIENCE!"
    CASE (23)
      WRITE(PU,*)"Am I accompanied by a PARENT or GUARDIAN?"
    CASE (24)
      WRITE(PU,*)"Am I in GRADUATE SCHOOL yet?"
    CASE (25)
      WRITE(PU,*)"Are we live or on tape?"
    CASE (26)
      WRITE(PU,*)"As a FAD follower my BEVERAGE choices are rich and fulfilling!"
    CASE (27)
      WRITE(PU,*)"BELA LUGOSI is my co-pilot.."
    CASE (28)
      WRITE(PU,*)"Can I have an IMPULSE ITEM instead?"
    CASE (29)
      WRITE(PU,*)"I used to be a FUNDAMENTALIST but then I heard about the HIGH RADIATION LEVELS and bought an ENCYCLOPEDIA!!"
    CASE (30)
      WRITE(PU,*)"Everywhere I look I see NEGATIVITY and ASPHALT..."
    CASE (31)
      WRITE(PU,*)"Half a mind is a terrible thing to waste!"
    CASE (32)
      WRITE(PU,*)"He is the MELBA-BEING...  the ANGEL CAKE... XEROX him...  XEROX him --"
    CASE (33)
      WRITE(PU,*)"Hmmm...  an arrogant bouquet with a subtle suggestion of POLYVINYL CHLORIDE..."
    CASE (34)
      WRITE(PU,*)"I can't decide which WRONG TURN to make first!!"
    CASE (35)
      WRITE(PU,*)"I guess you guys got BIG MUSCLES from doing too much STUDYING!"
    CASE (36)
      WRITE(PU,*)"I had pancake makeup for brunch!"
    CASE (37)
      WRITE(PU,*)"I just heard the SEVENTIES were over!!  And I was just getting in touch with my LEISURE SUIT!!"
    CASE (38)
      WRITE(PU,*)"I think I am an overnight sensation right now!!"
    CASE (39)
      WRITE(PU,*)"I think I'd better go back to my DESK and toy with a few common MISAPPREHENSIONS..."
    CASE (40)
      WRITE(PU,*)"I will SHAVE and buy JELL-O and bring my MARRIAGE MANUAL!!"
    CASE (41)
      WRITE(PU,*)"Impudent..  Yet possessing a certain ALUMINUM SILICATE overbite....Needs REDDY-WHIP!!"
    CASE (42)
      WRITE(PU,*)"If I had a Q-TIP, I could prevent th'collapse of NEGOTIATIONS!!"
    CASE (43)
      WRITE(PU,*)"If this was a SWEDISH MOVIE, I'd take off your GO-GO BOOTS!!"
    CASE (44)
      WRITE(PU,*)"Hand me a pair of leather pants and a CASIO keyboard -- I'm living for today!"
    CASE (45)
      WRITE(PU,*)"I'm using my X-RAY VISION to obtain a rare glimpse of the INNER WORKINGS of this POTATO!!"
    END SELECT
  END SUBROUTINE ZippySez
END MODULE ZippyQuote
