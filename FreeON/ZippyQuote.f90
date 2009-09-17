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

    CHARACTER(LEN=*), DIMENSION(45), PARAMETER :: quotes = (/ &
      "Are we having fun yet?", &
      "I am demographically correct.", &
      "I just became one with my browser software.", &
      "Virtual reality isn't what it used to be.", &
      "I want a mega-meal in a mega-mall.", &
      "Adopt my lifestyle or I'll have to press charges.", &
      "If you can't say something nice say something surrealistic.", &
      "I'm afraid! I need something in a heavy cream sauce.", &
      "I can silence Joan Rivers with a single slice of Kraft cheese.", &
      "I am protected by the power of stain-reistant Scotchguard", &
      "I just accepted provolone into my life.", &
      "Frivolity is a stern taskmaster.", &
      "All life is a blur of Republicans and meat.", &
      "I'm Zippy the Pinhead and I'm totally committed to the festive mode.", &
      "I just felt a paradigm shift.", &
      "My boxer shorts just went on a rampage through a Long Island bowling alley.", &
      "Glazed donuts are the building blocks of the universe.", &
      "Nobody brings small problems into a laundromat.", &
      "Consciousness is vastly overrated.", &
      "I hope my sensitive female side is wearing sensible leather pumps.", &
      "Reality distorts my sense of television.", &
      "AIEEEEE!  I am having an UNDULATING EXPERIENCE!", &
      "Am I accompanied by a PARENT or GUARDIAN?", &
      "Am I in GRADUATE SCHOOL yet?", &
      "Are we live or on tape?", &
      "As a FAD follower my BEVERAGE choices are rich and fulfilling!", &
      "BELA LUGOSI is my co-pilot..", &
      "Can I have an IMPULSE ITEM instead?", &
      "I used to be a FUNDAMENTALIST but then I heard about the HIGH RADIATION LEVELS and bought an ENCYCLOPEDIA!!", &
      "Everywhere I look I see NEGATIVITY and ASPHALT...", &
      "Half a mind is a terrible thing to waste!", &
      "He is the MELBA-BEING...  the ANGEL CAKE... XEROX him...  XEROX him --", &
      "Hmmm...  an arrogant bouquet with a subtle suggestion of POLYVINYL CHLORIDE...", &
      "I can't decide which WRONG TURN to make first!!", &
      "I guess you guys got BIG MUSCLES from doing too much STUDYING!", &
      "I had pancake makeup for brunch!", &
      "I just heard the SEVENTIES were over!!  And I was just getting in touch with my LEISURE SUIT!!", &
      "I think I am an overnight sensation right now!!", &
      "I think I'd better go back to my DESK and toy with a few common MISAPPREHENSIONS...", &
      "I will SHAVE and buy JELL-O and bring my MARRIAGE MANUAL!!", &
      "Impudent..  Yet possessing a certain ALUMINUM SILICATE overbite....Needs REDDY-WHIP!!", &
      "If I had a Q-TIP, I could prevent th'collapse of NEGOTIATIONS!!", &
      "If this was a SWEDISH MOVIE, I'd take off your GO-GO BOOTS!!", &
      "Hand me a pair of leather pants and a CASIO keyboard -- I'm living for today!", &
      "I'm using my X-RAY VISION to obtain a rare glimpse of the INNER WORKINGS of this POTATO!!" &
    /)

    CALL MondoLogPlain("Zippy sez: "//TRIM(quotes(Random((/ 1, SIZE(quotes) /)))))

  END SUBROUTINE ZippySez
END MODULE ZippyQuote
