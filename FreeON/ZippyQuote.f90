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

! The following Zippy quotes are taken from the "Vim global plugin for Zippy The Pinhead quotes".

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

    CHARACTER(LEN=*), DIMENSION(807), PARAMETER :: quotes = (/ &
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
      "I'm using my X-RAY VISION to obtain a rare glimpse of the INNER WORKINGS of this POTATO!!", &
      "A can of ASPARAGUS, 73 pigeons, some LIVE ammo, and a FROZEN DAQUIRI!!", &
      "A dwarf is passing out somewhere in Detroit!", &
      "A GRAM??  A BRAM...  A GROOM...  A BROOM...  Oh, Yeh!!  Wash the ROOM!!", &
      "...A housewife is wearing a polypyrene jumpsuit!!", &
      "A shapely CATHOLIC SCHOOLGIRL is FIDGETING inside my costume..", &
      "A wide-eyed, innocent UNICORN, poised delicately in a MEADOW filled with LILACS, LOLLIPOPS & small CHILDREN at the HUSH of twilight??", &
      "Actually, what I'd like is a little toy spaceship!!", &
      "After this, I'm going to BURN some RUBBER!!", &
      "After THIS, let's go to PHILADELPHIA and have TRIPLETS!!", &
      "AIEEEEE!  I am having an UNDULATING EXPERIENCE!", &
      "ALFRED JARRY!  Say something about th' DEATH of DISCO!!", &
      "All I can think of is a platter of organic PRUNE CRISPS being trampled by an army of swarthy, Italian LOUNGE SINGERS...", &
      "All of a sudden, I want to THROW OVER my promising ACTING CAREER, grow a LONG BLACK BEARD and wear a BASEBALL HAT!!  ...  Although I don't know WHY!!", &
      "All of life is a blur of Republicans and meat!", &
      "All right, you degenerates!  I want this place evacuated in 20 seconds!", &
      "All this time I've been VIEWING a RUSSIAN MIDGET SODOMIZE a HOUSECAT!", &
      "Alright, you!!  Imitate a WOUNDED SEAL pleading for a PARKING SPACE!!", &
      "Am I accompanied by a PARENT or GUARDIAN?", &
      "Am I elected yet?", &
      "..Am I in a SOAP OPERA??", &
      "Am I in GRADUATE SCHOOL yet?", &
      "Am I SHOPLIFTING?", &
      "America!!  I saw it all!!  Vomiting!  Waving!  JERRY FALWELLING into your void tube of UHF oblivion!!  SAFEWAY of the mind --", &
      "An air of FRENCH FRIES permeates my nostrils!!", &
      "An INK-LING?  Sure -- TAKE one!!  Did you BUY any COMMUNIST UNIFORMS??", &
      "An Italian is COMBING his hair in suburban DES MOINES!", &
      "And furthermore, my bowling average is unimpeachable!!!", &
      "ANN JILLIAN'S HAIR makes LONI ANDERSON'S HAIR look like RICARDO MONTALBAN'S HAIR!", &
      "Are BOTH T.V.S on??", &
      "..  are the STEWED PRUNES still in the HAIR DRYER?", &
      "..Are we having FUN yet...?", &
      "Are we live or on tape?", &
      "Are we on STRIKE yet?", &
      "Are we THERE yet?", &
      "Are we THERE yet?  My MIND is a SUBMARINE!!", &
      "Are we THERE yet?!", &
      "Are we THERE yet??", &
      "Are you guys lined up for the METHADONE PROGRAM or FOOD STAMPS??", &
      "Are you mentally here at Pizza Hut??", &
      "Are you selling NYLON OIL WELLS??  If so, we can use TWO DOZEN!!", &
      "Are you still an ALCOHOLIC?", &
      "Are you still SEXUALLY ACTIVE?  Did you BRING th' REINFORCEMENTS?", &
      "As a FAD follower, my BEVERAGE choices are rich and fulfilling!", &
      "As President I have to go vacuum my coin collection!", &
      "Ask me the DIFFERENCE between PHIL SoILVERS and ALEXANDER HAIG!!", &
      "Awright, which one of you hid my PENIS ENVY?", &
      "Bagels...", &
      "BARBARA STANWYCK makes me nervous!!", &
      "Barbie says, Take quaaludes in gin and go to a disco right away! But Ken says, WOO-WOO!!  No credit at ``Mr. Liquor''!!", &
      "BARRY..  That was the most HEART-WARMING rendition of ``I DID IT MY WAY'' I've ever heard!!", &
      "BEEP-BEEP!!  I'm a '49 STUDEBAKER!!", &
      "Being a BALD HERO is almost as FESTIVE as a TATTOOED KNOCKWURST.", &
      "BELA LUGOSI is my co-pilot..", &
      "BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-BI-", &
      "...  Blame it on the BOSSA NOVA!!!", &
      "..  bleakness....  desolation....  plastic forks...", &
      "Bo Derek ruined my life!", &
      "Boy, am I glad it's only 1971...", &
      "Boys, you have ALL been selected to LEAVE th' PLANET in 15 minutes!!", &
      "BRILL CREAM is CREAM O' WHEAT in another DIMENSION..", &
      "But was he mature enough last night at the lesbian masquerade?", &
      "By MEER biz doo SCHOIN..", &
      "C'MON, everybody!!  I've flown in LESLIE GORE and two dozen KOSHER BUTCHERS!  They'll be doing intricate MILITARY MANEUVERS to the soundtrack from 'OKLAHOMA'!!", &
      "CALIFORNIA is where people from IOWA or NEW YORK go to subscribe to CABLE TELEVISION!!", &
      "Can I have an IMPULSE ITEM instead?", &
      "Can you MAIL a BEAN CAKE?", &
      "Catsup and Mustard all over the place!  It's the Human Hamburger!", &
      "CHUBBY CHECKER just had a CHICKEN SANDWICH in downtown DULUTH!", &
      "CHUBBY CHECKER owns my BUILDING!", &
      "Civilization is fun!  Anyway, it keeps me busy!!", &
      "Clear the laundromat!!  This whirl-o-matic just had a nuclear meltdown!!", &
      "Concentrate on th'cute, li'l CARTOON GUYS! Remember the SERIAL NUMBERS!!  Follow the WHIPPLE AVE EXIT!! Have a FREE PEPSI!!  Turn LEFT at th'HOLIDAY INN!! JOIN the CREDIT WORLD!!  MAKE me an OFFER!!!", &
      "CONGRATULATIONS!  Now should I make thinly veiled comments about DIGNITY, self-esteem and finding TRUE FUN in your RIGHT VENTRICLE??", &
      "Content:  80% POLYESTER, 20% DACRON..  The waitress's UNIFORM sheds TARTAR SAUCE like an 8'' by 10'' GLOSSY..", &
      "Could I have a drug overdose?", &
      "'DARK SHADOWS' is on!!  Hey, I think the VAMPIRE forgot his UMBRELLA!!", &
      "Darling, my ELBOW is FLYING over FRANKFURT, Germany..", &
      "Dehydrated EGGS are STREWN across ROULETTE TABLES..", &
      "Did an Italian CRANE OPERATOR just experience uninhibited sensations in a MALIBU HOT TUB?", &
      "Did I do an INCORRECT THING??", &
      "Did I say I was a sardine?  Or a bus???", &
      "Did I SELL OUT yet??", &
      "Did we bring enough BEEF JERKY?", &
      "Did YOU find a DIGITAL WATCH in YOUR box of VELVEETA?", &
      "Did you find a DIGITAL WATCH in YOUR box of VELVEETA??", &
      "Did you GAIN WEIGHT in th' past 5 MINUTES or am I just DREAMING of two BROCCOLI FLORETS lying in an empty GAS TANK?", &
      "Did you move a lot of KOREAN STEAK KNIVES this trip, Dingy?", &
      "DIDI...  is that a MARTIAN name, or, are we in ISRAEL?", &
      "Didn't I buy a 1951 Packard from you last March in Cairo?", &
      "Didn't KIRKEGAARD wear out his TIRES in VIENNA during a SNOWSTORM of FREUD's unpaid DENTAL BILLS?", &
      "Disco oil bussing will create a throbbing naugahide pipeline running straight to the tropics from the rug producing regions and devalue the dollar!", &
      "Dizzy, are we 'REAL PEOPLE' or 'AMAZING ANIMALS'?", &
      "Do I have a lifestyle yet?", &
      "Do I hear th' SPINNING of various WHIRRING, ROUND, and WARM WHIRLOMATICS?!", &
      "Do you guys know we just passed thru a BLACK HOLE in space?", &
      "Do you have exactly what I want in a plaid poindexter bar bat??", &
      "..  Do you like ``TENDER VITTLES?''?", &
      "Do you need any MOUTH-TO-MOUTH resuscitation?", &
      "Do you think the ``Monkees'' should get gas on odd or even days?", &
      "Does someone from PEORIA have a SHORTER ATTENTION span than me?", &
      "Does that mean I'm not a well-adjusted person??", &
      "..  does your DRESSING ROOM have enough ASPARAGUS?", &
      "DON'T go!!  I'm not HOWARD COSELL!!  I know POLISH JOKES...  WAIT!! Don't go!!  I AM Howard Cosell!...  And I DON'T know Polish jokes!!", &
      "Don't hit me!!  I'm in the Twilight Zone!!!", &
      "Don't SANFORIZE me!!", &
      "Don't worry, nobody really LISTENS to lectures in MOSCOW, either! ..  FRENCH, HISTORY, ADVANCED CALCULUS, COMPUTER PROGRAMMING, BLACK STUDIES, SOCIOBIOLOGY!..  Are there any QUESTIONS??", &
      "Edwin Meese made me wear CORDOVANS!!", &
      "Eisenhower!!  Your mimeograph machine upsets my stomach!!", &
      "Either CONFESS now or we go to ``PEOPLE'S COURT''!!", &
      "Everybody gets free BORSCHT!", &
      "Everybody is going somewhere!!  It's probably a garage sale or a disaster Movie!!", &
      "..Everything is....FLIPPING AROUND!!", &
      "Everything will be ALL RIGHT if we can just remember things about ALGEBRA.. or SOCCER..  or SOCIALISM..", &
      "Everywhere I look I see NEGATIVITY and ASPHALT...", &
      "Excuse me, but didn't I tell you there's NO HOPE for the survival of OFFSET PRINTING?", &
      "FEELINGS are cascading over me!!!", &
      "Fell th' WHIRLING BUFFERS buffing away all that stress... Years of ROAD TAR gently washing away...", &
      "Finally, Zippy drives his 1958 RAMBLER METROPOLITAN into the faculty dining room.", &
      "FIRST, I was in a TRUCK...THEN, I was in a DINER...", &
      "FIRST, I'm covering you with OLIVE OIL and PRUNE WHIP!!", &
      "First, I'm going to give you all the ANSWERS to today's test.. So just plug in your SONY WALKMANS and relax!!", &
      "FISH-NET-FISH-NET-FISH-NET-FISH-NET-FISH!!", &
      "Fold, fold, FOLD!!  FOLDING many items!!", &
      "FOOLED you!  Absorb EGO SHATTERING impulse rays, polyester poltroon!!", &
      "Four thousand different MAGNATES, MOGULS & NABOBS are romping in my gothic solarium!!", &
      "FROZEN ENTREES may be flung by members of opposing SWANSON SECTS..", &
      "FUN is never having to say you're SUSHI!!", &
      "Gee, I feel kind of LIGHT in the head now, knowing I can't make my satellite dish PAYMENTS!", &
      "...Get me a GIN and TONIC!!...make it HAIR TONIC!!", &
      "Gibble, Gobble, we ACCEPT YOU ---", &
      "Give them RADAR-GUIDED SKEE-BALL LANES and VELVEETA BURRITOS!!", &
      "Go on, EMOTE!  I was RAISED on thought balloons!!", &
      "GOOD-NIGHT, everybody..  Now I have to go administer FIRST-AID to my pet LEISURE SUIT!!", &
      "Ha ha   Ha ha  Ha ha   Ha  Ha  Ha  Ha  -- When will I EVER stop HAVING FUN?!!", &
      "HAIR TONICS, please!!", &
      "Half a mind is a terrible thing to waste!", &
      "Hand me a pair of leather pants and a CASIO keyboard -- I'm living for today!", &
      "Has everybody got HALVAH spread all over their ANKLES??... Now, it's time to ``HAVE A NAGEELA''!!", &
      "Have my two-tone, 1958 Nash METRO brought around..", &
      "..  he dominates the DECADENT SUBWAY SCENE.", &
      "He is the MELBA-BEING...  the ANGEL CAKE...  XEROX him...  XEROX him --", &
      "He probably just wants to take over my CELLS and then EXPLODE inside me like a BARREL of runny CHOPPED LIVER!  Or maybe he'd like to PSYCHOLOGICALLY TERRORIZE ME until I have no objection to a RIGHT-WING MILITARY TAKEOVER of my apartment!!  I guess I should call AL PACINO!", &
      "HELLO KITTY gang terrorizes town, family STICKERED to death!", &
      "HELLO, everybody, I'm a HUMAN!!", &
      "Hello, GORRY-O!!  I'm a GENIUS from HARVARD!!", &
      "HELLO, little boys!   Gimme a MINT TULIP!!  Let's do the BOSSA NOVA!!", &
      "--Hello, POLICE?  I've got ABBOTT & COSTELLO here on suspicion of HIGHWAY ROBBERY!!", &
      "Hello.  I know the divorce rate among unmarried Catholic Alaskan females!!", &
      "Hello.  Just walk along and try NOT to think about your INTESTINES being almost FORTY YARDS LONG!!", &
      "Hello...  IRON CURTAIN?  Send over a SAUSAGE PIZZA! World War III?  No thanks!", &
      "Hello?  Enema Bondage?  I'm calling because I want to be happy, I guess..", &
      "Here I am at the flea market but nobody is buying my urine sample bottles..", &
      "..  here I am in 53 B.C. and all I want is a dill pickle!!", &
      "Here I am in the POSTERIOR OLFACTORY LOBULE but I don't see CARL SAGAN anywhere!!", &
      "Here is my refrigerator full of FLANK STEAK...and over there is my UPHOLSTERED CANOE...I don't know WHY I OWN them!!", &
      "Here we are in America...  when do we collect unemployment?", &
      "HERE!!  Put THIS on!!  I'm in CHARGE!!", &
      "Hey!!  Let's watch the' ELEVATOR go UP and DOWN at th' HILTON HOTEL!!", &
      "Hey, I LIKE that POINT!!", &
      "Hey, LOOK!!  A pair of SIZE 9 CAPRI PANTS!!  They probably belong to SAMMY DAVIS, JR.!!", &
      "Hey, wait a minute!!  I want a divorce!!..  you're not Clint Eastwood!!", &
      "Hey, waiter!  I want a NEW SHIRT and a PONY TAIL with lemon sauce!", &
      "Hiccuping & trembling into the WASTE DUMPS of New Jersey like some drunken CABBAGE PATCH DOLL, coughing in line at FIORUCCI'S!!", &
      "Hmmm..  a CRIPPLED ACCOUNTANT with a FALAFEL sandwich is HIT by a TROLLEY-CAR..", &
      "Hmmm..  A hash-singer and a cross-eyed guy were SLEEPING on a deserted  island, when...", &
      "Hmmm...  a PINHEAD, during an EARTHQUAKE, encounters an ALL-MIDGET FIDDLE ORCHESTRA...  ha..  ha..", &
      "Hmmm...  an arrogant bouquet with a subtle suggestion of POLYVINYL CHLORIDE...", &
      "Hold the MAYO & pass the COSMIC AWARENESS...", &
      "HOORAY, Ronald!!  Now YOU can marry LINDA RONSTADT too!!", &
      "HOW could a GLASS be YELLING??", &
      "How do I get HOME?", &
      "How do you explain Wayne Newton's POWER over millions? It's th' MOUSTACHE...  Have you ever noticed th' way it radiates SINCERITY, HONESTY & WARMTH?  It's a MOUSTACHE you want to take HOME and introduce to NANCY SINATRA!", &
      "How many retired bricklayers from FLORIDA are out purchasing PENCIL SHARPENERS right NOW??", &
      "How's it going in those MODULAR LOVE UNITS??", &
      "How's the wife?  Is she at home enjoying capitalism?", &
      "..  hubub, hubub, HUBUB, hubub, hubub, hubub, HUBUB, hubub, hubub, hubub.", &
      "HUGH BEAUMONT died in 1982!!", &
      "HUMAN REPLICAS are inserted into VATS of NUTRITIONAL YEAST...", &
      "Hydraulic pizza oven!!  Guided missile!  Herring sandwich!  Styrofoam! Jayne Mansfield!  Aluminum siding!  Borax!  Pedal pushers!  Jukebox!", &
      "I always have fun because I'm out of my mind!!!", &
      "I always liked FLAG DAY!!", &
      "I always wanted a NOSE JOB!!", &
      "I am a jelly donut.  I am a jelly donut.", &
      "I am a traffic light, and Alan Ginzberg kidnapped my laundry in 1927!", &
      "I am covered with pure vegetable oil and I am writing a best seller!", &
      "I am deeply CONCERNED and I want something GOOD for BREAKFAST!", &
      "I am having a CONCEPTION--", &
      "I am having a pleasant time!!", &
      "I am having FUN...  I wonder if it's NET FUN or GROSS FUN?", &
      "I am KING BOMBA of Sicily!..I will marry LUCILLE BALL next Friday!", &
      "I am NOT a nut....", &
      "I appoint you ambassador to Fantasy Island!!!", &
      "I brought my BOWLING BALL - and some DRUGS!!", &
      "I call it a 'SARDINE ON WHEAT'!", &
      "-- I can do ANYTHING ... I can even ... SHOPLIFT!!", &
      "I can see you GUYS an' GALS need a LOT of HELP...You're all very STUPID!!  I used to be STUPID, too..before I started watching UHF-TV!!", &
      "I can't decide which WRONG TURN to make first!! I wonder if BOB GUCCIONE has these problems!", &
      "I can't think about that.  It doesn't go with HEDGES in the shape of LITTLE LULU -- or ROBOTS making BRICKS...", &
      "I decided to be JOHN TRAVOLTA instead!!", &
      "I demand IMPUNITY!", &
      "I didn't order any WOO-WOO...  Maybe a YUBBA..  But no WOO-WOO!", &
      "I don't believe there really IS a GAS SHORTAGE..  I think it's all just a BIG HOAX on the part of the plastic sign salesmen-- ..  to sell more numbers!!", &
      "..  I don't know why but, suddenly, I want to discuss declining I.Q. LEVELS with a blue ribbon SENATE SUB-COMMITTEE!", &
      "I don't know WHY I said that..  I think it came from the FILLINGS in my rear molars..", &
      "I don't think you fellows would do so much RAPING and PILLAGING if you played more PINBALL and watched CABLE TELEVISION!!", &
      "..  I don't understand the HUMOR of the THREE STOOGES!!", &
      "I feel better about world problems now!", &
      "I feel like a wet parking meter on Darvon!", &
      "I feel like I am sharing a ``CORN-DOG'' with NIKITA KHRUSCHEV..", &
      "I feel like I'm in a Toilet Bowl with a thumbtack in my forehead!!", &
      "I feel partially hydrogenated!", &
      "I feel real SOPHISTICATED being in FRANCE!", &
      "..  I feel..  JUGULAR..", &
      "I fill MY industrial waste containers with old copies of the ``WATCHTOWER'' and then add HAWAIIAN PUNCH to the top..  They look NICE in the yard--", &
      "I FORGOT to do the DISHES!!", &
      "I guess it was all a DREAM..  or an episode of HAWAII FIVE-O...", &
      "I guess we can live on his POT FARM in HADES!!", &
      "I guess you guys got BIG MUSCLES from doing too much STUDYING!", &
      "I had a lease on an OEDIPUS COMPLEX back in '81...", &
      "I had pancake makeup for brunch!", &
      "I have a TINY BOWL in my HEAD", &
      "I HAVE a towel.", &
      "I have a very good DENTAL PLAN.  Thank you.", &
      "..  I have a VISION!  It's a RANCID double-FISHWICH on an ENRICHED BUN!!", &
      "I have accepted Provolone into my life!", &
      "I have many CHARTS and DIAGRAMS..", &
      "I have no actual hairline...", &
      "I have nostalgia for the late Sixties!  In 1969 I left my laundry with a hippie!!  During an unauthorized tupperware party it was chopped & diced!", &
      "--- I have seen the FUN ---", &
      "I have seen these EGG EXTENDERS in my Supermarket.. ..  I have read the INSTRUCTIONS...", &
      "I have the power to HALT PRODUCTION on all TEENAGE SEX COMEDIES!!", &
      "I HAVE to buy a new ``DODGE MISER'' and two dozen JORDACHE JEANS because my viewscreen is ``USER-FRIENDLY''!!", &
      "I haven't been married in over six years, but we had sexual counseling every day from Oral Roberts!!", &
      "I HIJACKED a 747 to get here!!  I hope those fabulous CONEHEADS are at HOME!!", &
      "I hope I bought the right relish...  zzzzzzzzz...", &
      "I hope something GOOD came in the mail today so I have a REASON to live!!", &
      "I hope the ``Eurythmics'' practice birth control...", &
      "I hope you millionaires are having fun!  I just invested half your life savings in yeast!!", &
      "I invented skydiving in 1989!", &
      "I joined scientology at a garage sale!!", &
      "I just bought FLATBUSH from MICKEY MANTLE!", &
      "I just forgot my whole philosophy of life!!!", &
      "I just got my PRINCE bumper sticker.. But now I can't remember WHO he is...", &
      "I just had a MAJOR CONTRACT DISPUTE with SUZANNE SOMERS!!", &
      "I just had a NOSE JOB!!", &
      "I just had my entire INTESTINAL TRACT coated with TEFLON!", &
      "I just heard the SEVENTIES were over!!  And I was just getting in touch with my LEISURE SUIT!!", &
      "I just put lots of the EGG SALAD in the SILK SOCKS --", &
      "I just remembered something about a TOAD!", &
      "..I just walked into th' HOUSE OF REPRESENTATIVES with fourteen WET DOLPHINS and an out-of-date MARRIAGE MANUAL...", &
      "I KAISER ROLL?!  What good is a Kaiser Roll without a little COLE SLAW on the SIDE?", &
      "I Know A Joke", &
      "I know how to do SPECIAL EFFECTS!!", &
      "I know how to get the hostesses released!  Give them their own television series!", &
      "I know th'MAMBO!!  I have a TWO-TONE CHEMISTRY SET!!", &
      "I know things about TROY DONAHUE that can't even be PRINTED!!", &
      "I left my WALLET in the BATHROOM!!", &
      "I LIKE Aisle 7a.", &
      "I like the IMPUDENT NOSE on that car..  Are you a TEEN-AGER?  ", &
      "I like the way ONLY their mouths move..  They look like DYING OYSTERS", &
      "I like your SNOOPY POSTER!!", &
      "... I live in a FUR-LINE FALLOUT SHELTER", &
      "I love FRUIT PICKERS!!", &
      "--``I love KATRINKA because she drives a PONTIAC.  We're going away now.  I fed the cat. - Zippy''", &
      "I love ROCK 'N ROLL!  I memorized the all WORDS to ``WIPE-OUT'' in 1965!!", &
      "..I must be a VETERINARIAN..", &
      "I need 'RONDO'.", &
      "I need to discuss BUY-BACK PROVISIONS with at least six studio SLEAZEBALLS!!", &
      "I once decorated my apartment entirely in ten foot salad forks!!", &
      "I own seven-eighths of all the artists in downtown Burbank!", &
      "I OWN six pink HIPPOS!!", &
      "I predict that by 1993 everyone will live in and around LAS VEGAS and wear BEATLE HAIRCUTS!", &
      "I pretend I'm living in a styrofoam packing crate, high in th' SWISS ALPS, still unable to accept th' idea of TOUCH-TONE DIALING!!", &
      "I put aside my copy of ``BOWLING WORLD'' and think about GUN CONTROL legislation..", &
      "I represent a sardine!!", &
      "I request a weekend in Havana with Phil Silvers!", &
      "..  I see TOILET SEATS...", &
      "I selected E5...  but I didn't hear ``Sam the Sham and the Pharaohs''!", &
      "I smell a RANCID CORN DOG!", &
      "I smell like a wet reducing clinic on Columbus Day!", &
      "I think I am an overnight sensation right now!!", &
      "..  I think I'd better go back to my DESK and toy with a few common MISAPPREHENSIONS...", &
      "I think I'll do BOTH if I can get RESIDUALS!!", &
      "..  I think I'll KILL myself by leaping out of this 14th STOREY WINDOW while reading ERICA JONG'S poetry!!", &
      "I think I'll make SCRAMBLED EGGS!!  They're each in LITTLE SHELLS..", &
      "...I think I'm having an overnight sensation right now!!", &
      "I think my career is ruined!", &
      "I think my CAREER is RUINED!!", &
      "I used to be a FUNDAMENTALIST, but then I heard about the HIGH RADIATION LEVELS and bought an ENCYCLOPEDIA!!", &
      "..  I want a COLOR T.V. and a VIBRATING BED!!!", &
      "I want a VEGETARIAN BURRITO to go..  with EXTRA MSG!!", &
      "I want a WESSON OIL lease!!", &
      "I want another RE-WRITE on my CAESAR SALAD!!", &
      "I want DUSTIN HOFFMAN!! .. I want LIBRACE!!  YOW!!", &
      "I want EARS!  I want two ROUND BLACK EARS to make me feel warm 'n secure!!", &
      "..  I want FORTY-TWO TRYNEL FLOATATION SYSTEMS installed within SIX AND A HALF HOURS!!!", &
      "I want the presidency so bad I can already taste the hors d'oeuvres.", &
      "I want to dress you up as TALLULAH BANKHEAD and cover you with VASELINE and WHEAT THINS..", &
      "I want to kill everyone here with a cute colorful Hydrogen Bomb!!", &
      "..  I want to perform cranial activities with Tuesday Weld!!", &
      "I want to read my new poem about pork brains and outer space...", &
      "I want to so HAPPY, the VEINS in my neck STAND OUT!!", &
      "I want to TAKE IT HOME and DRESS IT UP in HOT PANTS!!", &
      "I want you to MEMORIZE the collected poems of EDNA ST VINCENT MILLAY.. BACKWARDS!!", &
      "I want you to organize my PASTRY trays...  my TEA-TINS are gleaming in formation like a ROW of DRUM MAJORETTES -- please don't be FURIOUS with me --", &
      "I was born in a Hostess Cupcake factory before the sexual revolution!", &
      "I was giving HAIR CUTS to th' SAUCER PEOPLE ..  I'm CLEAN!!", &
      "I was in a HOT TUB!  I was NORMAL!  I was ITALIAN!!  I enjoyed th' EARTHQUAKE!", &
      "I was in EXCRUCIATING PAIN until I started reading JACK AND JILL Magazine!!", &
      "I was making donuts and now I'm on a bus!", &
      "I will establish the first SHOPPING MALL in NUTLEY, New Jersey...", &
      "I will invent 'TIDY BOWL'...", &
      "I will SHAVE and buy JELL-O and bring my MARRIAGE MANUAL!!", &
      "I wish I was a sex-starved manicurist found dead in the Bronx!!", &
      "I wish I was on a Cincinnati street corner holding a clean dog!", &
      "I wonder if I could ever get started in the credit world?", &
      "..  I wonder if I ought to tell them about my PREVIOUS LIFE as a COMPLETE STRANGER?", &
      "I wonder if I should put myself in ESCROW!!", &
      "I wonder if there's anything GOOD on tonight?", &
      "I would like to urinate in an OVULAR, porcelain pool --", &
      "I'd like MY data-base JULIENNED and stir-fried!", &
      "I'd like some JUNK FOOD...  and then I want to be ALONE --", &
      "I'd like TRAINED SEALS and a CONVERTIBLE on my doorstep by NOON!!", &
      "I'll clean your ROOM!!  I know some GOOD stories, too!!  All about ROAD Island's, HUSH Puppies, and how LUKE finds GOLD on his LAND!!", &
      "I'll eat ANYTHING that's BRIGHT BLUE!!", &
      "I'LL get it!!  It's probably a FEW of my ITALIAN GIRL-FRIENDS!!", &
      "..I'll make you an ASHTRAY!!", &
      "I'll show you MY telex number if you show me YOURS...", &
      "I'll take ROAST BEEF if you're out of LAMB!!", &
      "I'm a fuschia bowling ball somewhere in Brittany", &
      "I'm a GENIUS!  I want to dispute sentence structure with SUSAN SONTAG!!", &
      "I'm a nuclear submarine under the polar ice cap and I need a Kleenex!", &
      "I'm also against BODY-SURFING!!", &
      "I'm also pre-POURED pre-MEDITATED and pre-RAPHAELITE!!", &
      "I'm an East Side TYPE..", &
      "I'm ANN LANDERS!!  I can SHOPLIFT!!", &
      "I'm ANN LANDERS!!  I can SHOPLIFT!!", &
      "I'm changing the CHANNEL..  But all I get is commercials for ``RONCO MIRACLE BAMBOO STEAMERS''!", &
      "I'm continually AMAZED at th'breathtaking effects of WIND EROSION!!", &
      "I'm CONTROLLED by the CIA!!  EVERYONE is controlled by the CIA!!", &
      "I'm definitely not in Omaha!", &
      "I'm DESPONDENT...  I hope there's something DEEP-FRIED under this miniature DOMED STADIUM...", &
      "I'm dressing up in an ill-fitting IVY-LEAGUE SUIT!!  Too late...", &
      "I'm EMOTIONAL now because I have MERCHANDISING CLOUT!!", &
      "I'm encased in the lining of a pure pork sausage!!", &
      "I'm EXCITED!!  I want a FLANK STEAK WEEK-END!!  I think I'm JULIA CHILD!!", &
      "I'm GLAD I remembered to XEROX all my UNDERSHIRTS!!", &
      "I'm gliding over a NUCLEAR WASTE DUMP near ATLANTA, Georgia!!", &
      "I'm having a BIG BANG THEORY!!", &
      "I'm having a MID-WEEK CRISIS!", &
      "I'm having a RELIGIOUS EXPERIENCE..  and I don't take any DRUGS", &
      "I'm having a tax-deductible experience!  I need an energy crunch!!", &
      "I'm having an emotional outburst!!", &
      "I'm having an EMOTIONAL OUTBURST!!  But, uh, WHY is there a WAFFLE in my PAJAMA POCKET??", &
      "I'm having BEAUTIFUL THOUGHTS about the INSIPID WIVES of smug and wealthy CORPORATE LAWYERS..", &
      "I'm having fun HITCHHIKING to CINCINNATI or FAR ROCKAWAY!!", &
      "..  I'm IMAGINING a sensuous GIRAFFE, CAVORTING in the BACK ROOM of a KOSHER DELI --", &
      "I'm in a twist contest!!  I'm in a bathtub!  It's on Mars!!  I'm in tip-top condition!", &
      "I'm in ATLANTIC CITY riding in a comfortable ROLLING CHAIR...", &
      "I'm in direct contact with many advanced fun CONCEPTS.", &
      "I'm in DISGUISE as a BAGGAGE CHECKER....I can watch the house, if it's ORANGE...", &
      "I'm in LOVE with DON KNOTTS!!", &
      "I'm into SOFTWARE!", &
      "I'm losing my hair..did it go to ATLANTIC CITY??", &
      "I'm meditating on the FORMALDEHYDE and the ASBESTOS leaking into my PERSONAL SPACE!!", &
      "I'm MENTALLY here..  but PHYSICALLY I'm purchasing NAUGAHYDE furniture in the' SUBURBS of PHOENIX!!", &
      "I'm mentally OVERDRAWN!  What's that SIGNPOST up ahead? Where's ROD STERLING when you really need him?", &
      "I'm not an Iranian!!  I voted for Dianne Feinstein!!", &
      "I'm not available for comment..", &
      "I'm pretending I'm pulling in a TROUT!  Am I doing it correctly??", &
      "I'm pretending that we're all watching PHIL SILVERS instead of RICARDO MONTALBAN!", &
      "I'm protected by a ROLL-ON I rented from AVIS..", &
      "I'm QUIETLY reading the latest issue of ``BOWLING WORLD'' while my wife and two children stand QUIETLY BY..", &
      "I'm rated PG-34!!", &
      "I'm receiving a coded message from EUBIE BLAKE!!", &
      "I'm RELIGIOUS!!  I love a man with a HAIRPIECE!! Equip me with MISSILES!!", &
      "I'm reporting for duty as a modern person.  I want to do the Latin Hustle now!", &
      "I'm shaving!!  I'M SHAVING!!", &
      "I'm sitting on my SPEED QUEEN..  To me, it's ENJOYABLE.. I'm WARM..  I'm VIBRATORY..", &
      "I'm thinking about DIGITAL READ-OUT systems and computer-generated IMAGE FORMATIONS..", &
      "I'm totally DESPONDENT over the LIBYAN situation and the price of CHICKEN..", &
      "I'm using my X-RAY VISION to obtain a rare glimpse of the INNER WORKINGS of this POTATO!!", &
      "I'm wearing PAMPERS!!", &
      "I'm wet!  I'm wild!", &
      "I'm working under the direct orders of WAYNE NEWTON to deport consenting adults!", &
      "I'm young..  I'm HEALTHY..  I can HIKE THRU CAPT GROGAN'S LUMBAR REGIONS!", &
      "I'm ZIPPY the PINHEAD and I'm totally committed to the festive mode.", &
      "I'm ZIPPY!!  Are we having FUN yet??", &
      "I've been WRITING to SOPHIA LOREN every 45 MINUTES since JANUARY 1ST!!", &
      "I've got a COUSIN who works in the GARMENT DISTRICT...", &
      "I've got an IDEA!!  Why don't I STARE at you so HARD, you forget your SOCIAL SECURITY NUMBER!!", &
      "I've got to get these SNACK CAKES to NEWARK by DAWN!!", &
      "I've gotta GO, now!!  I wanta tell you you're a GREAT bunch of guys but you ought to CHANGE your UNDERWEAR more often!!", &
      "I've read SEVEN MILLION books!!", &
      "..  ich bin in einem dusenjet ins jahr 53 vor chr... ich lande im antiken Rom...  einige gladiatoren spielen scrabble... ich rieche PIZZA...", &
      "If a person is FAMOUS in this country, they have to go on the ROAD for MONTHS at a time and have their name misspelled on the SIDE of a GREYHOUND SCENICRUISER!!", &
      "If elected, Zippy pledges to each and every American a 55-year-old houseboy...", &
      "If I am elected no one will ever have to do their laundry again!", &
      "If I am elected, the concrete barriers around the WHITE HOUSE will be replaced by tasteful foam replicas of ANN MARGARET!", &
      ".. If I cover this entire WALL with MAZOLA, wdo I have to give my AGENT ten per cent??", &
      "If I felt any more SOPHISTICATED I would DIE of EMBARRASSMENT!", &
      "If I had a Q-TIP, I could prevent th'collapse of NEGOTIATIONS!!", &
      "..  If I had heart failure right now, I couldn't be a more fortunate man!!", &
      "If I have enough money to buy 5,000 CANS of NOODLE-RONI, can I get a VAT of MARSHMALLOW FLUFF free??", &
      "If I pull this SWITCH I'll be RITA HAYWORTH!!  Or a SCIENTOLOGIST!", &
      "- if it GLISTENS, gobble it!!", &
      "If our behavior is strict, we do not need fun!", &
      "If Robert Di Niro assassinates Walter Slezak, will Jodie Foster marry Bonzo??", &
      "If this is the DATING GAME I want to know your FAVORITE PLANET!  Do I get th' MICROWAVE MOPED?", &
      "If this was a SWEDISH MOVIE, I'd take off your GO-GO BOOTS!!", &
      "If you STAY in China, I'll give you 4,000 BUSHELS of 'ATOMIC MOUSE' pencil sharpeners!!", &
      "Imagine--a WORLD without POODLES...", &
      "Impudent..  Yet possessing a certain ALUMINUM SILICATE overbite....Needs REDDY-WHIP!!", &
      "-- In 1962, you could buy a pair of SHARKSKIN SLACKS, with a ``Continental Belt,'' for $10.99!!", &
      "In Newark the laundromats are open 24 hours a day!", &
      "In order to make PLANS for the WEEKEND...so that we can read RESTAURANT REVIEWS and decide to GO to that restaurant & then NEVER GO...so we can meet a FRIEND after work in a BAR and COMPLAIN about Interior Sect'y JAMES WATT until the SUBJECT is changed to NUCLEAR BLACKMAIL...and so our RELATIVES can FORCE us to listen to HOCKEY STATISTICS while we wait for them to LEAVE on the 7:48....", &
      "INSIDE, I have the same personality disorder as LUCY RICARDO!!", &
      "Inside, I'm already SOBBING!", &
      "Intra-mural sports results are filtering through th' plumbing...", &
      "Is a tattoo real, like a curb or a battleship? Or are we suffering in Safeway?", &
      "Is it 1974?  What's for SUPPER?  Can I spend my COLLEGE FUND in one wild afternoon??", &
      "Is it clean in other dimensions?", &
      "Is it FUN to be a MIDGET?", &
      "Is it NOUVELLE CUISINE when 3 olives are struggling with a scallop in a plate of SAUCE MORNAY?", &
      "Is something VIOLENT going to happen to a GARBAGE CAN?", &
      "Is the EIGHTIES when they had ART DECO and GERALD McBOING-BOING lunch boxes??", &
      "Is there something I should be DOING with a GLAZED DONUT??", &
      "Is this 'BIKINI BEACH'?", &
      "Is this 'BOOZE'?", &
      "Is this an out-take from the ``BRADY BUNCH''?", &
      "Is this ANYWHERE, USA?", &
      "Is this BOISE??", &
      "Is this going to involve RAW human ecstasy?", &
      "Is this my STOP??", &
      "Is this TERMINAL fun?", &
      "Is this the line for the latest whimsical  YUGOSLAVIAN drama which also makes you want to CRY and reconsider the VIETNAM WAR?", &
      "Is this where people are HOT and NICE and they give you TOAST for FREE??", &
      "Isn't this my STOP?!", &
      "It don't mean a THING if you ain't got that SWING!!", &
      "It was a JOKE!!  Get it??  I was receiving messages from DAVID LETTERMAN!!  YOW!!", &
      "It's 74 degrees, 12 minutes NORTH, and 41 degrees, 3 minutes EAST!! Soon, it will be TUESDAY!!", &
      "It's a lot of fun being alive...  I wonder if my bed is made?!?", &
      "It's hard being an ARTIST!!", &
      "It's NO USE..  I've gone to ``CLUB MED''!!", &
      "It's OBVIOUS..  The FURS never reached ISTANBUL..  You were an EXTRA in the REMAKE of ``TOPKAPI''..  Go home to your WIFE..  She's making FRENCH TOAST!", &
      "It's OKAY --- I'm an INTELLECTUAL, too.", &
      "...It's REAL ROUND..  And it's got a POINTY PART right in the MIDDLE!! The shape is SMOOTH..  ..And COLD.. It feels very COMFORTABLE on my CHEEK..  I'm getting EMOTIONAL..", &
      "It's so OBVIOUS!!", &
      "It's strange, but I'm only TRULY ALIVE when I'm covered in POLKA DOTS and TACO SAUCE...", &
      "It's the land of DONNY AND MARIE as promised in TV GUIDE!", &
      "It's the RINSE CYCLE!!  They've ALL IGNORED the RINSE CYCLE!!", &
      "It's today's SPECIAL!", &
      "JAPAN is a WONDERFUL planet -- I wonder if we'll ever reach their level of COMPARATIVE SHOPPING...", &
      "Jesus is my POSTMASTER GENERAL..", &
      "Join the PLUMBER'S UNION!!", &
      "...Just enough time to do my LIBERACE impression...", &
      "Just imagine you're entering a state-of-the-art CAR WASH!!", &
      "Just to have MORE FUN, I'll pretend I am JAMES CAGNEY and I am having a tense, UP-TIGHT EXPERIENCE!!", &
      "KARL MALDEN'S NOSE just won an ACADEMY AWARD!!", &
      "Kids, don't gross me off..  ``Adventures with MENTAL HYGIENE'' can be carried too FAR!", &
      "Kids, the seven basic food groups are GUM, PUFF PASTRY, PIZZA, PESTICIDES, ANTIBIOTICS, NUTRA-SWEET and MILK DUDS!!", &
      "Laundry is the fifth dimension!!  ...um...um...  th' washing machine is a black hole and the pink socks are bus drivers who just fell in!!", &
      "LBJ, LBJ, how many JOKES did you tell today??!", &
      "Leona, I want to CONFESS things to you.. I want to WRAP you in a SCARLET ROBE trimmed with POLYVINYL CHLORIDE.. I want to EMPTY your ASHTRAYS...", &
      "Let me do my TRIBUTE to FISHNET STOCKINGS...", &
      "Let's all show human CONCERN for REVEREND MOON's legal difficulties!!", &
      "Let's climb to the TOP of that MOUNTAIN and think about STRIP MINING!!", &
      "Let's go to CHURCH!", &
      "Let's send the Russians defective lifestyle accessories!", &
      "LIFE is a never-ending INFORMERCIAL!", &
      "Life is a POPULARITY CONTEST!  I'm REFRESHINGLY CANDID!!", &
      "Life is selling REVOLUTIONARY HAIR PRODUCTS!", &
      "..  Like I always say -- nothing can beat the BRATWURST here in DUSSELDORF!!", &
      "Loni Anderson's hair should be LEGALIZED!!", &
      "Look DEEP into the OPENINGS!!  Do you see any ELVES or EDSELS... or a HIGHBALL??...", &
      "Look into my eyes and try to forget that you have a Macy's charge card!", &
      "Look!  A ladder!  Maybe it leads to heaven, or a sandwich!", &
      "Look!!  Karl Malden!", &
      "LOOK!!  Sullen American teens wearing MADRAS shorts and ``Flock of Seagulls'' HAIRCUTS!", &
      "LOOK!!!  I'm WALKING in my SLEEP again!!", &
      "LOU GRANT froze my ASSETS!!", &
      "Make me look like LINDA RONSTADT again!!", &
      "Mary Tyler Moore's SEVENTH HUSBAND is wearing my DACRON TANK TOP in a cheap hotel in HONOLULU!", &
      "Maybe we could paint GOLDIE HAWN a rich PRUSSIAN BLUE--", &
      "MERYL STREEP is my obstetrician!", &
      "MMM-MM!!  So THIS is BIO-NEBULATION!", &
      "Mmmmmm-MMMMMM!!  A plate of STEAMING PIECES of a PIG mixed with the shreds of SEVERAL CHICKENS!!...  Oh BOY!!  I'm about to swallow a TORN-OFF section of a COW'S LEFT LEG soaked in COTTONSEED OIL and SUGAR!!  ..  Let's see.. Next, I'll have the GROUND-UP flesh of CUTE, BABY LAMBS fried in the MELTED, FATTY TISSUES from a warm-blooded animal someone once PETTED!!  ...  YUM!!  That was GOOD!! For DESSERT, I'll have a TOFU BURGER with BEAN SPROUTS on a stone-ground, WHOLE WHEAT BUN!!", &
      "Mr and Mrs PED, can I borrow 26.7% of the RAYON TEXTILE production of the INDONESIAN archipelago?", &
      "My Aunt MAUREEN was a military advisor to IKE & TINA TURNER!!", &
      "My BIOLOGICAL ALARM CLOCK just went off..  It has noiseless DOZE FUNCTION and full kitchen!!", &
      "My CODE of ETHICS is vacationing at famed SCHROON LAKE in upstate New York!!", &
      "My DIGITAL WATCH has an automatic SNOOZE FEATURE!!", &
      "My EARS are GONE!!", &
      "My ELBOW is a remote FRENCH OUTPOST!!", &
      "My face is new, my license is expired, and I'm under a doctor's care!!!!", &
      "My FAVORITE group is 'QUESTION MARK & THE MYSTERIANS'...", &
      "My forehead feels like a PACKAGE of moist CRANBERRIES in a remote FRENCH OUTPOST!!", &
      "My haircut is totally traditional!", &
      "MY income is ALL disposable!", &
      "My LESLIE GORE record is BROKEN..", &
      "My LIBRARY CARD expired...", &
      "My life is a patio of fun!", &
      "My mind is a potato field...", &
      "My mind is making ashtrays in Dayton....", &
      "My nose feels like a bad Ronald Reagan movie...", &
      "..  my NOSE is NUMB!", &
      "..  My pants just went on a wild rampage through a Long Island Bowling Alley!!", &
      "My pants just went to high school in the Carlsbad Caverns!!!", &
      "My polyvinyl cowboy wallet was made in Hong Kong by Montgomery Clift!", &
      "My TOYOTA is built like a ... BAGEL with CREAM CHEESE!!", &
      "My uncle Murray conquered Egypt in 53 B.C.  And I can prove it too!!", &
      "..  My vaseline is RUNNING...", &
      "NANCY!!  Why is everything RED?!", &
      "NATHAN...  your PARENTS were in a CARCRASH!! They're VOIDED - They COLLAPSED They had no CHAINSAWS...  They had no MONEY MACHINES... They did PILLS in SKIMPY GRASS SKIRTS... Nathan, I EMULATED them...  but they were OFF-KEY...", &
      "NEWARK has been REZONED!!  DES MOINES has been REZONED!!", &
      "Nice decor!", &
      "Not enough people play SKEE-BALL..  They're always thinking about COCAINE or and ALIEN BEINGS!!", &
      "NOT fucking!! Also not a PACKAGE of LOOSE-LEAF PAPER!!", &
      "Not SENSUOUS...  only ``FROLICSOME''... and in need of DENTAL WORK...  in PAIN!!!", &
      "NOW do I get to blow out the CANLDES??", &
      "Now I am depressed...", &
      "Now I can join WEIGHT WATCHERS!", &
      "Now I need a suntan, a tennis lesson, Annette Funicello and two dozen Day-Glo orange paper jumpsuits!!", &
      "..  Now I think I just reached the state of HYPERTENSION that comes JUST BEFORE you see the TOTAL at the SAFEWAY CHECKOUT COUNTER!", &
      "Now I understand the meaning of ``THE MOD SQUAD''!", &
      "Now I'm being INVOLUNTARILY shuffled closer to the CLAM DIP with the BROKEN PLASTIC FORKS in it!!", &
      "Now I'm concentrating on a specific tank battle toward the end of World War II!", &
      "Now I'm having INSIPID THOUGHTS about the beautiful, round wives of HOLLYWOOD MOVIE MOGULS encased in PLEXIGLASS CARS and being approached by SMALL BOYS selling FRUIT..", &
      "Now I'm telling MISS PIGGY about MONEY MARKET FUNDS!", &
      "..  Now KEN and BARBIE are PERMANENTLY ADDICTED to MIND-ALTERING DRUGS..", &
      "Now KEN is having a MENTAL CRISIS beacuse his 'R.V.' PAYMENTS are OVER-DUE!!", &
      "Now my EMOTIONAL RESOURCES are heavily committed to 23% of the SMELTING and REFINING industry of the state of NEVADA!!", &
      "Now that I have my ``APPLE,'' I comprehend COST ACCOUNTING!!", &
      "Now that we're in LOVE, you can BUY this GOLDFISH for a 48% DISCOUNT.", &
      "Now, I think it would be GOOD to buy FIVE or SIX STUDEBAKERS and CRUISE for ARTIFICIAL FLAVORING!!", &
      "NOW, I'm supposed to SCRAMBLE two, and HOLD th' MAYO!!", &
      "NOW, I'm taking the NEXT FLIGHT to ACAPULCO so I can write POEMS about BROKEN GUITAR STRINGS and sensuous PRE-TEENS!!", &
      "Now, let's SEND OUT for QUICHE!!", &
      "Now, my ENTIRE LIFE is flashing before my EYES as I park my DODGE DART in your EXXON service area for a COMPLETE LUBRICATION!!", &
      "O.K.!  Speak with a PHILADELPHIA ACCENT!!  Send out for CHINESE FOOD!! Hop a JET!", &
      "Of course, you UNDERSTAND about the PLAIDS in the SPIN CYCLE --", &
      "Oh my GOD -- the SUN just fell into YANKEE STADIUM!!", &
      "Oh, FISH sticks, CHEEZ WHIZ, GIN fizz, SHOW BIZ!!", &
      "Oh, I get it!!  ``The BEACH goes on,'' huh, SONNY??", &
      "OKAY!!  Turn on the sound ONLY for TRYNEL CARPETING, FULLY-EQUIPPED R.V.'S and FLOATATION SYSTEMS!!", &
      "Okay, BARBRA STREISAND, I recognize you now!!  Also EFREM ZIMBALIST, JUNIOR!!  And BEAUMONT NEWHALL!!  Everybody into th' BATHROOM!", &
      "Okay..  I'm going home to write the ``I HATE RUBIK's CUBE HANDBOOK FOR DEAD CAT LOVERS''..", &
      "OMNIVERSAL AWARENESS??  Oh, YEH!!  First you need 4 GALLONS of JELL-O and a BIG WRENCH!!...  I think you drop th'WRENCH in the JELL-O as if it was a FLAVOR, or an INGREDIENT...  ...or...I...um...  WHERE'S the WASHING MACHINES?", &
      "On SECOND thought, maybe I'll heat up some BAKED BEANS and watch REGIS PHILBIN..  It's GREAT to be ALIVE!!", &
      "On the other hand, life can be an endless parade of TRANSSEXUAL QUILTING BEES aboard a cruise ship to DISNEYWORLD if only we let it!!", &
      "On the road, ZIPPY is a pinhead without a purpose, but never without a POINT.", &
      "..  Once upon a time, four AMPHIBIOUS HOG CALLERS attacked a family of DEFENSELESS, SENSITIVE COIN COLLECTORS and brought DOWN their PROPERTY VALUES!!", &
      "Once, there was NO fun...  This was before MENU planning, FASHION statements or NAUTILUS equipment... Then, in 1985..  FUN was completely encoded in this tiny MICROCHIP.. It contain 14,768 vaguely amusing SIT-COM pilots!! We had to wait FOUR BILLION years but we finally got JERRY LEWIS, MTV and a large selection of creme-filled snack cakes!", &
      "..  One FISHWICH coming up!!", &
      "ONE:  I will donate my entire ``BABY HUEY'' comic book collection to the downtown PLASMA CENTER.. TWO:  I won't START a BAND called ``KHADAFY & THE HIT SQUAD''.. THREE:  I won't ever TUMBLE DRY my FOX TERRIER again!!", &
      "..  or were you driving the PONTIAC that HONKED at me in MIAMI last Tuesday?", &
      "Our father who art in heaven..  I sincerely pray that SOMEBODY at this table will PAY for my SHREDDED WHAT and ENGLISH MUFFIN.. and also leave a GENEROUS TIP...", &
      "..  over in west Philadelphia a puppy is vomiting..", &
      "OVER the underpass!  UNDER the overpass!  Around the FUTURE and BEYOND REPAIR!!", &
      "PARDON me, am I speaking ENGLISH?", &
      "Pardon me, but do you know what it means to be TRULY ONE with your BOOTH!", &
      "PEGGY FLEMING is stealing BASKET BALLS to feed the babies in VERMONT.", &
      "...PENGUINS are floating by...", &
      "PIZZA!!", &
      "Place me on a BUFFER counter while you BELITTLE several BELLHOPS in the Trianon Room!!  Let me one of your SUBSIDIARIES!", &
      "Please come home with me...  I have Tylenol!!", &
      "Psychoanalysis??  I thought this was a nude rap session!!!", &
      "PUNK ROCK!!  DISCO DUCK!!  BIRTH CONTROL!!", &
      "Put FIVE DOZEN red GIRDLES in each CIRCULAR OPENING!!", &
      "Quick, sing me the BUDAPEST NATIONAL ANTHEM!!", &
      "QUIET!!  I'm being CREATIVE!!  Is it GREAT yet?  It's s'posed to SMOKEY THE BEAR...", &
      "RELATIVES!!", &
      "RELAX!! ... This is gonna be a HEALING EXPERIENCE!!  Besides, I work for DING DONGS!", &
      "Remember, if you try to ESCAPE, many APARTMENT HOPPING ALCOHOLICS will SIMONIZE your HALLWAYS!!  This is your LAST WARNING!!", &
      "Remember, in 2039, MOUSSE & PASTA will be available ONLY by prescription!!", &
      "RHAPSODY in Glue!", &
      "SANTA CLAUS comes down a FIRE ESCAPE wearing bright blue LEG WARMERS..  He scrubs the POPE with a mild soap or detergent for 15 minutes, starring JANE FONDA!!", &
      "Send your questions to ``ASK ZIPPY'', Box 40474, San Francisco, CA 94140, USA", &
      "SHHHH!!  I hear SIX TATTOOED TRUCK-DRIVERS tossing ENGINE BLOCKS into empty OIL DRUMS..", &
      "Should I do my BOBBIE VINTON medley?", &
      "..  Should I get locked in the PRINCIPAL'S OFFICE today -- or have a VASECTOMY??", &
      "Should I start with the time I SWITCHED personalities with a BEATNIK hair stylist or my failure to refer five TEENAGERS to a good OCULIST?", &
      "Sign my PETITION.", &
      "So this is what it feels like to be potato salad", &
      "..  So, if we convert SUPPLY-SIDE SOYBEAN FUTURES into HIGH-YIELD T-BILL INDICATORS, the PRE-INFLATIONARY risks will DWINDLE to a rate of 2 SHOPPING SPREES per EGGPLANT!!", &
      "..  someone in DAYTON, Ohio is selling USED CARPETS to a SERBO-CROATIAN", &
      "Someone is DROOLING on my collar!!", &
      "Sometime in 1993 NANCY SINATRA will lead a BLOODLESS COUP on GUAM!!", &
      "Somewhere in DOWNTOWN BURBANK a prostitute is OVERCOOKING a LAMB CHOP!!", &
      "Somewhere in suburban Honolulu, an unemployed bellhop is whipping up a batch of illegal psilocybin chop suey!!", &
      "Somewhere in Tenafly, New Jersey, a chiropractor is viewing ``Leave it to Beaver''!", &
      "Sorry, wrong ZIP CODE!!", &
      "Spreading peanut butter reminds me of opera!!  I wonder why?", &
      "TAILFINS!!  ...click...", &
      "Talking Pinhead Blues: Oh, I LOST my ``HELLO KITTY'' DOLL and I get BAD reception on  channel TWENTY-SIX!! Th'HOSTESS FACTORY is closin' down and I just heard ZASU PITTS  has been DEAD for YEARS..  (sniff) My PLATFORM SHOE collection was CHEWED up by th'dog, ALEXANDER  HAIG won't let me take a SHOWER 'til Easter.. (snurf) So I went to the kitchen, but WALNUT PANELING whup me  upside mah HAID!! (on no, no, no..  Heh, heh)", &
      "TAPPING?  You POLITICIANS!  Don't you realize that the END of the ``Wash Cycle'' is a TREASURED MOMENT for most people?!", &
      "TATTOOED MIDGETS are using ALFREDO in their SALAMI FACTORY!", &
      "Tex SEX!  The HOME of WHEELS!  The dripping of COFFEE!!  Take me to Minnesota but don't EMBARRASS me!!", &
      "Th' MIND is the Pizza Palace of th' SOUL", &
      "Th' PINK SOCK... soaking... soaking... soaking... Th' PINK SOCK... washing... washing... washing... Th' PINK SOCK... rinsing... rinsing... rinsing...", &
      "Thank god!!..  It's HENNY YOUNGMAN!!", &
      "That's a decision that can only be made between you & SY SPERLING!!", &
      "The appreciation of the average visual graphisticator alone is worth the whole suaveness and decadence which abounds!!", &
      "The entire CHINESE WOMEN'S VOLLEYBALL TEAM all share ONE personality -- and have since BIRTH!!", &
      "The fact that 47 PEOPLE are yelling and sweat is cascading down my SPINAL COLUMN is fairly enjoyable!!", &
      "The FALAFEL SANDWICH lands on my HEAD and I become a VEGETARIAN...", &
      "..  the HIGHWAY is made out of LIME JELLO and my HONDA is a barbequed OYSTER!  Yum!", &
      "The Korean War must have been fun.", &
      "'THE LITTLE PINK FLESH SISTERS,' I saw them at th' FLUROESCENT BULB MAKERS CONVENTION...", &
      "The LOGARITHM of an ISOCELES TRIANGLE is TUESDAY WELD!!", &
      "..  the MYSTERIANS are in here with my CORDUROY SOAP DISH!!", &
      "The Osmonds!  You are all Osmonds!!  Throwing up on a freeway at dawn!!!", &
      "The PILLSBURY DOUGHBOY is CRYING for an END to BURT REYNOLDS movies!!", &
      "The PINK SOCKS were ORIGINALLY from 1952!! But they went to MARS around 1953!!", &
      "The SAME WAVE keeps coming in and COLLAPSING like a rayon MUU-MUU..", &
      "..The TENSION mounts as I MASSAGE your RIGHT ANKLE according to ancient Tibetan ACCOUNTING PROCEDURES..are you NEUROTIC yet??", &
      "Then, it's off to RED CHINA!!", &
      "There's a little picture of ED MCMAHON doing BAD THINGS to JOAN RIVERS in a $200,000 MALIBU BEACH HOUSE!!", &
      "There's a lot of BIG MONEY in MISERY if you have an AGENT!!", &
      "There's a SALE on STRETCH SOCKS down at the '7-11'!!", &
      "There's enough money here to buy 5000 cans of Noodle-Roni!", &
      "These PRESERVES should be FORCE-FED to PENTAGON OFFICIALS!!", &
      "They collapsed....  like nuns in the street... they had no teen appeal!", &
      "They don't hire PERSONAL PINHEADS, Mr. Toad!", &
      "This ASEXUAL PIG really BOILS my BLOOD...  He's so..so.....URGENT!!", &
      "This is a NO-FRILLS flight -- hold th' CANADIAN BACON!!", &
      "This is my WILLIAM BENDIX memorial CORNER where I worship William Bendix like a GOD!!", &
      "This is PLEASANT!", &
      "This MUST be a good party -- My RIB CAGE is being painfully pressed up against someone's MARTINI!!", &
      "..  this must be what it's like to be a COLLEGE GRADUATE!!", &
      "This PIZZA symbolizes my COMPLETE EMOTIONAL RECOVERY!!", &
      "This PORCUPINE knows his ZIPCODE..  And he has ``VISA''!!", &
      "This TOPS OFF my partygoing experience!  Someone I DON'T LIKE is talking to me about a HEART-WARMING European film..", &
      "Those aren't WINOS--that's my JUGGLER, my AERIALIST, my SWORD SWALLOWER, and my LATEX NOVELTY SUPPLIER!!", &
      "Thousands of days of civilians ...  have produced a... feeling for the aesthetic modules --", &
      "Three attractive BANK ROBBERS are discussing RELIGIOUS DIFFERENCES and MAKE-UP TECHNIQUE with them!!", &
      "Today, THREE WINOS from DETROIT sold me a framed photo of TAB HUNTER before his MAKEOVER!", &
      "Toes, knees, NIPPLES.  Toes, knees, nipples, KNUCKLES... Nipples, dimples, knuckles, NICKLES, wrinkles, pimples!! I don't like FRANK SINATRA or his CHILDREN.", &
      "TONY RANDALL!  Is YOUR life a PATIO of FUN??", &
      "Two LITTLE black dots and one BIG black dot...nice 'n' FLUFFY!!", &
      "Two with FLUFFO, hold th' BEETS..side of SOYETTES!", &
      "Uh-oh --  WHY am I suddenly thinking of a VENERABLE religious leader frolicking on a FORT LAUDERDALE weekend?", &
      "Uh-oh!!  I forgot to submit to COMPULSORY URINALYSIS!", &
      "UH-OH!!  I put on ``GREAT HEAD-ON TRAIN COLLISIONS of the 50's'' by mistake!!!", &
      "UH-OH!!  I think KEN is OVER-DUE on his R.V. PAYMENTS and HE'S having a NERVOUS BREAKDOWN too!!  Ha ha.", &
      "Uh-oh!!  I'm having TOO MUCH FUN!!", &
      "UH-OH!!  We're out of AUTOMOBILE PARTS and RUBBER GOODS!", &
      "...Um...Um...", &
      "Used staples are good with SOY SAUCE!", &
      "Vote for ME -- I'm well-tapered, half-cocked, ill-conceived and TAX-DEFERRED!", &
      "..Wait 'til those  ITALIAN TEENAGERS get back to their HONDAS & discover them to be FILLED to the BRIM with MAZOLA!!", &
      "Wait..  is this a FUN THING or the END of LIFE in Petticoat Junction??", &
      "Was my SOY LOAF left out in th'RAIN?  It tastes REAL GOOD!!", &
      "We are now enjoying total mutual interaction in an imaginary hot tub...", &
      "We have DIFFERENT amounts of HAIR --", &
      "We just joined the civil hair patrol!", &
      "We place two copies of PEOPLE magazine in a DARK, HUMID mobile home. 45 minutes later CYNDI LAUPER emerges wearing a BIRD CAGE on her head!", &
      "Well, here I am in AMERICA..  I LIKE it.  I HATE it.  I LIKE it.  I HATE it.  I LIKE it.  I HATE it.  I LIKE it.I HATE it.  I LIKE..  EMOTIONS are SWEEPING over me!!", &
      "Well, I'm a classic ANAL RETENTIVE!!  And I'm looking for a way to VICARIOUSLY experience some reason to LIVE!!", &
      "Well, I'm INVISIBLE AGAIN..  I might as well pay a visit to the LADIES ROOM...", &
      "Well, I'm on the right planet---everyone looks like me!!!", &
      "Well, O.K.  I'll compromise with my principles because of EXISTENTIAL DESPAIR!", &
      "Were these parsnips CORRECTLY MARINATED in TACO SAUCE?", &
      "What a COINCIDENCE!  I'm an authorized ``SNOOTS OF THE STARS'' dealer!!", &
      "What GOOD is a CARDBOARD suitcase ANYWAY?", &
      "What I need is a MATURE RELATIONSHIP with a FLOPPY DISK...", &
      "What I want to find out is -- do parrots know much about Astro-Turf?", &
      "What PROGRAM are they watching?", &
      "What UNIVERSE is this, please??", &
      "What's the MATTER Sid?..  Is your BEVERAGE unsatisfactory?", &
      "When I met th'POPE back in '58, I scrubbed him with a MILD SOAP or DETERGENT for 15 minutes.  He seemed to enjoy it..", &
      "When this load is DONE I think I'll wash it AGAIN..", &
      "When you get your PH.D. will you get able to work at BURGER KING?", &
      "When you said ``HEAVILY FORESTED'' it reminded me of an overdue CLEANING BILL..  Don't you SEE?  O'Grogan SWALLOWED a VALUABLE COIN COLLECTION and HAD to murder the ONLY MAN who KNEW!!", &
      "Where do your SOCKS go when you lose them in th' WASHER?", &
      "Where does it go when you flush?", &
      "Where's my SOCIAL WORKER?", &
      "Where's SANDY DUNCAN?", &
      "Where's th' DAFFY DUCK EXHIBIT??", &
      "Where's the Coke machine?  Tell me a joke!!", &
      "While I'm in LEVITTOWN I thought I'd like to see the NUCLEAR FAMILY!!", &
      "While my BRAINPAN is being refused service in BURGER KING, Jesuit priests are DATING CAREER DIPLOMATS!!", &
      "While you're chewing, think of STEVEN SPIELBERG'S bank account..  This will have the same effect as two ``STARCH BLOCKERS''!", &
      "WHO sees a BEACH BUNNY sobbing on a SHAG RUG?!", &
      "Who wants some OYSTERS with SEN-SEN an' COOL WHIP?", &
      "WHOA!!  I'm having a RELIGIOUS EXPERIENCE right NOW!!", &
      "WHOA!!  Ken and Barbie are having TOO MUCH FUN!!  It must be the NEGATIVE IONS!!", &
      "Why am I in this ROOM in DOWNTOWN PHILADELPHIA?", &
      "Why are these athletic shoe salesmen following me??", &
      "WHY are we missing KOJAK?", &
      "Why don't you ever enter and CONTESTS, Marvin?? Don't you know your own ZIPCODE?", &
      "Why is everything made of Lycra Spandex?", &
      "Why is it that when you DIE, you can't take your HOME ENTERTAINMENT CENTER with you??", &
      "Why was I BORN?", &
      "Will it improve my CASH FLOW?", &
      "Will the third world war keep ``Bosom Buddies'' off the air?", &
      "Will this never-ending series of PLEASURABLE EVENTS never cease?", &
      "With this weapon I can expose fictional characters and bring about sweeping reforms!!", &
      "With YOU, I can be MYSELF..  We don't NEED Dan Rather..", &
      "World War Three can be averted by adherence to a strictly enforced dress code!", &
      "Wow!  Look!!  A stray meatball!!  Let's interview it!", &
      "Xerox your lunch and file it under ``sex offenders!''", &
      "Yes, but will I see the EASTER BUNNY in skintight leather at an IRON MAIDEN concert?", &
      "Yes, Private DOBERMAN!!", &
      "You can't hurt me!!  I have an ASSUMABLE MORTGAGE!!", &
      "You mean now I can SHOOT YOU in the back and further BLUR th' distinction between FANTASY and REALITY?", &
      "You mean you don't want to watch WRESTLING from ATLANTA?", &
      "You must be a CUB SCOUT!!  Have you made your MONEY-DROP today??", &
      "YOU PICKED KARL MALDEN'S NOSE!!", &
      "You should all JUMP UP AND DOWN for TWO HOURS while I decide on a NEW CAREER!!", &
      "You were s'posed to laugh!", &
      "YOU!!  Give me the CUTEST, PINKEST, most charming little VICTORIAN DOLLHOUSE you can find!!  An make it SNAPPY!!", &
      "YOU'D cry too if it happened to YOU!!", &
      "Your CHEEKS sit like twin NECTARINES above a MOUTH that knows no BOUNDS --", &
      "Youth of today!  Join me in a mass rally for traditional mental attitudes!", &
      "Yow!", &
      "Yow!  Am I cleansed yet?!", &
      "Yow!  Am I having fun yet?", &
      "Yow!  Am I in Milwaukee?", &
      "Yow!  Am I JOGGING yet??", &
      "Yow!  And then we could sit on the hoods of cars at stop lights!", &
      "Yow!  Are we in the perfect mood?", &
      "Yow!  Are we laid back yet?", &
      "Yow!  Are we wet yet?", &
      "Yow!  Are you the self-frying president?", &
      "Yow!  Did something bad happen or am I in a drive-in movie??", &
      "YOW!  I can see 1987!!  PRESIDENT FORD is doing the REMAKE of 'PAGAN LOVE SONG'...he's playing ESTHER WILLIAMS!!", &
      "Yow!  I forgot my PAIL!!", &
      "Yow!  I just went below the poverty line!", &
      "Yow!  I like my new DENTIST...", &
      "Yow!  I threw up on my window!", &
      "Yow!  I want my nose in lights!", &
      "Yow!  I want to mail a bronzed artichoke to Nicaragua!", &
      "Yow!  I'm having a quadraphonic sensation of two winos alone in a steel mill!", &
      "Yow!  I'm imagining a surfer van filled with soy sauce!", &
      "Yow!  I'm out of work...I could go into shock absorbers...or SCUBA GEAR!!", &
      "Yow!  I'm UNEMPLOYED!", &
      "Yow!  Is my fallout shelter termite proof?", &
      "Yow!  Is this sexual intercourse yet??  Is it, huh, is it??", &
      "Yow!  It's a hole all the way to downtown Burbank!", &
      "Yow!  It's some people inside the wall!  This is better than mopping!", &
      "Yow!  Maybe I should have asked for my Neutron Bomb in PAISLEY--", &
      "Yow!  Now I get to think about all the BAD THINGS I did to a BOWLING BALL when I was in JUNIOR HIGH SCHOOL!", &
      "Yow!  Now we can become alcoholics!", &
      "Yow!  STYROFOAM..", &
      "Yow!  Those people look exactly like Donnie and Marie Osmond!!", &
      "Yow!  We're going to a new disco!", &
      "YOW!!", &
      "Yow!!  'Janitor trapped in sewer uses ESP to find decayed burger'!!", &
      "YOW!!  Everybody out of the GENETIC POOL!", &
      "YOW!!  I am having FUN!!", &
      "YOW!!  I'm in a very clever and adorable INSANE ASYLUM!!", &
      "Yow!!  It's LIBERACE and TUESDAY WELD!!  High on a HILL... driving a LITTLE CAR...  I wanna be in that LITTLE CAR, too!!  I wanna drive off with LIBBY and TUESDAY!", &
      "YOW!!  Now I understand advanced MICROBIOLOGY and th' new TAX REFORM laws!!", &
      "YOW!!  Now I'm playing with my HOLOGRAPHIC ATOMIC SIMULATION LASER pinball machine!!  WORLD PEACE is in the BALANCE!!", &
      "Yow!!  That's a GOOD IDEA!!  Eating a whole FIELD of COUGH MEDICINE should make you feel MUCH BETTER!!", &
      "YOW!!  The land of the rising SONY!!", &
      "YOW!!  Up ahead!  It's a DONUT HUT!!", &
      "YOW!!  What should the entire human race DO??  Consume a fifth of CHIVAS REGAL, ski NUDE down MT. EVEREST, and have a wild SEX WEEKEND!", &
      "YOW!!!  I am having fun!!!", &
      "YUGGA-HUGGA-BUGGA-TUGGA!!  HEY-HEY!!  A TRAIN STATION!!  No, a POST OFFICE!!  An OCEAN LINER!!  No, I think it's a CAFETERIA!!!", &
      "Zippy's brain cells are straining to bridge synapses..." &
    /)

    CALL MondoLogPlain("Zippy sez: "//TRIM(quotes(Random((/ 1, SIZE(quotes) /)))))

  END SUBROUTINE ZippySez
END MODULE ZippyQuote
