%%
%% This is file `revtex4.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% revtex4.dtx  (with options: `class')
%% 
%% This file is part of the APS files in the REVTeX 4 distribution.
%% For the version number, search on the string 
%% Original version by David Carlisle
%% Modified by Arthur Ogawa (mailto:ogawa@teleport.com)
%% 
%% Copyright (c) 1999 The American Physical Society.
%% http://publish.aps.org/revtex4/
%% mailto:revtex4@aps.org
%% 
%% See the REVTeX 4 README file for restrictions and more information.
%% 
\NeedsTeXFormat{LaTeX2e}[1996/06/01]
\ProvidesClass{revtex4}
          [1999/09/13 v4.0beta 2 RevTeX class (DPC,AO,APS)]%FileInfo
\def\class@name{revtex4}%
\if@compatibility
  \ifx\every@math@size\@undefined
  \else
    \let\old@expast\@expast
    \def\@expast#1{\old@expast{#1}\let\@tempa\reserved@a}
  \fi
  \input{revtex.sty}
  \expandafter\endinput
\fi
\def\REVTEX@err#1{\ClassError{\class@name}{#1}\@eha}%
\def\REVTEX@warn#1{\ClassWarningNoLine{\class@name}{#1}}%
\def\REVTEX@info#1{\ClassInfo{\class@name}{#1}}%
\def\obsolete@command#1{%
 \REVTEX@warn@end{Command \string#1\space is obsolete.^^JPlease remove from your document.}%
 \global\let#1\@empty
 #1%
}%
\def\replace@command#1#2{%
 \REVTEX@warn@end{Command \string#1\space is obsolete;^^JUse \string#2\space instead.}%
 \global\let#1#2%
 #1%
}%
\def\REVTEX@warn@end#1{\AtEndDocument{\REVTEX@warn{#1}}}%
\def\halignt@{\halign to}%
\chardef\f@ur=4\relax
\def\let@environment#1#2{%
 \expandafter\let
 \csname#1\expandafter\endcsname\csname#2\endcsname
 \expandafter\let
 \csname end#1\expandafter\endcsname\csname end#2\endcsname
 }%
\newcommand\tracingplain{%
 \tracingonline\z@\tracingcommands\z@\tracingstats\z@
 \tracingpages\z@\tracingoutput\z@\tracinglostchars\@ne
 \tracingmacros\z@\tracingparagraphs\z@\tracingrestores\z@
 \showboxbreadth5\showboxdepth3\relax %\errorstopmode
 }%
\newcommand\traceoutput{%
 \appdef\@resetactivechars{\showoutput}%
}%
\def\prepdef#1#2{%
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}%
 \toks@ii{#2}%
 \edef#1{\the\toks@ii\the\toks@}%
}%
\def\appdef#1#2{%
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}%
 \toks@ii{#2}%
 \edef#1{\the\toks@\the\toks@ii}%
}%
\def\gappdef#1#2{%
 \@ifxundefined#1{\toks@{}}{\toks@\expandafter{#1}}%
 \toks@ii{#2}%
 \global\edef#1{\the\toks@\the\toks@ii}%
}%
\def\appdef@val#1#2{%
 \appdef#1{{#2}}%
}%
\def\appdef@e#1#2{%
 \expandafter\appdef
 \expandafter#1%
 \expandafter{#2}%
}%
\def\appdef@eval#1#2{%
 \expandafter\appdef@val
 \expandafter#1%
 \expandafter{#2}%
}%
\toksdef\toks@ii=\tw@
\long\def\@ifxundefined#1{\@ifx{\undefined#1}}%
\long\def\@ifnotrelax#1#2#3{\@ifx{\relax#1}{#3}{#2}}%
\long\def\@argswap#1#2{#2#1}%
\long\def\@argswap@val#1#2{#2{#1}}%
\def\@ifxundefined@cs#1{\expandafter\@ifx\expandafter{\csname#1\endcsname\relax}}%
\def\@boolean#1#2{%
  \long\def#1{%
    #2% \if<something>
      \expandafter\boolean@true
    \else
      \expandafter\boolean@false
    \fi
  }%
}%
\def\@boole@def#1#{\@boolean{#1}}% Implicit #2
\def\@booleantrue#1{\let#1\boolean@true}%
\def\@booleanfalse#1{\let#1\boolean@false}%
\@boole@def\@ifx#1{\ifx#1}%
\@boole@def\@ifx@empty#1{\ifx\@empty#1}%
\@boole@def\@ifempty#1{\if!#1!}%
\@boole@def\@if@sw#1{\csname if#1\endcsname}%
\@boole@def\@ifnum#1{\ifnum#1}%
\@boole@def\@ifodd#1{\ifodd#1}%
\@boole@def\@ifdim#1{\ifdim#1}%
\@boole@def\@ifvoid#1{\ifvoid#1}%
\@boole@def\@ifhmode{\ifhmode}%
\@boole@def\@ifvmode{\ifvmode}%
\long\def\boolean@true#1#2{#1}%
\long\def\boolean@false#1#2{#2}%
\prepdef\document{%
 \endgroup
 \set@typesize@hook
 \normalsize
 \set@pica@hook
 \boolean@true{}%
}%
\def\set@typesize@hook{}%
\def\set@pica@hook{}%
\def\eqnarray@broken{%
   \stepcounter{equation}%
   \def\@currentlabel{\p@equation\theequation}%
   \global\@eqnswtrue
   \m@th
   \global\@eqcnt\z@
   \tabskip\@centering
   \let\\\@eqncr
   $$\everycr{}\halign to\displaywidth\bgroup
       \hskip\@centering$\displaystyle\tabskip\z@skip{##}$\@eqnsel
      &\global\@eqcnt\@ne\hskip \tw@\arraycolsep \hfil${##}$\hfil
      &\global\@eqcnt\tw@ \hskip \tw@\arraycolsep
         $\displaystyle{##}$\hfil\tabskip\@centering
      &\global\@eqcnt\thr@@ \hb@xt@\z@\bgroup\hss##\egroup
         \tabskip\z@skip
      \cr
}
\long\def\eqnarray@fleqn@fixed{%
 \stepcounter{equation}\def\@currentlabel{\p@equation\theequation}%
 \global\@eqnswtrue\m@th\global\@eqcnt\z@
 \tabskip\mathindent
 \let\\=\@eqncr
 \setlength\abovedisplayskip{\topsep}%
 \ifvmode\addtolength\abovedisplayskip{\partopsep}\fi
 \addtolength\abovedisplayskip{\parskip}%
 \setlength\belowdisplayskip{\abovedisplayskip}%
 \setlength\belowdisplayshortskip{\abovedisplayskip}%
 \setlength\abovedisplayshortskip{\abovedisplayskip}%
 $$%
 \everycr{}%
 \halignt@\linewidth\bgroup
  \hskip\@centering$\displaystyle\tabskip\z@skip{##}$\@eqnsel
  &\global\@eqcnt\@ne
   \hskip\tw@\eqncolsep
   \hfil${{}##{}}$\hfil
  &\global\@eqcnt\tw@
   \hskip\tw@\eqncolsep
   $\displaystyle{##}$\hfil\tabskip\@centering
  &\global\@eqcnt\thr@@\hb@xt@\z@\bgroup\hss##\egroup
   \tabskip\z@skip
  \cr
}%
\@ifx{\eqnarray\eqnarray@broken}{%
 \expandafter\PackageInfo\expandafter{\class@name}%
  {Repairing broken {eqnarray} definition}%
 \let\eqnarray\eqnarray@fleqn@fixed
 \newlength\eqncolsep
 \setlength\eqncolsep\z@
 \let\eqnarray@broken\relax
 \let\eqnarray@fleqn@fixed\relax
}{}%
\def\mathindent{\@centering}%
\def\set@eqnarray@skips{}%
\AtBeginDocument{%
 \prepdef\caption{\minipagefootnotes}%
}%
\def\minipagefootnotes{%
    \par
    \ifvoid\@mpfootins\else
      \unvbox\@mpfootins
    \fi}
\let\REVTEX@xfloat\@xfloat
\def\@xfloat#1[#2]{%
  \par
  \let\footnote\REVTEX@footnote
  \def\@mpfn{mpfootnote}%
  \def\thempfn{\thempfootnote}%
  \c@mpfootnote\z@
  \let\@footnotetext\@mpfootnotetext
  \let\@makefntext\@mpmakefntext
  \@ifxundefined\floats@sw{\global\@booleantrue\floats@sw}{}%
  \floats@sw{%
   \REVTEX@xfloat#1[#2]%
  }{%
   \def\@captype{#1}%
   \addvspace\intextsep
   \vbox\bgroup
    \def\end@float{%
      \minipagefootnotes
      \par\vskip\z@skip      %% \par\vskip\z@ added 15 Dec 87
      \egroup
      \par
      \addvspace\intextsep
    }%
  }%
}%
\def\@mpmakefntext#1{%
 \parindent=1em
 \noindent
 \hb@xt@1em{\hss\@makefnmark}%
 #1%
}%
\prepdef\end@float{\minipagefootnotes}%
\def\do@if@floats#1#2{%
 \@ifxundefined\floats@sw{\global\@booleantrue\floats@sw}{}%
 \floats@sw{}{%
  \expandafter\newwrite
              \csname#1write\endcsname
  \expandafter\def
              \csname#1@stream\endcsname{\jobname#2}%
  \expandafter\immediate
  \expandafter\openout
              \csname#1write\endcsname
              \csname#1@stream\endcsname\relax
  \@ifxundefined\REVTEX@float{%
   \let\REVTEX@float\@float
   \let\REVTEX@dblfloat\@dblfloat
   \let\@float\write@float
   \let\@dblfloat\write@floats
  }{}%
  \let@environment{#1@float}{#1}%
  \let@environment{#1@floats}{#1*}%
  \@ifxundefined@cs{#1@write}{}{%
   \let@environment{#1}{#1@write}%
  }%
 }%
}%
\def\print@float#1{%
 \@ifxundefined@cs{#1write}{}{%
  \begingroup
   \@ifxundefined@cs{#1sname}{}{%
    \expandafter\section
    \expandafter*%
    \expandafter{%
                 \csname#1sname\endcsname
                }%
   }%
   \raggedbottom
   \let\@float\REVTEX@float
   \let\@dblfloat\REVTEX@dblfloat
   \let@environment{#1}{#1@float}%
   \let@environment{#1*}{#1@floats}%
   \expandafter\immediate
   \expandafter\closeout
               \csname#1write\endcsname
   \input{\csname#1@stream\endcsname}%
  \endgroup
  \expandafter\let\csname#1write\endcsname\relax
 }%
}%
\def\write@float#1{\write@@float{#1}{#1}}%
\def\endwrite@float{\@Esphack}%
\def\write@floats#1{\write@@float{#1*}{#1}}%
\def\endwrite@floats{\@Esphack}%
\def\write@@float#1#2{%
  \ifhmode
     \@bsphack
  \fi
  \chardef\@tempc\csname#2write\endcsname
  \toks@{\begin{#1}}%
  \def\@tempb{#1}%
  \expandafter\let\csname end#1\endcsname\endwrite@float
  \catcode`\^^M\active
  \@makeother\{\@makeother\}\@makeother\%
  \write@floatline
}%
\begingroup
 \catcode`\[\the\catcode`\{\catcode`\]\the\catcode`\}\@makeother\{\@makeother\}%
 \gdef\float@end@tag#1\end{#2}#3\@nul[%
  \def\@tempa[#2]%
  \@ifx[\@tempa\@tempb][\end[#2]][\write@floatline]%
 ]%
 \obeylines%
 \gdef\write@floatline#1^^M[%
  \begingroup%
   \newlinechar`\^^M%
   \toks@\expandafter[\the\toks@#1]\immediate\write\@tempc[\the\toks@]%
  \endgroup%
  \toks@[]%
  \float@end@tag#1\end{}\@nul%
 ]%
\endgroup
\def\@alph#1{\ifcase#1\or a\or b\or c\or d\else\@ialph{#1}\fi}
\def\@ialph#1{\ifcase#1\or \or \or \or \or e\or f\or g\or h\or i\or j\or
  k\or l\or m\or n\or o\or p\or q\or r\or s\or t\or u\or v\or w\or x\or
  y\or z\or aa\or bb\or cc\or dd\or ee\or ff\or gg\or hh\or ii\or jj\or
  kk\or ll\or mm\or nn\or oo\or pp\or qq\or rr\or ss\or tt\or uu\or
  vv\or ww\or xx\or yy\or zz\else\@ctrerr\fi}
\def\@fnsymbol#1{\ensuremath{\ifcase#1\or *\or \dagger\or \ddagger\or
  \mathchar "278\or **\or \dagger\dagger
  \or \ddagger\ddagger \or \mathsection \mathsection\or
  *{*}*\or\dagger\dagger\dagger\or\ddagger\ddagger\ddagger\or
  \mathsection \mathsection \mathsection
  \else\@ctrerr\fi}}
\def\@startsection#1#2#3#4#5#6{%
 \@startsection@hook
 \if@noskipsec \leavevmode \fi
 \par
 \@tempskipa #4\relax
 \@afterindenttrue
 \ifdim \@tempskipa <\z@
  \@tempskipa -\@tempskipa \@afterindentfalse
 \fi
 \if@nobreak
  \everypar{}%
 \else
  \addpenalty\@secpenalty\addvspace\@tempskipa
 \fi
 \@ifstar
  {\@dblarg{\@ssect{#1}{#2}{#3}{#4}{#5}{#6}}}%
  {\@dblarg{\@sect {#1}{#2}{#3}{#4}{#5}{#6}}}%
}%
\def\@startsection@hook{}%
\REVTEX@info
  {Repairing broken definition of \string\@sect.}%
\def\@sect#1#2#3#4#5#6[#7]#8{%
  \@ifnum{#2>\c@secnumdepth}{%
    \let\@svsec\@empty
  }{%
    \refstepcounter{#1}%
    \protected@edef\@svsec{{#1}}%
    \@ifundefined{@#1cntformat}{%
     \prepdef\@svsec\@seccntformat
    }{%
     \expandafter\prepdef
     \expandafter\@svsec
                 \csname @#1cntformat\endcsname
    }%
  }%
  \@tempskipa #5\relax
  \@ifdim{\@tempskipa>\z@}{%
    \begingroup
      \interlinepenalty \@M
      #6{%
       \@ifundefined{@hangfrom@#1}{\@hangfrom}{\csname @hangfrom@#1\endcsname}%
       {\hskip#3\relax\@svsec}{#8}%
      }%
      \@@par
    \endgroup
    \@ifundefined{#1mark}{\@gobble}{\csname #1mark\endcsname}{#7}%
    \addcontentsline{toc}{#1}{%
      \ifnum #2>\c@secnumdepth \else
        \protect\numberline{\csname the#1\endcsname}%
      \fi
      #8}%
  }{%
    \def\@svsechd{%
      #6{%
       \@ifundefined{@runin@to@#1}{\@runin@to}{\csname @runin@to@#1\endcsname}%
       {\hskip#3\relax\@svsec}{#8}%
      }%
      \@ifundefined{#1mark}{\@gobble}{\csname #1mark\endcsname}{#7}%
      \addcontentsline{toc}{#1}{%
        \ifnum #2>\c@secnumdepth \else
          \protect\numberline{\csname the#1\endcsname}%
        \fi
        #8}%
    }%
  }%
  \@xsect{#5}%
}%
\def\@runin@to#1#2{#1#2}%
\REVTEX@info
  {Repairing broken definition of \string\@ssect.}%
\def\@ssect#1#2#3#4#5#6[#7]#8{%
  \@tempskipa #5\relax
  \@ifdim{\@tempskipa>\z@}{%
    \begingroup
      \interlinepenalty \@M
      #6{%
       \@ifundefined{@hangfroms@#1}{\@hangfrom}{\csname @hangfroms@#1\endcsname}%
       {\hskip#3\relax}{#8}%
      }%
      \@@par
    \endgroup
    \@ifundefined{#1smark}{\@gobble}{\csname #1smark\endcsname}{#7}%
    \addcontentsline{toc}{#1}{#8}%
  }{%
    \def\@svsechd{%
      #6{%
       \@ifundefined{@runin@tos@#1}{\@runin@to}{\csname @runin@tos@#1\endcsname}%
       {\hskip#3\relax}{#8}%
      }%
      \@ifundefined{#1smark}{\@gobble}{\csname #1smark\endcsname}{#7}%
      \addcontentsline{toc}{#1}{#8}%
    }%
  }%
  \@xsect{#5}%
}%
\def\sec@upcase#1{\relax{#1}}%
\AtBeginDocument{%
 \prepdef\endtabular{\endtabular@hook}%
 \prepdef\endarray{\endarray@hook}%
 \providecommand\array@hook{}%
 \@ifxundefined\endarray@hook{\let\endarray@hook\@empty}{}%
 \providecommand\tabular@hook{}%
 \@ifxundefined\endtabular@hook{\let\endtabular@hook\@empty}{}%
 \@ifpackageloaded{array}{%
  \let\@array@sw\@array@sw@array
 }{%
  \let\@array@sw\@array@sw@LaTeX
 }%
}%
\def\@array@sw@LaTeX{\@ifx{\\\@tabularcr}}%
\def\@array@sw@array{\@ifx{\d@llarbegin\begingroup}}%
\prepdef\@tabular{\tabular@font}%
\newcommand*\@gobble@opt@one[2][]{}%
\def\@starttoc#1{%
  \begingroup
    \toc@pre
    \makeatletter
    \@input{\jobname.#1}%
    \if@filesw
      \expandafter\newwrite\csname tf@#1\endcsname
      \immediate\openout \csname tf@#1\endcsname \jobname.#1\relax
    \fi
    \@nobreakfalse
    \toc@post
  \endgroup
}%
\def\toc@pre{}%
\def\toc@post{}%
\def\toc@@font{}%{\footnotesize\rmfamily}%
\def\@dotsep{\z@}%{5.5pt}%
\let\tocdim@section       \leftmargini
\let\tocdim@subsection    \leftmarginii
\let\tocdim@subsubsection \leftmarginiii
\let\tocdim@paragraph     \leftmarginiv
\let\tocdim@appendix      \leftmarginv
\let\tocdim@pagenum       \leftmarginvi
\def\toc@pre@auto{%
  \toc@@font
  \@tempdima\z@
  \toc@setindent\@tempdima{section}%
  \toc@setindent\@tempdima{subsection}%
  \toc@setindent\@tempdima{subsubsection}%
  \toc@setindent\@tempdima{paragraph}%
  \toc@letdimen{appendix}%
  \toc@letdimen{pagenum}%
}%
\def\toc@post@auto{%
  \if@filesw
   \begingroup
    \toc@writedimen{section}%
    \toc@writedimen{subsection}%
    \toc@writedimen{subsubsection}%
    \toc@writedimen{paragraph}%
    \toc@writedimen{appendix}%
    \toc@writedimen{pagenum}%
   \endgroup
  \fi
}%
\def\toc@setindent#1#2{%
 \csname tocdim@#2\endcsname\z@
 \@ifundefined{tocmax@#2}{\@namedef{tocmax@#2}{\z@}}{}%
 \advance#1\@nameuse{tocmax@#2}\relax
 \expandafter\edef\csname tocleft@#2\endcsname{\the#1}%
}%
\def\toc@letdimen#1{%
 \csname tocdim@#1\endcsname\z@
 \@ifundefined{tocmax@#1}{\@namedef{tocmax@#1}{\z@}}{}%
 \expandafter\let\csname tocleft@#1\expandafter\endcsname\csname tocmax@#1\endcsname
}%
\def\toc@writedimen#1{%
 \immediate\write\@auxout{%
  \gdef\expandafter\string\csname tocmax@#1\endcsname{%
   \expandafter\the\csname tocdim@#1\endcsname
  }%
 }%
}%
\def\l@@sections#1#2#3#4{%
 % #1 - superior section
 % #2 - this section
 % #3 - content, including possible \numberline
 % #4 - page number
 \expandafter\let\expandafter\dimen@temp\csname tocdim@#2\endcsname
 \setbox\z@\hbox{\ignorespaces#4}\ifdim\tocdim@pagenum<\wd\z@\tocdim@pagenum\wd\z@\fi
 \begingroup
  \everypar{}%
  \leftskip\csname tocleft@#2\endcsname\relax
  \parindent-\leftskip\advance\parindent\csname tocleft@#1\endcsname\relax
  \rightskip\tocleft@pagenum plus 1fil\relax
  \skip@\parfillskip\parfillskip\z@
  \let\numberline\numberline@@sections
  \ignorespaces#3\unskip\nobreak\hskip\skip@
  \hb@xt@\rightskip{\hfil\unhbox\z@}\hskip-\rightskip\hskip\z@skip
  \par
  \aftergroup\dimen@temp
  \expandafter
 \endgroup\the\dimen@temp\relax
}%
\def\numberline@@sections#1{%
 \leavevmode\hb@xt@-\parindent{%
  \aftergroup\dimen@temp
  \hfil
  \setbox\z@\hbox{#1.\kern\@dotsep}\dimen@\wd\z@
  \unhbox\z@
  \@ifdim{\dimen@<\dimen@temp}{\dimen@\dimen@temp}{}%
  \expandafter
 }\the\dimen@\relax
 \ignorespaces
}%
\AtBeginDocument{%
 \@ifxundefined\twocolumn@sw{\@booleanfalse\twocolumn@sw}{}%
 \@ifxundefined\multicol@sw{\@booleanfalse\multicol@sw}{}%
}%
\AtBeginDocument{%
 \@ifxundefined\footinbib@sw{\@booleanfalse\footinbib@sw}{}%
}%
\AtBeginDocument{%
 \@ifxundefined\preprintsty@sw{\@booleanfalse\preprintsty@sw}{}%
}%
\AtBeginDocument{%
 \@ifxundefined\galley@sw{\@booleanfalse\galley@sw}{}%
}%
\AtBeginDocument{%
 \@ifxundefined\showPACS@sw{\@booleanfalse\showPACS@sw}{}%
}%
\AtBeginDocument{%
 \@ifxundefined\twoside@sw{\@booleanfalse\twoside@sw}{}%
 \twoside@sw{\@twosidetrue}{\@twosidefalse}%
}%
\AtBeginDocument{%
 \@ifxundefined\draft@sw{\@booleanfalse\draft@sw}{}%
 \draft@sw{\overfullrule 5\p@}{\overfullrule\z@}%
}%
\AtBeginDocument{%
 \@ifxundefined\tightenlines@sw{\@booleanfalse\tightenlines@sw}{}%
}%
\AtBeginDocument{%
 \@ifxundefined\lengthcheck@sw{\@booleanfalse\lengthcheck@sw}{}%
}%
\AtBeginDocument{%
 \@ifxundefined\eqsecnum@sw{\@booleanfalse\eqsecnum@sw}{}%
 \eqsecnum@sw{%
  \@addtoreset{equation}{section}%
  \def\theequation@prefix{\thesection.}%
 }{}%
}%
\AtBeginDocument{%
 \@ifxundefined\byrevtex@sw{\@booleanfalse\byrevtex@sw}{}%
}%
\DeclareOption{preprint}{%
  \@booleantrue\preprintsty@sw
  \def\@pointsize{12}%
}%
\DeclareOption{manuscript}{%
 \REVTEX@warn{Document class option manuscript is obsolete; use preprint instead}%
 \ExecuteOptions{preprint}%
}%
\DeclareOption{showpacs}{%
  \@booleantrue\showPACS@sw
}%
\DeclareOption{noshowpacs}{%
  \@booleanfalse\showPACS@sw
}%
\DeclareOption{hyperref}{%
  \AtEndOfClass{%
    \let\REVTEX@tempa\@clsextension
    \let\@clsextension\@empty
    \RequirePackage{hyperref}%
    \let\@clsextension\REVTEX@tempa
  }%
  \def\@pointsize{10}%
}%
\DeclareOption{10pt}{\def\@pointsize{10}}
\DeclareOption{11pt}{\def\@pointsize{11}}
\DeclareOption{12pt}{\def\@pointsize{12}}
\DeclareOption{a4paper}
   {\setlength\paperheight {297mm}%
    \setlength\paperwidth  {210mm}}
\DeclareOption{a5paper}
   {\setlength\paperheight {210mm}%
    \setlength\paperwidth  {148mm}}
\DeclareOption{b5paper}
   {\setlength\paperheight {250mm}%
    \setlength\paperwidth  {176mm}}
\DeclareOption{letterpaper}
   {\setlength\paperheight {11in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{legalpaper}
   {\setlength\paperheight {14in}%
    \setlength\paperwidth  {8.5in}}
\DeclareOption{executivepaper}
   {\setlength\paperheight {10.5in}%
    \setlength\paperwidth  {7.25in}}
\DeclareOption{landscape}
   {\setlength\@tempdima   {\paperheight}%
    \setlength\paperheight {\paperwidth}%
    \setlength\paperwidth  {\@tempdima}}
\DeclareOption{bibnotes}{\def\author@note{\rtx@bibnote}}%
\DeclareOption{nobibnotes}{\def\author@note{\footnote}}%
\AtBeginDocument{\@ifxundefined\author@note{\def\author@note{\footnote}}{}}%
\DeclareOption{footinbib}{\@booleantrue\footinbib@sw}
\DeclareOption{nofootinbib}{\@booleanfalse\footinbib@sw}
\DeclareOption{superbib}{%
 \def\place@bibnumber{\place@bibnumber@sup}%
}%
\DeclareOption{twoside}{\@booleantrue\twoside@sw \@mparswitchfalse}
\DeclareOption{oneside}{\@booleanfalse\twoside@sw\@mparswitchtrue}
\AtBeginDocument{%
 \@ifxundefined\twocolumn@sw{%
  \@booleanfalse\twocolumn@sw
 }{}%
 \twocolumn@sw{\@twocolumntrue}{\@twocolumnfalse}%
}%
\AtEndDocument{%
 \class@enddocumenthook
}%
\AtEndOfClass{%
 \appdef\class@enddocumenthook{%
  \close@column
 }%
 \AtBeginDocument{%
  \@ifxundefined\close@column{%
   \let\close@column\close@column@default
  }{}%
 }%
}%
\DeclareOption{twocolumn}{%
  \@booleantrue\twocolumn@sw
  \@booleantrue\multicol@sw
}%
\DeclareOption{twocolumngrid}{%
  \@booleantrue\twocolumn@sw
  \@booleanfalse\multicol@sw
}%
\DeclareOption{onecolumn}{%
  \@booleanfalse\twocolumn@sw
}%
\DeclareOption{galley}{%
  \ExecuteOptions{onecolumn}%
  \@booleantrue\galley@sw
  \appdef\set@pica@hook{%
    \advance\textwidth-\columnsep
    \textwidth.5\textwidth
  }%
}%
\DeclareOption{tightenlines}{\@booleantrue\tightenlines@sw}
\DeclareOption{lengthcheck}{\@booleantrue\lengthcheck@sw}
\DeclareOption{draft}{%
 \@booleantrue\draft@sw
}%
\DeclareOption{final}{%
 \@booleanfalse\draft@sw
}%
\DeclareOption{eqsecnum}{\@booleantrue\eqsecnum@sw}
\DeclareOption{secnumarabic}{\@booleantrue\secnumarabic@sw}
\DeclareOption{fleqn}{\input{fleqn.clo}}
\DeclareOption{floats}{\@booleantrue\floats@sw}
\DeclareOption{endfloats}{\@booleanfalse\floats@sw}
\AtBeginDocument{%
 \@ifxundefined\floats@sw{\@booleantrue\floats@sw}{}%
}%
\DeclareOption{titlepage}{\@booleantrue\titlepage@sw}
\DeclareOption{notitlepage}{\@booleanfalse\titlepage@sw}
\appdef\defaults@hook{%
 \@ifxundefined\titlepage@sw{\@booleanfalse\titlepage@sw}{}%
}%
\def\change@society#1{%
 \def\@tempa{#1}%
 \@ifxundefined\@society{}{%
  \@ifx\@tempa\@society{}{%
   \REVTEX@warn{Conflicting society \@tempa<>\@society.}%
  }%
 }%
 \let\@society\@tempa
}%
\DeclareOption{aps}{\change@society{aps}\def\@journal{}}%
\DeclareOption{pra}{\change@society{aps}\def\@journal{pra}}%
\DeclareOption{prb}{\change@society{aps}\def\@journal{prb}}%
\DeclareOption{prc}{\change@society{aps}\def\@journal{prc}}%
\DeclareOption{prd}{\change@society{aps}\def\@journal{prd}}%
\DeclareOption{pre}{\change@society{aps}\def\@journal{pre}}%
\DeclareOption{prl}{\change@society{aps}\def\@journal{prl}}%
\DeclareOption{prstab}{\change@society{aps}\def\@journal{prstab}}%
\DeclareOption{rmp}{\change@society{aps}\def\@journal{rmp}}%
\DeclareOption{amsfonts}{%
  \def\REVTEX@amsfonts{\RequirePackage{amsfonts}}%
}%
\DeclareOption{noamsfonts}{%
  \let\REVTEX@amsfonts\@empty
}%
\AtEndOfClass{\@ifxundefined\REVTEX@amsfonts{}{\REVTEX@amsfonts}}%
\DeclareOption{amssymb}{%
  \def\REVTEX@amssymb{\RequirePackage{amssymb}}%
}%
\DeclareOption{noamssymb}{%
  \let\REVTEX@amssymb\@empty
}%
\AtEndOfClass{\@ifxundefined\REVTEX@amssymb{}{\REVTEX@amssymb}}%
\DeclareOption{groupedaddress}{\clo@groupedaddress}%
\def\clo@groupedaddress{%
 \@booleantrue\groupauthors@sw
 \@booleantrue\@affils@sw
 \@booleanfalse\runinaddress@sw
}%
\AtBeginDocument{%
 \@ifxundefined\groupauthors@sw{\clo@groupedaddress}{}%
}%
\DeclareOption{unsortedaddress}{\clo@unsortedaddress}%
\def\clo@unsortedaddress{%
 \@booleantrue\groupauthors@sw
 \@booleanfalse\@affils@sw
 \@booleanfalse\runinaddress@sw
}%
\DeclareOption{runinaddress}{\clo@runinaddress}%
\def\clo@runinaddress{%
 \@booleantrue\groupauthors@sw
 \@booleantrue\@affils@sw
 \@booleantrue\runinaddress@sw
}%
\DeclareOption{superscriptaddress}{\clo@superscriptaddress}%
\def\clo@superscriptaddress{%
 \@booleanfalse\groupauthors@sw
 \@booleantrue\@affils@sw
 \@booleanfalse\runinaddress@sw
}%
\DeclareOption{byrevtex}{\@booleantrue\byrevtex@sw}%
\DeclareOption*{%
 \IfFileExists{\CurrentOption\REVTEX@society@ext}{%
  \expandafter\change@society\expandafter{\CurrentOption}%
 }{%
  \OptionNotUsed
 }%
}%
\def\REVTEX@society@ext{.rtx}%
\AtEndOfClass{%
 \@ifxundefined\@society{%
  \REVTEX@warn{No Society specified, using generic APS style}%
  \def\@society{aps}\def\@journal{}%
 }{}%
 \expandafter\input\expandafter{\@society\REVTEX@society@ext}%
}%
\ExecuteOptions{10pt}
\ProcessOptions\relax
\def\ps@headings@twoside{%
    \let\@mkboth\markboth
    \let\@oddfoot\@empty
    \let\@evenfoot\@empty
    \def\@evenhead{\thepage\hfil \slshape \leftmark}%
    \def\@oddhead{{\slshape \rightmark}\hfil \thepage}%
    \def\sectionmark##1{%
      \markboth{%
       \MakeTextUppercase{%
        \@ifnum{\c@secnumdepth >\z@}{\thesection\hskip 1em\relax}{}%
         ##1%
       }%
       }{}%
    }%
    \def\subsectionmark##1{%
      \markright {%
        \@ifnum{\c@secnumdepth >\@ne}{\thesubsection\hskip 1em\relax}{}%
         ##1%
      }%
    }%
}%
\def\ps@headings@oneside{%
    \let\@mkboth\markboth
    \let\@oddfoot\@empty
    \let\@evenfoot\@empty
    \def\@oddhead{{\slshape \rightmark}\hfil\thepage}%
    \def\sectionmark##1{%
     \markright {%
      \MakeTextUppercase{%
       \@ifnum{\c@secnumdepth >\z@}{\thesection\hskip 1em\relax}{}%
       ##1%
      }%
     }%
    }%
}%
\AtBeginDocument{%
 \twoside@sw{%
  \let\ps@headings\ps@headings@twoside
 }{%
  \let\ps@headings\ps@headings@oneside
 }%
}%
\def\ps@myheadings{%
  \let\@mkboth\@gobbletwo
  \def\@oddhead{{\slshape\rightmark}\hfil \thepage}%
  \let\@oddfoot\@empty
  \def\@evenhead{\thepage\hfil\slshape\leftmark}%
  \let\@evenfoot\@empty
  \let\sectionmark\@gobble
  \let\subsectionmark\@gobble
}%
\def\ps@article{%
  \def\@evenhead{\thepage\hfil{\leftmark}}%
  \def\@oddhead{{\rightmark}\hfil \thepage}%
  \def\@oddfoot{\byrevtex@sw{\hfil Typeset by \revtex\hfil}{}}%
  \let\@evenfoot\@oddfoot
  \let\@mkboth\@gobbletwo
  \let\sectionmark\@gobble
  \let\subsectionmark\@gobble
}%
\def\lastpage@putlabel{%
 \if@filesw
  \begingroup
    \advance\c@page\m@ne
    \immediate\write\@auxout{\string\newlabel{LastPage}{{}{\thepage}}}%
  \endgroup
 \fi
}%
\AtEndDocument{%
 \clearpage\lastpage@putlabel
}%
\RequirePackage{url}%
\renewcommand\title{\@dblarg\REVTEX@title}%
\def\REVTEX@title[#1]#2{%
  \def\@title{#2}%
  \def\@shorttitle{#1}%
}%
\def\@title{\REVTEX@warn{No title}}%
\def\@shorttitle{}%
\newcounter{affil}
\renewcommand\author[1]{%
  \@author@sw{%
    \@addaffils
  }{}%
  \move@author
  \def\@author{{}{}{#1}{}}%
  \let\thanks\thanks@author
}%
\def\move@author{%
 \@author@sw{}{%
  \collaboration@sw{%
    \@booleanfalse\collaboration@sw
    \expandafter\appdef
    \expandafter\@authors@curr
    \expandafter{%
    \expandafter \@collaboration@opr
                 \@author
                }%
  }{%
    \expandafter\appdef
    \expandafter\@authors@curr
    \expandafter{%
    \expandafter \@author@opr
                 \@author
                }%
  }%
 }%
}%
\@booleanfalse\collaboration@sw
\def\collaboration#1{%
 \groupauthors@sw{%
  \REVTEX@warn{%
   \string\collaboration\space not allowed with groupaddress option.
   Use superscriptaddress instead
  }%
 }{%
  \@booleantrue\collaboration@sw
 }%
 \author{#1}%
}%
\def\@collaboration@opr#1#2#3#4{%
  \par
  {%
   \let\@affilID@temp\@empty
   \@author@opr{#1}{#2}{(\ignorespaces#3\unskip)}{#4}%
  }%
  \par
}%
\def\@collaboration@count{%
 \appdef@eval\@tempa{\the\@tempcnta}%
 \@tempcnta\z@
 \@gobblefour
}%
\def\@collaboration@present#1#2#3#4{%
 \par
 \@author@opr{#1}{#2}{(\ignorespaces#3\unskip)}{#4}\par
 \smallskip
 \reset@authorcount
}%
\def\@author@cleared{{}{}{}{}}%
\def\@author@init{\let\@author\@author@cleared\let\thanks\thanks@warning}%
\def\@author@sw{\@ifx{\@author\@author@cleared}}%
\AtBeginDocument{%
 \@author@init
}%
\def\@author@join@i#1#2#3#4#5{%
  \def\@author{{\@join{\@separator}{#1}{#5}}{#2}{#3}{#4}}%
}%
\def\@author@join@ii#1#2#3#4#5{%
  \def\@author{{#1}{\@join{\@separator}{#2}{#5}}{#3}{#4}}%
}%
\def\@author@join@iii#1#2#3#4#5{%
  \def\@author{{#1}{#2}{\@join{\@separator}{#3}{#5}}{#4}}%
}%
\def\@author@join@iv#1#2#3#4#5{%
  \def\@author{{#1}{#2}{#3}{\@join{\@separator}{#4}{#5}}}%
}%
\def\@author@finish{%
 \@author@sw{}{%
  \move@author
  \REVTEX@warn{No affiliation specified for \@author}%
 }%
 \@addaffils
}%
\def\@join#1#2#3{%
  \@ifempty{#2}{#3}{#2#1#3}%
}%
\def\@separator{;\space}%
\let\surname\@firstofone
\let\firstname\@firstofone
\let\@authors@curr\@empty
\renewcommand\and{\REVTEX@err{\protect\and\space is not supported}}
\def\cat@comma@active{\catcode`\,\active}%
{\cat@comma@active\gdef,{\active@comma}}%
\def\active@comma{,\penalty-400\relax}%
\newcommand*\affiliation[1][]{\begingroup\cat@comma@active\@affiliation{#1}}%
\def\@affiliation#1#2{%
  \endgroup
  \move@author\@author@init
  \@booleanfalse\temp@sw
  \let\@tempd\@empty
  \@affils@sw{%
    \def\@tempa{#2}\let\@affil@opr\@affil@match\@affil@list
  }{}\temp@sw
  {%
   \def\@tempa{#1}%
   \@ifx{\@tempd\@tempa}{}{%
    \REVTEX@warn{%
     New footnote specified for \string\affiliation{#2}.
     Please put it on the *first* instance thereof.
    }%
   }%
  }%
  {%
    \stepcounter{affil}%
    \edef\@tempc{\theaffil}%
    \appdef     \@affil@list{\@affil@opr}%
    \appdef@eval\@affil@list\@tempc
    \appdef     \@affil@list{{#1}{#2}}%
  }%
  \appdef@eval\@affilID@curr\@tempc
}%
\newcommand\noaffiliation{%
  \move@author\@author@init
  \edef\@affilID@curr{\@affilID@curr}%
}%
\let\@authors\@empty
\def\@addaffils{%
    \appdef     \@authors{\@affilID@opr}%
    \appdef@eval\@authors{\@affilID@curr}%
    \appdef@e   \@authors{\@authors@curr}%
    \@booleanfalse\temp@sw
    \let\affils@opr\affils@match\affils@list
    \temp@sw{}{%
      \appdef     \affils@list{\affils@opr}%
      \appdef@eval\affils@list{\@affilID@curr}%
    }%
    \let\@authors@curr\@empty
    \let\@affilID@curr\@empty
}%
\let\affils@list\@empty
\let\@affil@list\@empty
\let\@affilID@curr\@empty
\def\@affil@match#1#2#3{%
  \def\@tempb{#3}%
  \@ifx{\@tempa\@tempb}{%
    \@booleantrue\temp@sw
    \def\@tempd{#2}%
    \def\@tempc{#1}%
  }{}%
}%
\def\@tempaffil@script#1#2{%
 ,\,#1%
}%
\def\@tempaffil@group#1#2{%
 \def\@tempc{#1}%
 \def\@tempb{#2}%
 \begingroup
  \let\@affil@opr\@affil@footnote\@affil@list
 \endgroup
}%
\def\@affil@footnote#1#2#3{%
  \def\@tempa{#1}%
  \@ifx{\@tempa\@tempc}{%
   \comma@space
   \footnote{\@tempb#3}%
   \footnote{#2}%
  }{}%
}%
\def\@affil@group#1#2#3{%
  \@ifnum{\count@=#1\relax}{%
    #3%
    \@ifempty{#2}{}{\footnote{#2}}%
    \after@address
    \advance\@tempcnta\m@ne
  }{}%
}%
\newcommand*\email[2][Electronic address: ]{%
  \expandafter\@author@join@ii\@author{#1\href{mailto:#2}{#2}}%
}%
\newcommand\homepage{%
 \begingroup
  \@makeother\~%
  \@makeother\%%
  \@homepage
}%
\newcommand*\@homepage[2][URL:~]{%
 \endgroup
 \expandafter\@author@join@ii\@author{#1\url{#2}}%
}%
\AtBeginDocument{%
 \providecommand\url{\texttt}%
 \providecommand\href[2]{\url{#2}}%
}%
\def\thanks@warning{%
 \REVTEX@warn{\string\thanks\space is valid only following \string\author.^^JPlease move it.}%
 \@gobble
}%
\let\thanks\thanks@warning
\def\thanks@author{%
  \expandafter\@author@join@ii\@author
}%
\newcommand*\altaffiliation[2][]{%
  \expandafter\@author@join@ii\@author{#1#2}%
}%
\def\@listand{\@ifnum{\@tempcnta=\tw@}{\andname}{}\space}%
\def\@listcomma{\@ifnum{\@tempcnta>\@ne}{,}{}}%
\def\@author@present@group#1#2#3#4{%
  \gdef\comma@space{\gdef\comma@space{\textsuperscript{,\,}}}%
  \doauthor{#1}{#2}{#3}{#4}%
  \advance\@tempcnta\m@ne
}%
\def\@author@present@script#1#2#3#4{%
  \gdef\comma@space{\textsuperscript{,\,}}%
  \doauthor{#1}{#2}{#3}{%
   \textsuperscript{\expandafter\@affilcomma\@affilID@temp\relax\relax}%
   #4%
  }%
  \advance\@tempcnta\m@ne
}%
\def\@affilcomma#1#2{%
  #1%
  \@ifx{\relax#2}{}{%
    ,\,\@affilcomma
  }%
  #2%
}%
\def\@doauthor@count{\advance\@tempcnta\@ne\@gobblefour}%
\def\doauthor#1#2#3#4{%
  \ignorespaces#3\unskip\@listcomma
  \@ifx{\@affilID@temp\@empty}{}{#4}%
  \@ifempty{#1#2}{}{%
   \@ifx{\@affilID@temp\@empty}{}{\comma@space}%
   \author@note{%
    % \ttfamily
    #1%
    \@ifempty{#1}{}{\@ifempty{#2}{}{, }}%
    #2%
   }%
  }%
  \space\@listand
}%#4
\def\@affilID@def#1{%
 \def\@affilID@temp{#1}%
}%
\let\@affilID@temp\@empty
\let\@affilID@opr\@affilID@def
\def\@affilID@count#1{%
  \def\@tempb{#1}%
  \@ifx{\@tempa\@tempb}{%
    \let\@author@opr\@doauthor@count
  }{%
    \let\@author@opr\@gobblefour
  }%
}%
\def\@affilID@match#1{%
 \def\@tempb{#1}%
 \@ifx{\@tempa\@tempb}{%
   \let\@author@opr\@author@present
 }{%
   \let\@author@opr\@gobblefour
 }%
}%
\def\@thanks{\comma@space\footnote}
\def\@affil@marked#1#2#3{%
  \par
  \begingroup
   \frontmatter@affiliationfont
   \textsuperscript{\normalfont#1}%
   #3\par
  \endgroup
}%
\def\affils@match#1{%
  \def\@tempa{#1}%
  \@ifx{\@affilID@curr\@tempa}{%
    \@booleantrue\temp@sw
  }{}%
}%
\def\do@affil@fromgroup#1#2{%
  \@ifx{\relax#2}{}{%
    \count@#2\relax
    #1\do@affil@fromgroup#1%
  }%
}%
\def\keywords#1{%
  \gdef\@keywords{#1}}
\def\@keywords{}
\renewcommand*\date[2][\Dated@name]{\def\@date{#1#2}}%
\def\@date{}%
\newcommand*\received[2][\Received@name]{\def\@received{#1#2}}%
\def\@received{}%
\newcommand*\revised[2][\Revised@name]{\def\@revised{#1#2}}%
\def\@revised{}%
\newcommand*\accepted[2][\Accepted@name]{\def\@accepted{#1#2}}%
\def\@accepted{}%
\newcommand*\published[2][\Published@name]{\def\@published{#1#2}}%
\def\@published{}%
\def\pacs#1{%
  \REVTEX@aftermaketitle{\pacs}%
  \gdef\@pacs{#1}%
}%
\def\@pacs{}%{\REVTEX@warn{No PACS numbers}}%
\def\preprint#1{\gappdef\@preprint{\preprint{#1}}}%
\let\@preprint\@empty
\newbox\absbox
\newenvironment{abstract}{%
  \REVTEX@aftermaketitle{\begin{abstract}}%
  \lrbox\absbox
   \minipage{\frontmatter@abstractwidth}%
    \frontmatter@abstractheading
    \frontmatter@abstractfont
}{%
   \endminipage
  \endlrbox
}%
\def\frontmatter@abstractheading{%
 \begin{center}\large\abstractname\end{center}%
 \vspace{-1pc}%
}%
\def\frontmatter@abstractfont{}%
\def\frontmatter@abstractwidth{\textwidth}
\def\frontmatter@abstract@produce{%
  \addvspace\frontmatter@preabstractspace
  \leavevmode\box\absbox\par
  \addvspace\frontmatter@postabstractspace
}%
\def\frontmatter@preabstractspace{.5\baselineskip}
\def\frontmatter@postabstractspace{.5\baselineskip}
\newenvironment{titlepage}{%
      \let\wastwocol@sw\twocolumn@sw
      \wastwocol@sw{\onecolumn}{\newpage}%
      \thispagestyle{titlepage}%
      \setcounter{page}\@ne
}{%
     \wastwocol@sw{\twocolumn}{\newpage}%
     \twoside@sw{}{%
        \setcounter{page}\@ne
     }%
}%
\def\maketitle{%
  \@author@finish
  \title@column\titleblock@produce
  \suppressfloats[t]%
  \let\and\relax
  \let\affiliation\@gobble@opt@one
  \let\author\@gobble
  \@author@init
  \let\@authors\@empty
  \let\@authors@curr\@empty
  \let\@affil@list\@empty
  \let\keywords\@gobble
  \let\@keywords\@empty
  \let\email\@gobble
  \let\@address\@empty
  \let\maketitle\relax
  \let\thanks\@gobble
  \titlepage@sw{%
    \clearpage
  }{}%
}%
\def\titleblock@produce{%
 \begingroup
  \groupauthors@sw{%
   \let\@author@present\@author@present@group
   \let\frontmatter@author@produce\frontmatter@author@produce@group
   \let\@tempaffil\@tempaffil@group
  }{%
   \let\@author@present\@author@present@script
   \let\frontmatter@author@produce\frontmatter@author@produce@script
   \let\@tempaffil\@tempaffil@script
  }%
  \let\footnote\REVTEX@footnote
  \let\thefootnote\frontmatter@thefootnote
  \global\c@footnote\@ne % FIXME: what if the title had a footnote?
  \let\@makefnmark\frontmatter@makefnmark
  \frontmatter@setup
  \thispagestyle{titlepage}\label{FirstPage}%
  \frontmatter@title@produce
  \frontmatter@author@produce
  \frontmatter@RRAPformat{%
   \expandafter\produce@RRAP\expandafter{\@date}%
   \expandafter\produce@RRAP\expandafter{\@received}%
   \expandafter\produce@RRAP\expandafter{\@revised}%
   \expandafter\produce@RRAP\expandafter{\@accepted}%
   \expandafter\produce@RRAP\expandafter{\@published}%
  }%
  \frontmatter@abstract@produce
  \@pacs@produce\@pacs
  \@keywords@produce\@keywords
  \par
  \addvspace\frontmatter@finalspace
  \par
 \endgroup
}%
\def\frontmatter@title@produce{%
 \begingroup\@title\par\endgroup
 \addvspace{\baselineskip}%
}%
\def\frontmatter@author@produce@group{%
  \begingroup
    \frontmatter@authorformat
    \let\@collaboration@opr\@gobblefour
    \let\affils@opr\affils@present\affils@list
    \par
  \endgroup
}%
\def\frontmatter@author@produce@script{%
  \begingroup
    \frontmatter@authorformat
    \let\@tempa\@empty
    \@tempcnta\z@
    \let\@affilID@opr\@gobble
    \let\@author@opr\@doauthor@count
    \let\@collaboration@opr\@collaboration@count
    \@authors
    \appdef@eval\@tempa{\the\@tempcnta}%
    \reset@authorcount
    \let\@affilID@opr\@affilID@def
    \let\@author@opr\@author@present
    \let\@collaboration@opr\@collaboration@present
    \@authors\unskip\unskip
    \par
    \frontmatter@above@affiliation
    \let\@affil@opr\@affil@marked\@affil@list
    \par
  \endgroup
}%
\def\affils@present#1{%
  \def\@tempa{#1}%
  \let\@affilID@opr\@affilID@count
  \begingroup
   \@tempcnta\z@
   \@authors
   \@ifnum{\@tempcnta=\z@}{}{%
    \begingroup
     \@ifnum{\@tempcnta=\tw@}{%
        \let\@listcomma\relax
     }{}%
     \let\@affilID@opr\@affilID@match\@authors
    \endgroup
    \begingroup
     \par
     \frontmatter@above@affiliation
     \frontmatter@affiliationfont
     \let\\\frontmatter@addressnewline
     \runinaddress@sw{%
       \@tempcnta\z@
       \@tfor\@tempa:=#1\do{\advance\@tempcnta\@ne}%
       \@ifnum{\@tempcnta=\tw@}{%
          \let\@listcomma\relax
       }{}%
       \def\after@address{\@listcomma\ \@listand}%
     }{%
      \tightenlines@sw{}{%
       \parskip\z@ %space between affilations in grouped mode
      }%
      \let\after@address\par
     }%
     \let\@affil@opr\@affil@group
     \do@affil@fromgroup\@affil@list#1\relax
    \endgroup
   }%
   \par
  \endgroup
 \smallskip
}%
\def\reset@authorcount{%
  \expandafter\@reset@ac\@tempa0\@reset@ac
  \@ifnum{\@tempcnta=\tw@}{%
    \let\@listcomma\relax
  }{%
    \def\@listcomma{\@ifnum{\@tempcnta>\@ne}{,}{}}%
  }%
}%
\def\@reset@ac#1#2\@reset@ac{%
  \@tempcnta#1%
  \def\@tempa{#2}%
}%
\def\@pacs@produce#1{%
 \showPACS@sw{\parbox[t]{\frontmatter@abstractwidth}{\@pacs@name#1}}{}%
}%
\def\@keywords@produce#1{%
}%
\def\frontmatter@thefootnote{\@fnsymbol\c@footnote}%
\def\frontmatter@makefnmark{\@textsuperscript{\normalfont\@thefnmark}}%
\long\def\frontmatter@makefntext#1{%
 \parindent 1em
 \noindent
 \@textsuperscript{\normalfont\@thefnmark}#1%
}%
\def\frontmatter@setup{}%
\def\frontmatter@RRAPformat#1{\begingroup#1\par\endgroup}%
\def\punct@RRAP{; }%
\def\produce@RRAP#1{%
  \@ifempty{#1}{}{%
   \@ifvmode{\leavevmode}{\unskip\punct@RRAP\ignorespaces}%
   #1%
  }%
}%
\def\frontmatter@authorformat{}%
\def\frontmatter@above@affiliation{}%
\def\frontmatter@affiliationfont{}%
\def\frontmatter@finalspace{2\baselineskip}
\def\frontmatter@addressnewline{%
  \@ifhmode{\skip@\lastskip\unskip\unpenalty\break\hskip\skip@}{}%
  % was: \vskip-.5ex
}%
\def\REVTEX@aftermaketitle#1{%
  \@ifx{\maketitle\relax}{%
    \REVTEX@err{\protect#1 must be used before \protect\maketitle}%
  }{}%
}%
\def\ps@titlepage{\ps@empty}%
\def\volumeyear#1{\gdef\@volumeyear{#1}}%
\def\@volumeyear{}%
\def\volumenumber#1{\gdef\@volumenumber{#1}}%
\def\@volumenumber{}%
\def\issuenumber#1{\gdef\@issuenumber{#1}}%
\def\@issuenumber{}%
\def\eid#1{\gdef\@eid{#1}}%
\def\@eid{}%
\def\startpage#1{\gdef\@startpage{#1}\c@page#1\relax}%
\def\@startpage{\pageref{FirstPage}}%
\def\endpage#1{\gdef\@endpage{#1}}%
\def\@endpage{\pageref{LastPage}}%
\def\print@toc#1{%
 \begingroup
  \expandafter\section
  \expandafter*%
  \expandafter{%
              \csname#1name\endcsname
              }%
  \let\appendix\appendix@toc
  \@starttoc{#1}%
 \endgroup
}%
\def\appendix@toc{}%
\def\labelenumi{\theenumi.}
\def\theenumi{\arabic{enumi}}
\def\labelenumii{(\theenumii)}
\def\theenumii{\alph{enumii}}
\def\p@enumii{\theenumi}
\def\labelenumiii{\theenumiii.}
\def\theenumiii{\roman{enumiii}}
\def\p@enumiii{\theenumi(\theenumii)}
\def\labelenumiv{\theenumiv.}
\def\theenumiv{\Alph{enumiv}}
\def\p@enumiv{\p@enumiii\theenumiii}
\def\labelitemi{\textbullet}
\def\labelitemii{\normalfont\bseries\textendash}
\def\labelitemiii{\textasteriskcentered}
\def\labelitemiv{\textperiodcentered}
\pagenumbering{arabic}
\setcounter{topnumber}{2}
\def\topfraction{.7}
\setcounter{bottomnumber}{1}
\def\bottomfraction{.3}
\setcounter{totalnumber}{3}
\def\textfraction{.2}
\def\floatpagefraction{.5}
\setcounter{dbltopnumber}{2}
\def\dbltopfraction{.7}
\def\dblfloatpagefraction{.5}
\newenvironment{verse}{%
  \let\\=\@centercr
  \list{}{%
    \itemsep\z@ \itemindent -1.5em\listparindent \itemindent
    \rightmargin\leftmargin\advance\leftmargin 1.5em}\item[]%
}{%
 \endlist
}%
\newenvironment{quotation}{%
  \list{}{%
    \listparindent 1.5em
    \itemindent\listparindent
    \rightmargin\leftmargin \parsep \z@ \@plus\p@}\item[]%
}{%
 \endlist
}%
\newenvironment{quote}{%
  \list{}{%
  \rightmargin\leftmargin}\item[]%
}{%
 \endlist
}%
\def\descriptionlabel#1{%
  \hspace\labelsep \normalfont\bfseries #1%
}%
\newenvironment{description}{%
 \list{}{%
    \labelwidth\z@ \itemindent-\leftmargin
    \let\makelabel\descriptionlabel
 }%
}{%
 \endlist
}%
\newcounter{part}
\newcounter{section}
\newcounter{subsection}[section]
\newcounter{subsubsection}[subsection]
\newcounter{paragraph}[subsubsection]
\newcounter{subparagraph}[paragraph]
\def\thepart          {\Roman{part}} %
\def\thesection       {\Roman{section}}
\def\p@section        {}
\def\thesubsection    {\Alph{subsection}}
\def\p@subsection     {\thesection\,}
\def\thesubsubsection {\arabic{subsubsection}}
\def\p@subsubsection  {\thesection\,\thesubsection\,}
\def\theparagraph     {\alph{paragraph}}
\def\thesubparagraph  {\theparagraph.\arabic{subparagraph}}
\newenvironment{acknowledgments}{%
 \section*{\acknowledgmentsname}%
}{%
 \par
}%
\def\part{\par
   \addvspace{4ex}
   \@afterindentfalse
   \secdef\@part\@spart}
\def\@part[#1]#2{%
 \@ifnum{\c@secnumdepth >\m@ne}{%
        \refstepcounter{part}
        \addcontentsline{toc}{part}{\thepart\hspace{1em}#1}%
 }{%
      \addcontentsline{toc}{part}{#1}%
 }%
 \begingroup
    \parindent \z@ \raggedright
    \interlinepenalty\@M
    \@ifnum{\c@secnumdepth >\m@ne}{%
      \Large \bf \partname~\thepart%
      \par\nobreak
    }{}%
    \huge \bf
    #2%
    \markboth{}{}\par
 \endgroup
   \nobreak
   \vskip 3ex
   \@afterheading
}%
\def\@spart#1{{\parindent \z@ \raggedright
    \interlinepenalty\@M
    \huge \bf
    #1\par}
    \nobreak
    \vskip 3ex
    \@afterheading}
\def\section{%
  \@startsection
    {section}%
    {1}%
    {\z@}%
    {0.8cm \@plus1ex \@minus .2ex}%
    {0.5cm \@plus1ex \@minus.2ex}%
    {%
      \normalfont\small\bfseries
      \centering
      \MakeTextUppercase
    }%
}%
\def\@sectioncntformat#1{\csname the#1\endcsname.\quad}%
\def\subsection{%
  \@startsection
    {subsection}%
    {2}%
    {\z@}%
    {.8cm \@plus1ex \@minus .2ex}%
    {.5cm \@plus1ex \@minus.2ex}%
    {\normalfont\small\bfseries\centering}}
\def\subsubsection{%
  \@startsection
    {subsubsection}%
    {3}%
    {\z@}%
    {.8cm \@plus1ex \@minus .2ex}%
    {.5cm \@plus1ex \@minus.2ex}%
    {\normalfont\small\itshape\centering}}
\def\paragraph{%
  \@startsection
    {paragraph}%
    {4}%
    {\parindent}%
    {\z@}%
    {-1em}%
    {\normalfont\normalsize\itshape}%
}%
\def\subparagraph{%
  \@startsection
    {subparagraph}%
    {5}%
    {\parindent}%
    {3.25ex \@plus1ex \@minus .2ex}%
    {-1em}%
    {\normalfont\normalsize\bseries}%
}%
\def\theequation{%
 \theequation@prefix\arabic{equation}%
}%
\def\theequation@prefix{}%
\AtEndOfClass{%
 \expandafter\input\expandafter{\@pointsize pt\REVTEX@society@ext}%
}%
\setcounter{secnumdepth}{4}
\lineskip 1pt
\normallineskip 1pt
\def\baselinestretch{1}%
\@lowpenalty   51
\@medpenalty  151
\@highpenalty 301
\@beginparpenalty -\@lowpenalty
\@endparpenalty   -\@lowpenalty
\@itempenalty     -\@lowpenalty
\arraycolsep 3pt
\tabcolsep 2pt
\arrayrulewidth .4pt
\doublerulesep 2pt
\skip\@mpfootins = \skip\footins
\fboxsep  = 3.0pt
\fboxrule = 0.4pt
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{0\p@}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{#1: #2}%
  \@ifdim{\wd\@tempboxa >\hsize}{%
    #1: #2\par
  }{%
    \global \@minipagefalse
    \hb@xt@\hsize{\unhbox\@tempboxa\hfil}%
  }%
  \vskip\belowcaptionskip
}%
\newcounter{figure}
\renewcommand \thefigure {\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\AtBeginDocument{\do@if@floats{figure}{.fgx}}%
\appdef\class@enddocumenthook{%
 \printfigures
}%
\newcommand\printfigures{\print@float{figure}}%
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@Roman\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\def\tabular@font{\small}%
\def\squeezetable{\def\tabular@font{\scriptsize}}%
\AtBeginDocument{%
 \floats@sw{}{%
  \let@environment{longtable@float}{longtable}%
  \let@environment{longtable}{longtable@write}%
 }%
 \do@if@floats{table}{.tbx}%
}%
\appdef\class@enddocumenthook{%
 \printtables
}%
\newenvironment{longtable@write}{%
 \write@@float{longtable}{table}%
}{%
 \endwrite@float
}%
\newcommand\printtables{%
 \let@environment{longtable}{longtable@float}%
 \print@float{table}%
}%
\def\tabular@hook{}%
\def\endtabular@hook{}%
\AtBeginDocument{%
 \@ifpackageloaded{dcolumn}{%
  \newcolumntype{d}{D{.}{.}{-1}}%
 }{}%
}%
\def\toprule{\hline\hline}%
\def\colrule{\hline}%
\def\botrule{\hline\hline}%
\newcommand\@makefntext[1]{%
  \def\baselinestretch{1}%
  \reset@font\footnotesize
  \parindent 1em%
  \noindent
  \hb@xt@1.8em{\hss\@makefnmark}%
  #1\par
}%
\newenvironment{thebibliography}{}{}%
\let\@listi\@empty
\RequirePackage{natbib}%
\let\NAT@thebibliography\thebibliography
\let\NAT@endthebibliography\endthebibliography
\let\LATEX@bibliographystyle\bibliographystyle
\def\bibliographystyle{\def\@bibstyle}%
\def\bibsection{\section*{\refname}}%
\def\REV@bibsetup#1{\NAT@bibsetup{\ref{LastBibItem}}}%
\let\bibpreamble\@empty
\def\bibfont{%
 \bibsep\z@\relax
 \let\NAT@bibsetup\@bibsetup
 \let\@bibsetup\REV@bibsetup
 \@ifxundefined\place@bibnumber{%
  \let\@biblabel\place@bibnumber@inl
 }{%
  \let\@biblabel\place@bibnumber
 }%
 \def\newblock{\ }%
}%
\let\NAT@citeyear\citeyear
\def\NAT@onlinecite#1{%
 \begingroup
  \let\@cite\NAT@citenum
  \let\NAT@mbox\mbox
  \let\citeyear\NAT@citeyear
  \def\NAT@space{ }%
  \cite{#1}%
 \endgroup
}%
\AtBeginDocument{%
 \@ifxundefined\onlinecite{\def\onlinecite{\citealp}}{}%
 \@ifxundefined\textcite{\def\textcite{\citet}}{}%
}%
\renewenvironment{thebibliography}[1]{%
 \NAT@thebibliography{#1}%
 \@bibnotes
}{%
 \@endnotesinbib
 \edef\@currentlabel{\arabic{NAT@ctr}}%
 \label{LastBibItem}%
 \NAT@endthebibliography
}%
\def\place@bibnumber@sup#1{\textsuperscript{#1}}%
\def\place@bibnumber@inl#1{[#1]}%
\def\rtx@bibnote#1{%
  \stepcounter\@mpfn
  \protected@xdef\@thefnmark{\thempfn}%
  \protected@xdef\@@thefnmark{\textsuperscript{\normalfont\thempfn}}%
  \@footnotemark
  \expandafter\g@addto@macro
  \expandafter\@bibnotes
  \expandafter{%
  \expandafter \item
  \expandafter [\@@thefnmark]#1\par
              }%
}%
\let\@bibnotes\@empty
\def\endnote{%
 \@ifnextchar[{%
  \@xendnote
 }{%
  \stepcounter{footnote}%
  \protected@xdef\@thefnmark{\thefootnote}%
  \@footnotemark
  \@endnotetext
 }%
}%
\def\@xendnote[#1]{%
   \begingroup
     \c@footnote#1\relax
     \unrestored@protected@xdef\@thefnmark{\thefootnote}%
   \endgroup
   \@footnotemark
   \@endnotetext
}%
\def\REVTEX@endnote@ext{.end}%
\long\def\@endnotetext#1{%
  \@ifx{\@endnoteout\@undefined}{%
    \newwrite\@endnoteout
    \gdef\endnote@stream{\jobname\REVTEX@endnote@ext}%
    \immediate\openout\@endnoteout\endnote@stream\relax
  }{}%
  \begingroup
    \let\label\relax \let\index\relax \let\glossary\relax
    \let\cite\relax \let\ref\relax \let\pageref\relax
    \let\\\relax
    \let\protect\noexpand
    \newlinechar`\^^M%
    \immediate\write\@endnoteout{\noexpand\REVdoendnote{\@thefnmark}{#1}}%
  \endgroup
}%
\def\endnotetext{%
  \@ifnextchar[{%
   \@xendnotenext
  }{%
   \protected@xdef\@thefnmark{\thefootnote}\@endnotetext
  }%
}%
\def\@xendnotenext[#1]{%
  \begingroup
     \c@footnote#1\relax
     \unrestored@protected@xdef\@thefnmark{\thefootnote}%
  \endgroup
  \@endnotetext
}%
\def\printendnotes{%
 \@ifx{\@endnoteout\@undefined}{}{%
  \begingroup
   \c@secnumdepth-\maxdimen
   \section{%
    \notesname
   }%
   \immediate\closeout\@endnoteout
   \global\let\@endnoteout\@undefined
   \begin{enumerate}%
   \makeatletter
   \input{\endnote@stream}%
   \end{enumerate}%
  \endgroup
 }%
}%
\appdef\class@enddocumenthook{%
 \printendnotes
}%
\chardef\REVfirstnote\z@
\AtBeginDocument{%
 \global\c@footnote\REVfirstnote
 \appdef\maketitle{\global\c@footnote\REVfirstnote}%
}%
\def\set@firstnote#1{%
 \@ifnum{\REVfirstnote=#1\relax}{}{%
  \REVTEX@warn@end{Endnote numbers changed: rerun LaTeX}%
 }%
 \immediate\write\@mainaux{%
   \global\chardef\string\REVfirstnote#1\relax
 }%
}%
\def\@endnotesinbib{%
 \set@firstnote\z@
}%
\def\@endnotesinbibliography{%
 \expandafter\set@firstnote\expandafter{\the\c@NAT@ctr}%
 \@ifx{\@endnoteout\@undefined}{}{%
  \immediate\closeout\@endnoteout
  \global\let\@endnoteout\@undefined
  \makeatletter
  \input{\endnote@stream}%
 }%
}%
\AtBeginDocument{%
 \let\REVTEX@footnote\footnote
 \footinbib@sw{%
  \let\footnote\endnote
  \let\printendnotes\relax
  \let\@endnotesinbib\@endnotesinbibliography
 }{}%
}%
\def\REVdoendnote#1#2{\bibitem{endnote#1}#2}
\let\@bibnotes\@empty
\AtBeginDocument{%
 \ps@article
 \pagenumbering{arabic}%
 \normalsize
 \tabbingsep \labelsep
 \leftmargin\leftmargini
 \labelwidth\leftmargin\advance\labelwidth-\labelsep
 \let\@listi\@listI
 \@listi
 \@ifxundefined\REVTEX@settypespec{}{\REVTEX@settypespec}%
 \twoside@sw{}{%
  \raggedbottom
 }%
 \multicol@sw{}{%
  \twocolumn@sw{%
   \twocolumn
   \sloppy
   \raggedbottom
  }{%
   \onecolumn
   \sloppy
   \raggedbottom
  }%
 }%
}%
\def\appendix{%
  \par
  \setcounter{section}\z@
  \setcounter{subsection}\z@
  \setcounter{subsubsection}\z@
  \def\thesection{\Alph{section}}%
  \def\thesubsection{\arabic{subsection}}%
  \def\thesubsubsection{\alph{subsubsection}}%
  \@addtoreset{equation}{section}%
  \def\theequation@prefix{\thesection}%
  \addtocontents{toc}{\protect\appendix}%
}%
\def\title@column@multicol#1{%
 \@twocolumntrue
 \let\@currenvir@save\@currenvir
 \def\@currenvir{multicols}%
 \multicols{2}[#1]\relax
}%
\def\close@column@multicol{%
 \endmulticols
 \let\@currenvir\@currenvir@save
 \hfuzz\maxdimen
}%
\newenvironment{widetext@multicol}{%
  \endgroup
  \def\@tempa{multicols}%
  \@ifx{\@currenvir\@tempa}{%
    \csname end\@currenvir\endcsname
  }{%
    \REVTEX@err{%
      widetext environment must not be inside any environment: \@currenvir
    }%
    \def\endwidetext{\begingroup\def\@currenvir{widetext}}%
  }%
}{%
  \multicols{2}\relax
  \begingroup
  \def\@currenvir{widetext}%
  \@endpetrue
}%
\def\title@column@default#1{#1}%
\def\close@column@default{}%
\def\galley@outdent{\rightmargin-\columnwidth\advance\rightmargin-\columnsep}%
\let\widetext@outdent\@empty
\newenvironment{widetext@galley}{%
  \list{}{%
    \listparindent \parindent
    \itemindent    \parindent
    \leftmargin    \z@
    \parsep        \z@\@plus\p@
    \widetext@outdent
    \relax
  }%
  \item\relax
}{
  \endlist
}%
\def\title@column@grid#1{%
 \@twocolumntrue
 \twocolumn[#1]%
}%
\def\close@column@grid{%
 \onecolumn
}%
\newenvironment{widetext@grid}{%
  \par\ignorespaces\twocolumn@grid
}{%
  \par\onecolumn@grid\global\@ignoretrue
  \@endpetrue
}%
\def\twocolumn@grid{\twocolumn\relax}%
\def\onecolumn@grid{\onecolumn}%
\AtEndOfClass{%
 \@ifxundefined\twocolumn@sw{}{%
  \twocolumn@sw{%
   \multicol@sw{%
     \let\REVTEX@tempa\@classoptionslist
     \let\@classoptionslist\@empty
     \RequirePackage{multicol}%
     \let\@classoptionslist\REVTEX@tempa
   }{}%
  }{}%
 }%
}%
\AtBeginDocument{%
 \twocolumn@sw{%
  \multicol@sw{%
   \let@environment{widetext}{widetext@multicol}%
   \let\title@column\title@column@multicol
   \let\close@column\close@column@multicol
  }{%
   \let@environment{widetext}{widetext@grid}%
   \let\title@column\title@column@grid
   \let\close@column\close@column@grid
  }%
 }{%
  \let@environment{widetext}{widetext@galley}%
  \preprintsty@sw{%
  }{%
   \galley@sw{%
    \let\widetext@outdent\galley@outdent
   }{%
   }%
  }%
 }%
}%
\AtBeginDocument{%
 \@ifxundefined\title@column{%
  \let\title@column\title@column@default
 }{}%
}%
\DeclareOldFontCommand{\rm}{\normalfont\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\normalfont\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\normalfont\ttfamily}{\mathtt}
\DeclareOldFontCommand{\bf}{\normalfont\bfseries}{\mathbf}
\DeclareOldFontCommand{\it}{\normalfont\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\normalfont\slshape}{\@nomath\sl}
\DeclareOldFontCommand{\sc}{\normalfont\scshape}{\@nomath\sc}
\DeclareRobustCommand*\cal{\@fontswitch\relax\mathcal}
\DeclareRobustCommand*\mit{\@fontswitch\relax\mathnormal}
\def\today{\ifcase\month\or
  January\or February\or March\or April\or May\or June\or
  July\or August\or September\or October\or November\or December\fi
  \space\number\day, \number\year}
\def\notesname{Notes}
\def\partname{Part}
\def\tocname{Contents}
\def\lofname{List of Figures}
\def\lotname{List of Tables}
\def\refname{References}
\def\indexname{Index}
\def\figurename{FIG.}
\def\figuresname{Figures}%
\def\tablename{TABLE}
\def\tablesname{Tables}%
\def\abstractname{Abstract}
\def\appendixesname{Appendixes}%
\def\appendixname{Appendix}%
\def\acknowledgmentsname{Acknowledgments}
\def\journalname{??}
\def\copyrightname{??}
\def\andname{and}
\def\@pacs@name{PACS numbers: }%
\def\ppname{pp}
\def\numbername{number}
\def\volumename{volume}
\def\Dated@name{dated: }%
\def\Received@name{received }%
\def\Revised@name{revised }%
\def\Accepted@name{accepted }%
\def\Published@name{published }%
\AtBeginDocument{%
  \providecommand\MakeTextUppercase{\MakeUppercase}%
}%
\def\address{\replace@command\address\affiliation}%
\def\altaddress{\replace@command\altaddress\altaffiliation}%
\newenvironment{references}{%
 \REVTEX@warn@end{The references environment is not supported; use thebibliography instead.}
 \gdef\references{\thebibliography{}}\references
}{%
 \endthebibliography
}%
\def\draft{%
 \REVTEX@warn@end{Command \string\draft\space is obsolete;^^JInvoke option draft instead.}%
 \@booleantrue\draft@sw
}%
\def\tighten{%
 \REVTEX@warn@end{Command \string\tighten\space is obsolete;^^JInvoke option tightenlines instead.}%
 \@booleantrue\tightenlines@sw
}%
\def\tableline{%
 \noalign{%
  \REVTEX@warn@end{Command \string\tableline\space is obsolete;^^JUse \string\colrule\space instead.}%
  \global\let\tableline\colrule
 }%
 \tableline
}%
\def\case{\replace@command\case\frac}%
\def\slantfrac{\replace@command\slantfrac\frac}%
\def\tablenote{\replace@command\tablenote\footnote}%
\def\tablenotemark{\replace@command\tablenotemark\footnotemark}%
\def\tablenotetext{\replace@command\tablenotetext\footnotetext}%
\DeclareRobustCommand\REV@text[1]{%
 \relax
 \ifmmode
  \mathchoice
   {\hbox{{\everymath{\displaystyle     }#1}}}%
   {\hbox{{\everymath{\textstyle        }#1}}}%
   {\hbox{{\everymath{\scriptstyle      }\let\f@size\sf@size\selectfont#1}}}%
   {\hbox{{\everymath{\scriptscriptstyle}\let\f@size\ssf@size\selectfont#1}}}%
  \glb@settings
 \else
  \mbox{#1}%
 \fi
}%
\DeclareRobustCommand\REV@bbox[1]{%
 \relax
 \ifmmode
  \mathchoice
   {\hbox{{\everymath{\displaystyle     }\boldmath$#1$}}}%
   {\hbox{{\everymath{\textstyle        }\boldmath$#1$}}}%
   {\hbox{{\everymath{\scriptstyle      }\boldmath$#1$}}}%
   {\hbox{{\everymath{\scriptscriptstyle}\boldmath$#1$}}}%
  \glb@settings
 \else
  \mbox{#1}%
 \fi
}%
\DeclareRobustCommand\REV@bm[1]{%
 \REVTEX@warn@end{To use \string\bm, please load the bm package!}%
 \global\let\bm\relax
}%
\def\FL{\obsolete@command\FL}%
\def\FR{\obsolete@command\FR}%
\def\narrowtext{\obsolete@command\narrowtext}%
\def\mediumtext{\obsolete@command\mediumtext}%
\newenvironment{quasitable}{%
 \let@environment{tabular}{longtable}%
}{%
}%
\AtBeginDocument{%
 \@ifxundefined\text{\let\text\REV@text}{}%
 \@ifxundefined\bm{\let\bm\REV@bm}{}%
 \providecommand\bibinfo[2]{#2}%
 \providecommand\eprint[1]{#1}%
 \providecommand\url[1]{#1}%
}%
\def\bbox#1{%
 \REVTEX@warn@end{\string\bbox\space is obsolete,^^Jload the bm package and use \string\bm\space instead.}%
 \global\let\bbox\relax
}%
\newenvironment{mathletters}{%
 \REVTEX@warn@end{Environment {mathletters} is obsolete;^^Jload the amsmath package and use {subequations}!}%
 \global\let\mathletters\@empty
}{%
}%
\def\eqnum#1{%
 \REVTEX@warn@end{\string\eqnum\space is obsolete, load the amsmath package and use \string\tag!}%
 \global\let\eqnum\@gobble
}%
\RequirePackage{revsymb}%
\AtEndOfClass{\InputIfFileExists{\jobname.rty}{}{}}%
\AtEndOfClass{%
 \@ifxundefined\secnumarabic@sw{\@booleanfalse\secnumarabic@sw}{}%
 \secnumarabic@sw{%
  \def\thesection       {\arabic{section}}%
  \def\thesubsection    {\thesection.\arabic{subsection}}%
  \def\thesubsubsection {\thesubsection.\arabic{subsubsection}}%
 }{}%
}%
\AtEndOfClass{%
 \AtBeginDocument{%
  \@ifxundefined\@bibstyle{}{%
   \expandafter\LATEX@bibliographystyle\expandafter{\@bibstyle}%
  }%
 }%
}%
\AtEndOfClass{%
 \AtBeginDocument{%
  \defaults@hook
 }%
}%
\@ifxundefined\defaults@hook{\let\defaults@hook\@empty}{}%
\endinput
%%
%% End of file `revtex4.cls'.
