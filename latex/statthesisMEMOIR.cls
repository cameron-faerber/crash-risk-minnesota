%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% statthesisMEMOIR.cls
% Created 18 November 2009.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Provides etd, onside, sections, skiptitle, noabstract, noacknowl,
% finaldraftwithcolorlinks, and nochapdots options.
%
% Commands and functionality are patterned after Memoir, the
% foundation of this class.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Last modified 01/15/11

\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{statthesisMEMOIR}[2009/11/18]

\typeout{}
\typeout{--------------------------------------------------------------}
\typeout{         BYU THESIS CLASS FILE }
\typeout{               November 2009}
\typeout{         written by Andrew Olsen (aolsen256@gmail.com) }
\typeout{--------------------------------------------------------------}
\typeout{}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         Appropriately configure options for the class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This will keep the paper version front matter all on a new page

\newif\if@paper
\@papertrue

% This is for the electronic document.

\newif\if@etd
\DeclareOption{etd}{
\PassOptionsToClass{openany}{memoir}
 \@etdtrue }

% This is for documents that do not contain chapters.

\newif\if@sec
\DeclareOption{sections}{\@sectrue}

% This is if no acknowledgment should be included in the automatic front matter.

\newif\if@noacknowl
\@noacknowltrue
\DeclareOption{noacknowl}{\@noacknowlfalse}

% This is if no abstract should be included in the automatic front matter.

\newif\if@noabstract
\@noabstracttrue
\DeclareOption{noabstract}{\@noabstractfalse}

% This is only needed if we don't want a table of contents on electronic documents.
% Otherwise, we disable this feature.
%\newif\if@toc
%\@toctrue

% This option allows you to skip the automatic generation of front matter.

\newif\if@skiptitle
\@skiptitletrue
\DeclareOption{skiptitle}{\@skiptitlefalse}

% This simply allows the oneside option of memoir to be used.

\DeclareOption{oneside}{ \@paperfalse \PassOptionsToClass{oneside}{memoir}}

% This option allows for the colored link boxes to not appear. The links are still created.

\newif\if@colorlink
\DeclareOption{finaldraftwithcolorlinks}{\@colorlinktrue} 

% This removes the chapter-to-page dotted leader on the table of contents 

\newif\if@nochapdots
\@nochapdotstrue
\DeclareOption{nochapdots}{\@nochapdotsfalse}

\ProcessOptions \relax

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Load the Memoir class and set the margins as required
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\LoadClass[12pt]{memoir}


% Sets the margins as desired

\setstocksize{11in}{8.5in}
\setulmarginsandblock{1in}{1in}{*}
\setlrmarginsandblock{1in}{1in}{*}

% This handles the electronic document options:

\if@etd
	\@paperfalse
	% Sets the margins
	\setlrmarginsandblock{1in}{1in}{*} 
	%\@tocfalse
	% This just fixes one last thing where the beginning of mainmatter previously had to be on an odd page. 
	%Now it can be on either for this type of document
	\renewcommand*{\@smemmain}{%
	 	 \@mainmattertrue
	 		 \setcounter{secnumdepth}{\value{maxsecnumdepth}}
	  			\ifartopt
	    		\if@twoside
	     		% \cleardoublepage
	   			 \else
	    		% \clearpage
	   			 \fi
	  			\else
	    		%\cleardoublepage
	    		\@memmain@floats
	  		\fi}
\fi

\checkandfixthelayout

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Set up special commands associated with the appendices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These commands are part of rewriting the commands to produce the word Appendix on the bookmarks and table of contents.
% Typically only the letter appears.  Later on, similar code will be prepared for the sections option.

% Create a storage location for the chapter title.

\def\mychaptertitle#1{\gdef\@mychaptertitle{#1}}
\def\@mychaptertitle{}

% Initialize more counters that are needed for the process.

\newcounter{appendixchapter}
\setcounter{appendixchapter}{0}
\addtodef{\chapter}{    
	\setcounter{figure}{0}
    \setcounter{table}{0} 
    \setcounter{equation}{0}
    \stepcounter{appendixchapter}}{} 
\newcounter{appendixsection}
\setcounter{appendixsection}{0}
\addtodef{\section}{\stepcounter{appendixsection}}{} 

% Copy the chapter definition and define a new chapter command, \mynewchapter, that will be used with appendices.

\let \tempchapter = \chapter
\newcommand{\mynewchapter}[1]{\setcounter{appendixsection}{0}  % Sets the new counter to zero
	\mychaptertitle{#1} % Stores the input as a different chapter name.
	\tempchapter{\appendixname \space \@Alph\c@appendixchapter :\space #1} \refstepcounter{chapter}} % Outputs the input with the word Appendix A: in front.

% This does a similar thing for the section command.  

\let \tempsection = \section
\newcommand{\mynewsection}[1]{\tempsection{\@Alph \c@appendixchapter .\@arabic \c@appendixsection \space  \space #1} \refstepcounter{section}}


% This redefines the appendix command to make things work nicely.

\renewcommand{\appendix}{\par
  \setcounter{chapter}{0}			 	% Turns out this won't matter unless we reenable the numbering depth we later disable.
  \setcounter{appendixchapter}{0} 		% This is what we'll use as the counter.
  \setcounter{section}{0}% 				% Same story as with chapter.
  \gdef\@chapapp{\@Alph\c@appendixchapter}% This just clears something that is used later in chapter definitions.
  \anappendixtrue 						% This is used in the background of the appendix command.
  
  % These commands control the formatting and production of the extra page announcing appendices.  
  
  \renewcommand*{\appendixpagename}{\normalfont \normalsize \MakeUppercase{Appendices} }
  \appendixpage
  
 % This is the part where we disable the numbering of the chapters (at least the numbering that appears).

  \setcounter{secnumdepth}{-1} 
  
 % This allows the equations to be numbered correctly. 
  \renewcommand{\theequation}{ \thechapter. \@arabic\c@equation}

 % These commands redefine the chapter and section commands as outlined outside of the appendix redefinition.

\let\chapter=\mynewchapter
\let\section=\mynewsection
  
  
% This is the chapter style that will be used in typical documents for the appendices.  

\makechapterstyle{thatcherappendix}{%
  \chapterstyle{default}
  \renewcommand*{\chapterheadstart}{}
  \renewcommand*{\printchaptername}{%
    \centerline{\chapnumfont{\appendixname \space \@Alph\c@appendixchapter}}} % Here we hardcode the Appendix A title above the solid line.
  \renewcommand*{\chapternamenum}{}
  \renewcommand*{\chapnumfont}{\normalfont\scshape\MakeLowercase}
  \renewcommand*{\printchapternum}{}
  \renewcommand*{\afterchapternum}{%
    \par\centerline{\parbox{1.5in}{\hrulefill}}\par}
  \renewcommand*{\printchapternonum}{%
    \vphantom{\chapnumfont \@chapapp 1}\par
    \parbox{1.5in}{}\par}
  \renewcommand*{\chaptitlefont}{\normalfont}
  \renewcommand*{\printchaptertitle}{%
    \centering \chaptitlefont\MakeUppercase{\@mychaptertitle}}} % This is where we specify the chapter title that appears on the chapter page.  It is the input from the chapter command--without all of the things we added.
    
 \renewcommand{\thechapter}{\@Alph\c@appendixchapter}  %This is still needed for figure and table captions.
 \renewcommand*{\thesection}{\@Alph\c@appendixchapter.\@arabic\c@appendixsection} % This may be needed too.
 
 % Here we define how a chapter heading should appear. Essentially, we just remove the part that suppresses the chapter number when chapter numbering is disabled.
 \renewcommand{\@makechapterhead}[1]{%
  \chapterheadstart%  \vspace*{50\p@}%
  {%\parindent \z@ \raggedright \normalfont
   \parindent \z@ \memRTLraggedright \normalfont
   %\ifm@m@And
     \printchaptername \chapternamenum \printchapternum
     \afterchapternum % \par\nobreak \vskip 20\p@
   %\else
   %  \printchapternonum
   %\fi
   \interlinepenalty\@M
   \printchaptertitle{} % \Huge \bfseries #1
   \afterchaptertitle % \par\nobreak \vskip 40\p@
  }}	
  
  % This actually calls the chapter style that we have defined.
  
  \chapterstyle{thatcherappendix}
  
 % This concludes the redefinition of the appendix command.
  
%  \addtodef{\chapter}{}{  \refstepcounter{chapter} }
%  \addtodef{\section}{}{  \refstepcounter{section} }

  	}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   		Load Desired Packages for this Class
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\RequirePackage{amsmath,natbib,graphicx,float,xcolor}
% Note that natbib goes with the ASA Bib Style. Float allows us to force captions for graphics and tables correctly.
\RequirePackage[pdftex,bookmarks=TRUE,bookmarksnumbered=TRUE,bookmarksopenlevel=6,bookmarksopen=true,
pdfpagemode=UseOutlines,colorlinks=true,citecolor=black,linkcolor=black,menucolor=black]{hyperref}
% This package is used to create the bookmarks that perfectly match the table of contents.
\RequirePackage{memhfixc}
% This package is required for running hyperref with memoir.

%\definecolor{citationscolor}{rgb}{.3647059,.2784314,.5450980}
%\definecolor{toccolor}{rgb}{.2705882,0.545098,.4549020}

\definecolor{citationscolor}{rgb}{.545,0,0}
\definecolor{toccolor}{rgb}{0,0,.545}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   		Customize the chapter and heading styles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This changes the font on the Table of Contents to not be bold.

\renewcommand{\cftchapterfont}{ \normalsize }
\renewcommand{\cftchapterformatpnum}{ \normalsize }

% This redefines the thatcher chapter definition to use a normal sized font for the chapter name, etc.
% See the memoir class for the original definition.

\makechapterstyle{thatcher}{%
  \chapterstyle{default}
  \renewcommand*{\chapterheadstart}{}
  \renewcommand*{\printchaptername}{%
    \centerline{\chapnumfont{\@chapapp\ \thechapter}}}
  \renewcommand*{\chapternamenum}{}
  \renewcommand*{\chapnumfont}{\normalfont\scshape\MakeLowercase} % Controls the font of the heading.
  \renewcommand*{\printchapternum}{}
  \renewcommand*{\afterchapternum}{%
    \par\centerline{\parbox{1.5in}{\hrulefill}}\par}
  \renewcommand*{\printchapternonum}{%
    \vphantom{\chapnumfont \@chapapp 1}\par
    \parbox{1.5in}{}\par}
  \renewcommand*{\chaptitlefont}{\normalfont}
  \renewcommand*{\printchaptertitle}[1]{%
    \centering \chaptitlefont\MakeUppercase{##1}}} % This is where the font really comes into play.


% This redefines the Bringhurst chapter style to better suit our non-chapter documents
% by adding the section numbers (which become chapters) and tweaking a couple of font sizes.

\makechapterstyle{bringhurst}{%
  \chapterstyle{default}
  \renewcommand*{\chapterheadstart}{}
  \renewcommand*{\printchaptername}{}
  \renewcommand*{\chapternamenum}{}
  \renewcommand*{\printchapternum}{ \thechapter. }
  \renewcommand*{\afterchapternum}{}
  \renewcommand*{\printchaptertitle}[1]{%
%%%    \raggedright\Large\scshape\MakeLowercase{##1}}
    \memRTLraggedright\MakeUppercase{##1}}
  \renewcommand*{\afterchaptertitle}{%
 % \vskip\onelineskip 
 % \hrule
  \vskip\onelineskip}}

% This defines the Bringhurst head style with a few slight font modifications to correspond well with the previously defined Thatcher chapter style

\makeheadstyles{bringhurst}{%
\chapterstyle{bringhurst}
  \setbeforesecskip{-1\onelineskip
                    \@plus -0.5\onelineskip \@minus -.5\onelineskip}%
  \setaftersecskip{1\onelineskip \@plus 0.1\onelineskip}%
%%%  \setsecheadstyle{\normalfont\raggedright\scshape\MakeLowercase}%
  \setsecheadstyle{\normalfont\memRTLraggedright\scshape}%
  \setbeforesubsecskip{-1.0\onelineskip
                       \@plus -0.25\onelineskip \@minus -0.25\onelineskip}%
  \setaftersubsecskip{1.0\onelineskip \@plus 0.1\onelineskip}%
%%%  \setsubsecheadstyle{\sethangfrom{\noindent ####1}\normalfont\itshape\raggedright}%
  \setsubsecheadstyle{\sethangfrom{\noindent ####1}\normalfont\itshape\memRTLraggedright}%
  \setbeforesubsubsecskip{1.0\onelineskip
                          \@plus 0.5\onelineskip \@minus 0.2\onelineskip}%
  \setaftersubsubsecskip{-1em}%
  \setsubsubsecheadstyle{\normalfont\normalsize\scshape\MakeLowercase}%
  \setbeforeparaskip{1.0\onelineskip
                     \@plus 0.5\onelineskip \@minus 0.2\onelineskip}%
  \setafterparaskip{-1em}%
  \setparaheadstyle{\normalfont\normalsize\itshape\addperiod}%
  \setsubparaindent{\parindent}%
  \setbeforesubparaskip{1.0\onelineskip
                        \@plus 0.5\onelineskip \@minus 0.2\onelineskip}%
  \setaftersubparaskip{-1em}%
  \setsubparaheadstyle{\normalfont\normalsize\itshape\addperiod}}
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  				 Set some desired defaults
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These commands control orphans and widows.

\widowpenalty=10000
\clubpenalty=10000
\raggedbottom

% This sets the paragraph indent space.

\setlength{\parindent}{.5in}

% This tells LaTeX to use the bringhurst headings we developed.

\headstyles{bringhurst}

% This command loads the thatcher chapter style which is generally the chosen style.

\chapterstyle{thatcher}

% These commands are for formatting bookmarks and the table of contents.

\setcounter{tocdepth}{1}			 % This controls what appears in the Table of Contents and bookmarks list. 1 = Section
\setcounter{secnumdepth}{1}   		 % This controls the level that should be numbered. 1 = Section

\bibliographystyle{asabyu}			 % This utilizes the ASA's method of citation, modified by Andrew Olsen to reflect new changes.

\bibpunct{(}{)}{;}{a}{}{,}			 % This adds the correct in-text citations according to ASA standards.

\DoubleSpacing	% This is the default line spacing style for the document.

\pagestyle{plain} % This pagestyle gives us just the page number at the bottom of regular pages.

% These lines force the captions of figures to be displayed below the plots and for captions of tables to be above.

% I've commented them out because they are causing an issue in how I want two figures to show up on one page by themselves.

%\floatstyle{plain}					
%\restylefloat{figure}

%\floatstyle{plaintop}
%\restylefloat{table}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	This controls what happens when chapters are not desired.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\if@sec 
\chapterstyle{bringhurst}		% Uses the Bringhurst style meant for this purpose

% These lines turn sections into chapters, subsections into sections, etc.

\let\tempsectionchap=\chapter			
\let\tempsec=\section
\let\tempsubsec=\subsection
\let\tempsubsubsec=\subsubsection
\renewcommand{\section}{\tempsectionchap}
\renewcommand{\subsection}{\tempsec}
\renewcommand{\subsubsection}{\tempsubsec}



% This defines the chapter style for sections documents to be utilized in the appendix.

\makechapterstyle{bringhurstappendix}{%
  \chapterstyle{default}
  \renewcommand*{\chapterheadstart}{}
  \renewcommand*{\printchaptername}{}
  \renewcommand*{\chapternamenum}{}
  \renewcommand*{\printchapternum}{}
  \renewcommand*{\afterchapternum}{}
  \renewcommand*{\printchaptertitle}[1]{%
%%%    \raggedright \scshape\MakeLowercase{##1}}
    \memRTLraggedright \MakeUppercase{\appendixname \space \@Alph \c@appendixchapter: \space \@mychaptertitle}} % Hardcode the chapter name.
    % Here we have the chapter name including the word Appendix A: Then the Title
  \renewcommand*{\afterchaptertitle}{%
%  \vskip\onelineskip 
%  \hrule
  \vskip\onelineskip}}
 

% Set the appendix counters to zero.
\setcounter{appendixchapter}{0}
\setcounter{appendixsection}{0}


% Copy the chapter definition and define a new chapter command, \mynewchapter, that will be used with appendices.
\let \tempchapter = \tempsectionchap
\renewcommand{\mynewchapter}[1]{\setcounter{appendixsection}{0}  % Sets the new counter to zero
	\mychaptertitle{#1} % Stores the input as a different chapter name.
	\tempchapter{\appendixname \space \@Alph\c@appendixchapter :\space #1} \refstepcounter{chapter} } % Outputs the input with the word Appendix A: in front.

% This does a similar thing for the section command.  

\let \tempsection = \tempsec
\renewcommand{\mynewsection}[1]{\tempsection{\@Alph \c@appendixchapter .\@arabic \c@appendixsection \space  \space #1} \refstepcounter{section}}


% This redefines the appendix command to make things work nicely.

\renewcommand{\appendix}{\par
  \setcounter{chapter}{0} 			% Turns out this won't matter unless we reenable the numbering depth we later disable.
  \setcounter{appendixchapter}{0}	% This is what we'll use as the counter.
  \setcounter{section}{0}% 			% Same story as with chapter.
  \gdef\@chapapp{}% 				% This just clears something that is used later in chapter definitions.
  \anappendixtrue 					% This is used in the background of the appendix command.
  
  % These commands control the formatting and production of the extra page announcing appendices.
    
  \renewcommand*{\appendixpagename}{\normalfont \normalsize \MakeUppercase{Appendices} }
  \appendixpage
  
 % This is the part where we disable the numbering of the chapters (at least the numbering that appears).

  \setcounter{secnumdepth}{-1} 
  
  % This allows equations to be numbered correctly.
  \renewcommand{\theequation}{ \thechapter. \@arabic\c@equation}

 % These commands redefine the chapter and section commands as outlined outside of the appendix redefinition.
	
	\let\section=\mynewchapter
	\let\subsection=\mynewsection
  
 % This just redefines a couple of things and takes care of business.  
  
 \renewcommand{\thechapter}{\@Alph\c@appendixchapter}  							  % This is important for figures 
 \renewcommand*{\thesection}{\@Alph\c@appendixchapter.\@arabic\c@appendixsection} % May be needed
 
 % Here we define how a chapter heading should appear. Essentially, we just remove the part that suppresses the chapter number when chapter numbering is disabled.
 \renewcommand{\@makechapterhead}[1]{%
  \chapterheadstart%  \vspace*{50\p@}%
  {%\parindent \z@ \raggedright \normalfont
   \parindent \z@ \memRTLraggedright \normalfont
   %\ifm@m@And
     \printchaptername \chapternamenum \printchapternum
     \afterchapternum % \par\nobreak \vskip 20\p@
   %\else
   %  \printchapternonum
   %\fi
   \interlinepenalty\@M
   \printchaptertitle{} % \Huge \bfseries #1
   \afterchaptertitle % \par\nobreak \vskip 40\p@
  }}	
  
  % This actually calls the chapter style that we have defined.
  
  \chapterstyle{bringhurstappendix}
  
 % This concludes the redefinition of the appendix command.
  
  	}

 % This concludes the changes to be made if the sections command is present.


\fi



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 		Define parameters needed for the front matter.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\title#1{\gdef\@title{#1}}
\def\@title{\emph{Primary Title Here}}
\def\author#1{\gdef\@author{#1}}
\def\@author{Name of Student}
\def\typeofwork#1{\gdef\@typeofwork{#1}}
\def\@typeofwork{\emph{type of work}}
\def\university#1{\gdef\@university{#1}}
\def\@university{Brigham Young University}
\def\uinitials#1{\gdef\@uinitials{#1}}
\def\@uinitials{BYU}
\def\degree#1{\gdef\@degree{#1}}
\def\@degree{Master of Science}
\def\chair#1{\gdef\@chair{#1}}
\def\@chair{Scott Grimshaw}
\def\committee#1{\gdef\@committee{#1}}
\def\@committee{Committee Member 1}
\def\committeetwo#1{\gdef\@committeetwo{#1}}
\def\@committeetwo{Committee Member 2}
\def\committeethree#1{\gdef\@committeethree{#1}}
\def\@committeethree{\oneline}
\def\department#1{\gdef\@department{#1}}
\def\@department{Department of Statistics}
\def\graddate#1{\gdef\@graddate{#1}}
\def\@graddate{Month Year}
\def\copyrightyear#1{\gdef\@copyrightyear{#1}}
\def\@copyrightyear{2009}
\def\keywords#1{\gdef\@keywords{#1}}
\def\@keywords{Enter keywords here}

\def\acknowledgments#1{\gdef\@acknowledgments{#1}}
\def\@acknowledgments{I'd like to thank the academy. . . Oops! \oneline Either \oneline (1) add an \texttt{acknowledgments} statement to your preamble to replace this warning, \oneline (2) specify the \texttt{noacknowl} option in the document class statement to delete this page, or \oneline (3) specify the \texttt{skiptitle} option and build the frontmatter from scratch using the \texttt{frontmatter} and \texttt{maketitle} commands along with the \texttt{theabstract} and \texttt{theacknowledgments} environments if desired. End with the \texttt{mainmatter} command (all commands are preceded with a backslash).}

\def\abstract#1{\gdef\@abstract{#1}}
\def\@abstract{The shorter version of my masterpiece . . . is gone! \oneline Either \oneline (1) add an \texttt{abstract} statement to your preamble to replace this warning, \oneline (2) specify the \texttt{noabstract} option in the document class statement to delete this page, or \oneline (3) specify the \texttt{skiptitle} option and build the frontmatter from scratch using the \texttt{frontmatter} and \texttt{maketitle} commands along with the \texttt{theabstract} and \texttt{theacknowledgments} environments if desired. End with the \texttt{mainmatter} command (all commands are preceded with a backslash). Note that the third option is required if you want multiple paragraphs in your abstract.}

% This is just used in defining the title.
\newcommand{\oneline}{
\vskip\onelineskip
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	This defines the titlepage--accessible with \maketitle
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\let \oldtitle \maketitle

% This command can be used if you want the former functionality of maketitle

% Defines a new title page perscribed by Graduate Studies, 2009

\renewcommand*{\maketitle}{
\pagestyle{empty}
\begin{SingleSpace}
\begin{center}
\@title
\vfill
\@author
\vfill
A selected \@typeofwork \space submitted to the faculty of
 
\@university

in partial fulfillment of the requirements for the degree of 
\oneline
\@degree
\vfill
\@chair, Chair

\@committee

\@committeetwo

\@committeethree
\vfill
\@department
\oneline
\@university
\oneline
\@graddate
\vfill
\vspace*{-12pt}
Copyright \copyright  \space \@copyrightyear \space \@author
\oneline
All Rights Reserved
\end{center}
\end{SingleSpace}
\newpage
\pagestyle{plain}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 			This defines an abstract environment.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newenvironment{theabstract}{
\pdfbookmark[0]{Abstract}{Abstract}
\pagestyle{empty}
\SingleSpacing
\renewcommand{\oneline}{\\} %Adds new functionality as of January 2011
\begin{center}
ABSTRACT
\vskip \onelineskip
\@title
\vskip \onelineskip
\@author
\oneline
\@department, \@uinitials
\oneline
\@degree
\vskip \onelineskip
\end{center}
\renewcommand{\oneline}{\vskip \onelineskip}
}{


\vfill
\oneline
\noindent Keywords: \space \@keywords
\vspace{.25in}
\clearpage 
\pagestyle{plain}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 		This defines an acknowledgments environment.	
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newenvironment{theacknowledgments}{
\pdfbookmark[0]{Acknowledgments}{Acknowledgments}
\pagestyle{empty}
\begin{center}
ACKNOWLEDGMENTS
\end{center}
\oneline 

}{
\clearpage
\pagestyle{plain}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	This modifies the document environment to automatically include the 
%	title, abstract (unless noabstract is specified), acknowledgments
%	(unless noacknowl is specified), and table of contents (if etd is 
%	not specified). The option skiptitle results in the regularly defined
%	document environment.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\if@skiptitle

\let\tempbegdoc\document
\let\tempthree\enddocument

\renewenvironment{document}{
\tempbegdoc
\frontmatter
\pdfbookmark[0]{Title}{Title}
\maketitle
\if@noabstract
\if@paper
\pagestyle{empty}
%\cleartooddpage
\fi
\begin{theabstract}
\hspace{0.5in}
\@abstract
\end{theabstract}
\fi
\if@noacknowl
\if@paper
\pagestyle{empty}
%\cleartooddpage
\fi
\begin{theacknowledgments}
\@acknowledgments
\end{theacknowledgments}
\fi
\if@paper
\pagestyle{empty}
%\cleartooddpage
\fi
\pagestyle{plain}
%\if@toc
\tableofcontents
%\fi
\clearpage

\mainmatter

}
{
\tempthree
}

\fi

% This adds functionality that many people are accustomed to.

\newcommand{\singlespace}{\SingleSpace}
\newcommand{\doublespace}{\DoubleSpace}
\newcommand{\titlebreak}{\oneline}

% This allows the nochapdots option to remove the chapter leaders from the table of contents.

\if@nochapdots
\renewcommand{\cftchapterleader}{\cftdotfill{\cftdotsep}}
\fi

% This allows for the colored links to be enabled when desired.

\if@colorlink
\hypersetup{citecolor=citationscolor,linkcolor=toccolor,menucolor=black}
\fi

% This completes the definition of the statthesisMEMOIR class.

% Thank you ! \smile

\endinput
