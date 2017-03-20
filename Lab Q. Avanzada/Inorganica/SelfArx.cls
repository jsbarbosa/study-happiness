%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stylish Article
% Class File
% Version 2.1 (1/10/15)
%
% This template has been downloaded from:
% http://www.LaTeXTemplates.com
%
% Original author:
% Mathias Legrand (legrand.mathias@gmail.com)
% With extensive modifications by:
% Vel (vel@latextemplates.com)
%
% License:
% CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{SelfArx}[25/01/2012, v1.0]
\RequirePackage{ifthen}
\RequirePackage{calc}
\AtEndOfClass{\RequirePackage{microtype}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions*
\LoadClass{article}
\RequirePackage{ifpdf}      % Needed to pick between latex and pdflatex

%----------------------------------------------------------------------
%	FONTS
%----------------------------------------------------------------------

\RequirePackage{times}      % Loads the Times-Roman Fonts
\RequirePackage{mathptmx}   % Loads the Times-Roman Math Fonts

%----------------------------------------------------------------------
%	VARIOUS USEFUL PACKAGES
%----------------------------------------------------------------------
\RequirePackage[english,spanish]{babel} % Specify a different language here - english by default
\selectlanguage{spanish}
\RequirePackage{lipsum} % Required to insert dummy text. To be removed otherwise
\RequirePackage[utf8]{inputenc}
\RequirePackage{amsmath,amsfonts,amssymb}
\RequirePackage{graphicx,xcolor}
\RequirePackage{booktabs}

%\DefineBibliographyStrings{english}{%
%	andothers = {\em et\addabbrvspace al\adddot}
%}
%----------------------------------------------------------------------
%	MARGINS	
%----------------------------------------------------------------------

\RequirePackage[left=2cm,%
right=2cm,%
top=2.25cm,%
bottom=2.25cm,%
headheight=11pt,%
letterpaper]{geometry}%

%----------------------------------------------------------------------
%	FIGURES AND TABLES CAPTIONS
%----------------------------------------------------------------------

\RequirePackage[labelfont={bf,sf,small},%
labelsep=period,%
justification=raggedright]{caption}
\setlength{\abovecaptionskip}{0pt}
\setlength{\belowcaptionskip}{0pt}

%----------------------------------------------------------------------
%	PAGE HEADER
%----------------------------------------------------------------------

\RequirePackage{fancyhdr}  % Needed to define custom headers/footers
\RequirePackage{lastpage}  % Number of pages in the document
\pagestyle{fancy}          % Enables the custom headers/footers
% Headers
\lhead{}%
\chead{}%
\rhead{\small\sffamily\bfseries\@PaperTitle\  --- \thepage/\pageref{LastPage}}
% Footers
\lfoot{}%
\cfoot{}%
\rfoot{}%
\renewcommand{\headrulewidth}{0pt}% % No header rule
\renewcommand{\footrulewidth}{0pt}% % No footer rule

%----------------------------------------------------------------------
%	SECTION/SUBSECTION/PARAGRAPH SET-UP
%----------------------------------------------------------------------

\RequirePackage[explicit]{titlesec}
\titleformat{\section}
{\color{color1}\large\sffamily\bfseries}
{}
{0em}
{\colorbox{color2!0}{\parbox{\dimexpr\linewidth-2\fboxsep\relax}{$\blacksquare\quad$\scshape #1}}}
[]
\titleformat{name=\section,numberless}
{\color{color1}\large\sffamily\bfseries}
{}
{0em}
{\colorbox{color2!0}{\parbox{\dimexpr\linewidth-2\fboxsep\relax}{$\blacksquare\quad$\scshape #1}}}
[]  
\titleformat{\subsection}
{\color{color1}\sffamily\bfseries}
{\thesubsection}
{0.5em}
{#1}
[]
\titleformat{\subsubsection}
{\sffamily\small\bfseries}
{\thesubsubsection}
{0.5em}
{#1}
[]    
\titleformat{\paragraph}[runin]
{\sffamily\small\bfseries}
{}
{0em}
{#1} 
\titlespacing*{\section}{0pc}{3ex \@plus4pt \@minus3pt}{5pt}
\titlespacing*{\subsection}{0pc}{2.5ex \@plus3pt \@minus2pt}{0pt}
\titlespacing*{\subsubsection}{0pc}{2ex \@plus2.5pt \@minus1.5pt}{0pt}
\titlespacing*{\paragraph}{0pc}{1.5ex \@plus2pt \@minus1pt}{10pt}

%----------------------------------------------------------------------
%	TABLEOFCONTENTS SET-UP
%----------------------------------------------------------------------
\newlength{\tocsep} 
\setlength\tocsep{1.5pc} % Sets the indentation of the sections in the table of contents
\setcounter{tocdepth}{3} % Three levels in the table of contents section: sections, subsections and subsubsections

\usepackage{titletoc}
\contentsmargin{0cm}
\titlecontents{section}[\tocsep]
{\addvspace{4pt}\small\bfseries\sffamily}
{\contentslabel[\thecontentslabel]{\tocsep}}
{}
{\hfill\thecontentspage}
[]
\titlecontents{subsection}[\tocsep]
{\addvspace{2pt}\sffamily}
{\contentslabel[\thecontentslabel]{\tocsep}}
{}
{\ \titlerule*[.5pc]{.}\ \thecontentspage}
[]
\titlecontents*{subsubsection}[\tocsep]
{\footnotesize\sffamily}
{}
{}
{}
[\ \textbullet\ ]  

%----------------------------------------------------------------------
%	MULTIPLE AUTHOR SET
%----------------------------------------------------------------------  

\newcount\@authcnt
\newcount\@tmpcnt\@tmpcnt\z@

\def\@affiliation{%
	\ifnum\@tmpcnt<\@authcnt
	\global\advance\@tmpcnt1
	\raggedright \csname @auth\romannumeral\the\@tmpcnt\endcsname\hfill\\%
	\let\next\@affiliation \vskip1pt
	\else
	\let\next\relax
	\fi
	\next}

\newcommand{\affiliation}[1]{%
	\global\advance\@authcnt1
	\expandafter\gdef\csname @auth\romannumeral\the\@authcnt\endcsname
	{#1}}

%----------------------------------------------------------------------
%	LIST CONTROL
%----------------------------------------------------------------------  

\RequirePackage{enumitem}
%\setlist{nolistsep} % Uncomment to remove spacing between items in lists (enumerate, itemize)

%----------------------------------------------------------------------
%	ABSTRACT+AUTHOR FRAME
%----------------------------------------------------------------------  

\newcommand{\PaperTitle}[1]{\def\@PaperTitle{#1}}
\newcommand{\Archive}[1]{\def\@Archive{#1}}
\newcommand{\Authors}[1]{\def\@Authors{#1}}
\newcommand{\JournalInfo}[1]{\def\@JournalInfo{#1}}
\newcommand{\Abstract}[1]{\def\@Abstract{#1}}
%\newcommand{\Keywords}[1]{\def\@Keywords{#1}}

% ---------------------------------------------------------------------

\renewcommand{\@maketitle}{%
	\twocolumn[{%
		\thispagestyle{empty}%
		\vskip-36pt%
		{\raggedleft\small\sffamily\bfseries\@JournalInfo\\\@Archive\par}%
		\vskip20pt%
		{\raggedright\sffamily\bfseries\fontsize{20}{25}\selectfont \scshape\@PaperTitle\par}%
		\vskip10pt%
		{\raggedright\sffamily\fontsize{12}{16}\selectfont  \@Authors\par}%
		\vskip10pt%
		%\begingroup%
		{\raggedright\sffamily\footnotesize \@affiliation\par}
		%\endgroup%%
		\vskip20pt%
		{%
			\parbox{\textwidth-2\fboxsep-2\fboxrule}{\centering%
				\colorbox{color2!100}{
					\parbox{\textwidth-4\fboxsep-2\fboxrule}{%
						{\color{color1}\hrule}
						\vskip10pt
						\ifx\@Keywords\@empty%
						\sffamily\textbf{\abstractname}\\\@Abstract%
						\else%
						{\sffamily\color{color1}\scshape\textbf{Abstract: }}\@Abstract%\\[4pt]%
						%\textbf{\keywordname}\\\@Keywords%
						%\fi%
						{\color{color1}\hrule}
					}%
				}%
			}%
		}%
		\vskip25pt%
	}]%
}


%----------------------------------------------------------------------
%	REFERENCES
%----------------------------------------------------------------------  

% Remove brackets from numbering in List of References
%\renewcommand{\@biblabel}[1]{\bfseries\color{color1}\textsuperscript{[#1]}}
%\setlength{\bibitemsep}{0cm}
%\let\oldbibliography\thebibliography
%\renewcommand{\thebibliography}[1]{%
%\addcontentsline{toc}{section}{\refname}%
%\oldbibliography{#1}%
%\setlength\itemsep{0pt}}%

%----------------------------------------------------------------------
%	SCHEME
%---------------------------------------------------------------------- 
\RequirePackage{newfloat}
\DeclareFloatingEnvironment[fileext=los, listname=Lista de Esquemas, name=Esquema, placement=tbhp]{scheme}
\RequirePackage[version=4]{mhchem}

%----------------------------------------------------------------------
%	LABORATORIO AVANZADO
%----------------------------------------------------------------------
\setlength{\columnsep}{0.55cm} % Distance between the two columns of text
\setlength{\fboxrule}{0.75pt} % Width of the border around the abstract

\definecolor{color1}{RGB}{34,80,169} % Color of the article title and sections
\definecolor{color2}{RGB}{255,243,201} % Color of the boxes behind the abstract and headings

\RequirePackage{hyperref} % Required for hyperlinks
\hypersetup{hidelinks,colorlinks,breaklinks=true,urlcolor=color2,citecolor=color1,linkcolor=color1,bookmarksopen=false,pdftitle={Title},pdfauthor={Author}}

\Authors{Andr\'es Perez$^\dagger$, Juan Barbosa $^\ddagger$} % Authors
\affiliation{ \textit{Departamento de Qu\'imica}. $^\dagger$\textit{Departamento de Ing. Qu\'imica}. $^\ddagger$\textit{Departamento de F\'isica}. Universidad de los Andes, Bogot\'a, Colombia.} % Author affiliation