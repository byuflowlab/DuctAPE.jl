\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{byuthesis}[2023/05/01 v1.0 BYU thesis]
% Written by Tim McLain starting in September 2020.
% Significant contributions by Andrew Ning (who suggested using the memoir class as the foundation for
% the template), Ed Alvarez, Judd Mehr, and Victoria Stephens. A tip of the hat to Max Opgenoord -- the
% format of his MIT PhD thesis inspired many of the elements of the fancy format class option.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------- set user options for byuthesis class -----------
% Options are declared by user when specifying the document class such as in the template.tex file
% For example, \documentclass[fancy,phd]{byuthesis}
%.          or \documentclass[simple,masters]{byuthesis}
% Format options: fancy or simple (simple is default)
\newif\if@simple   % defines logical variable "simple" and initializes it to false
\@simpletrue % sets simple = true as default
\newif\if@fancy    % defines logical variable "fancy" and initializes it to false
\DeclareOption{fancy}{\@fancytrue \@simplefalse}  % if the fancy option is passed, sets fancy = true

% Degree options: masters or phd (masters is default)
\newif\if@masters
\@masterstrue
\newif\if@phd
\DeclareOption{phd}{\@phdtrue \@mastersfalse}

% Pagination options: oneside or twoside (oneside is default)
% NOTE: Two-sided printing is not a requirement of the college or university
%	Some departments or students in our college like the option of printing
%	two-sided hard copies. The twoside option is intended to provide this
%	capability. However, this option may require some tweaking, such as adding
%	blank pages in places to ensure that abstract, acknowledgements, and first
%	pages of chapters start the recto (right-hand) page. This two-sided printing
%	capability is provided as a courtesy. The College of Engineering does not
%	provide user support for two-sided printing.
%
\newif\if@oneside
\@onesidetrue
\newif\if@twoside
\DeclareOption{twoside}{\@twosidetrue \@onesidefalse}
% -------------------------------------

% -------- send main options to class -------------
\ProcessOptions\relax
\if@twoside
	% two-sided printing
	\LoadClass[letterpaper,11pt,twoside,openright]{memoir}  % based off of memoir class
\else
	% one-sided printing (default)
	\LoadClass[letterpaper,11pt,oneside,openright]{memoir}  % based off of memoir class
\fi
\RequirePackage{memhfixc}  % fix for hyperref w/ memoir
% ----------------------------------------------

% ----------- page layout -------------
\setlrmarginsandblock{1.75in}{1.75in}{*}
\setulmarginsandblock{1in}{*}{1}
\checkandfixthelayout
% ---------------------------------------

% ----------- typesetting ---------------
\usepackage{microtype}		% this package makes micro adjustments to reduce need for hyphenation
% ---------------------------------------

% ---------- fonts ---------------
% Use T1 font encoding (8-bit) for loaded fonts
\RequirePackage[T1]{fontenc}

%%%% Serif fonts for body text
% newpxtext uses palatino for body text while newpxmath uses palatino form math mode for a consistent look
\RequirePackage{newpxtext,newpxmath}

%% This font is not part of the standard LaTeX package and requires installation on your computer
%%\usepackage[urw-garamond]{mathdesign}
%%\usepackage{garamondx}

%%%% Monospace fonts for \ttfamily (such as for displaying computer code or computer commands)
\RequirePackage[scaled=0.85]{beramono}

%% This font is not part of the standard LaTeX package and requires installation on your computer
% I like this font more than beramono, however (TWM)
%\RequirePackage{inconsolata}
%
%%%% Sans serif font selection for \sffamily (for chapter and section titles)
\RequirePackage[defaultsans,scale=0.95]{opensans} %% Alternatively
% use the option 'defaultsans' instead of 'default' to replace the
% sans serif font only.

%%%% Enable semi bold of opensans font
\DeclareRobustCommand{\sbshape}{%
  \not@math@alphabet\sbshape\relax
    \fontseries{sb}\selectfont
}
% Usage: \textsb{Words to be set with semibold font}
 \DeclareTextFontCommand\textsb{\sbshape}

% -------------------------------

% ---------- colors -------------
\RequirePackage{xcolor}
\definecolor{navy}{HTML}{002E5D}
%\definecolor{royal}{HTML}{003DA5}	% university brand guide definition of royal
\definecolor{royal}{HTML}{005CAB}	% royal that matches BYU Engineering logo
%\definecolor{royal}{HTML}{0062B8}	% royal definition from A. Ning
\definecolor{darkgray}{HTML}{141414}
%\definecolor{mediumgray}{HTML}{7C878E}	% slate gray from university brand guide
\definecolor{mediumgray}{HTML}{666666}	% medium gray definition from A. Ning
\definecolor{black}{HTML}{111111}
\definecolor{primary}{HTML}{005CAB}
\definecolor{lighttext}{HTML}{444444}
% -------------------------------

% -------- figures ------------
\RequirePackage{graphicx}
\RequirePackage{rotating}
\RequirePackage[export]{adjustbox} % inner option used to align full-width figures to inner margin
% ------------------------------

% -------- math ----------
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
% ------------------------------

% ----------- algorithms -------
\RequirePackage[chapter]{algorithm}  % algorithms
\RequirePackage{algpseudocode}  % algorithms
\renewcommand{\algorithmicrequire}{\textbf{Input:}}
\renewcommand{\algorithmicensure}{\textbf{return}}
\algrenewcommand\textproc{\texttt}  % don't all caps function name
\newcommand{\algoreturn}{\textbf{return} }
\algrenewcommand{\algorithmiccomment}[1]{\textcolor{primary}{\hfill\sffamily\scriptsize\itshape $\triangleright$ #1}}
\algnewcommand{\LineComment}[1]{\State \textcolor{primary}{\sffamily\scriptsize\itshape $\triangleright$ #1}}
\algnewcommand{\To}{\textbf{to\ }}   % For use in '\For{$i = 1$ \To $n$}' -- replaced $:\ $
\algnewcommand{\By}{\textbf{by\ }}   % For use in '\For{$i = 1$ \To $n$ \By $k$}' -- replaced $:\ $
% -----------------------------

% ------------- nomeclature -------------
\RequirePackage[intoc]{nomencl}
\makenomenclature
\renewcommand{\nomname}{Nomenclature}
\setlength{\nomitemsep}{0pt}  % relaxed spacing (nominal spacing is looser than this)
% ---------------------------------------

% ------------------------------

% -------- Sections -----------
\setsecnumdepth{subsection}  % number chapters, sections, and subsections
\maxsecnumdepth{subsection}

% set up byuchapter style
\newcommand{\hsp}{\hspace{20pt}}
\newcommand*{\colorsectitle}{}
\newcommand*{\colorsubsectitle}{}
\makechapterstyle{byuchapter}{%
	\setlength{\beforechapskip}{54pt}		% Sets vertical space from top of page to baseline of chapter title to 2.2 in.
	\setlength{\afterchapskip}{48pt}		% Sets vertical space from baseline of chapter title to baseline of following text to be 0.87 in.
	\renewcommand*{\printchaptername}{}		% We don't print the word "Chapter"
   	\renewcommand*{\chapternamenum}{}		% No space inserted after "Chapter" and before \chapnum
   	\renewcommand*{\chapnumfont}{\normalfont\LARGE\sffamily\itshape\colorchapnum}
	\renewcommand*{\chaptitlefont}{\normalfont\LARGE\sffamily\itshape\colorchaptitle}
    \renewcommand*{\afterchapternum}{}		% No vertical space after "Chapter" before \chaptertitle
%	\renewcommand*{\printchapternum}{\chapnumfont \thechapter\hsp}
	\renewcommand*{\printchapternum}{%
		\@hangfrom{\chapnumfont \thechapter\hsp}}%	% Aligns second line of chapter title with beginning
}													% of first line rather than chapter number
\renewcommand*{\colorchapnum}{\color{royal}}
\renewcommand*{\colorchaptitle}{\color{royal}}
\renewcommand*{\colorsectitle}{\color{royal}}

\chapterstyle{byuchapter}

% appendix formatting
\aliaspagestyle{part}{empty}
\renewcommand{\parttitlefont}{\normalfont\huge\sffamily\itshape\colorchapnum\raggedright}  % format the font
\RequirePackage{etoolbox}		% provides \appto command used below (append to)
% \appto\appendix{\appendixpage}

% section heading formatting
\setsecheadstyle{\normalfont\large\sffamily\colorsectitle\raggedright}%
\setbeforesecskip{-18pt plus -4pt minus -1pt}
\setaftersecskip{6pt plus 1pt}

% subsection heading formatting
\setsubsecheadstyle{\normalfont\large\sffamily\colorsectitle\raggedright}%
\setbeforesubsecskip{-12pt plus -3pt minus -1pt}
\setaftersubsecskip{3pt plus 1pt}

% subsection heading formatting
\setsubsubsecheadstyle{\normalfont\sffamily\colorsectitle\raggedright}%
\setbeforesubsubsecskip{-12pt plus -2pt minus -1pt}
\setaftersubsubsecskip{1pt plus 0pt}

% Table of contents formatting (and LoF and LoT)
\renewcommand{\contentsname}{Table of Contents}
\renewcommand{\listfigurename}{List of Figures}
\renewcommand{\listtablename}{List of Tables}

%% change title style for toc, lof, lot
% Bold face, size LARGE titles
%\renewcommand{\printtoctitle}[1]{\LARGE\bfseries #1}
%\renewcommand{\printloftitle}[1]{\LARGE\bfseries #1}
%\renewcommand{\printlottitle}[1]{\LARGE\bfseries #1}

% Color, sans serif titles
\renewcommand{\printtoctitle}[1]{\normalfont\LARGE\sffamily\colorchaptitle\itshape #1}
\renewcommand{\printloftitle}[1]{\normalfont\LARGE\sffamily\colorchaptitle\itshape #1}
\renewcommand{\printlottitle}[1]{\normalfont\LARGE\sffamily\colorchaptitle\itshape #1}

% TODO: add color to the numbers in ToC, etc.
% how to do this without overwriting \thechapter, etc. which messes up the assumption box coloring?

% remove the dots between titles and page numbers
\renewcommand{\cftdot}{}

% % or you can put the numbers right after
% % put page number right after in toc, lof, lot
% \newcommand{\cfthsp}{\hspace{8pt}}			% controls the spacing between the title and the page number in toc, lof, lot
% \renewcommand{\cftchapterleader}{\cfthsp}
% \renewcommand{\cftchapterafterpnum}{\cftparfillskip}
% \renewcommand\cftchapterformatpnum[1]{~{\cftchapterpagefont #1}}
% \renewcommand{\cftsectionleader}{\cfthsp}
% \renewcommand{\cftsectionafterpnum}{\cftparfillskip}
% \renewcommand\cftsectionformatpnum[1]{~{\cftsectionpagefont #1}}
% \renewcommand{\cftfigureleader}{\cfthsp}
% \renewcommand{\cftfigureafterpnum}{\cftparfillskip}
% \renewcommand\cftfigureformatpnum[1]{~{\cftsectionpagefont #1}}
% \renewcommand{\cfttableleader}{\cfthsp}
% \renewcommand{\cfttableafterpnum}{\cftparfillskip}
% \renewcommand\cfttableformatpnum[1]{~{\cftsectionpagefont #1}}

% unbold titles in toc
\renewcommand{\cftchapterfont}{}
%\renewcommand{\cftchapterpagefont}{\color{royal}}	% make page entries on ToC, LoF, LoT royal blue in color
%\renewcommand{\cftsectionpagefont}{\color{royal}}
%\renewcommand{\cftfigurepagefont}{\color{royal}}
%\renewcommand{\cfttablepagefont}{\color{royal}}
\renewcommand{\cftchapterpagefont}{\color{black}}   % make page entries on ToC, LoF, LoT black in color
\renewcommand{\cftsectionpagefont}{\color{black}}
\renewcommand{\cftfigurepagefont}{\color{black}}
\renewcommand{\cfttablepagefont}{\color{black}}
% ----------------------------------------

%----------- figure captions -------------
\captionnamefont{\normalfont\small}
\captiontitlefont{\normalfont\small}

% % You need a newsubfloat element to use subcaption
% \newsubfloat{figure}
% % Command to set subcaption styles
% \subcaptionfont{\normalfont\small}
% ----------------------------------------

%---------------- header and footer -----------------
%\setlength{\headwidth}{\textwidth}

% create page style for opening pages of chapter that puts
% page number in the lower center of the page
\makepagestyle{byuchaptertitle}
\makeoddfoot{byuchaptertitle}{}{\normalfont\small\slshape\thepage}{}
\makeevenfoot{byuchaptertitle}{}{\normalfont\small\slshape\thepage}{}

% create headers and footers
\makepagestyle{byu}
%\aliaspagestyle{chapter}{byuchaptertitle}
\aliaspagestyle{plain}{byuchaptertitle}
\makerunningwidth{byu}[\textwidth]{\textwidth}
\makeheadposition{byu}{flushright}{flushleft}{flushright}{flushleft}
\preto\tableofcontents{\pagestyle{byuchaptertitle}}	% sets page numbering format in frontmatter

\makepsmarks{byu}{%
\nouppercaseheads
%\createmark{chapter}{both}{nonumber}{}{}
\createmark{chapter}{left}{nonumber}{}{}
\createmark{section}{right}{shownumber}{}{. \space}
\createplainmark{toc}{both}{\contentsname}
\createplainmark{lof}{both}{\listfigurename}
\createplainmark{lot}{both}{\listtablename}
\createplainmark{bib}{both}{\bibname}
\createplainmark{index}{both}{\indexname}
\createplainmark{glossary}{both}{\glossaryname}}

% Set up headers and footers for one-sided (default) or two-sided pagination
\if@twoside
	% two-sided pagination
	% Place page number in top outer corner of all pages except chapter title pages
	% Place chapter title in top inner corner of all pages except chapter title pages
	\makeoddhead{byu}%
	    {\normalfont\small\slshape\leftmark}{}{\normalfont\small\slshape\thepage}
	\makeevenhead{byu}%
	    {\normalfont\small\slshape\thepage}{}{\normalfont\small\slshape\leftmark}
	\makeoddfoot{byu}%
	    {}{}{}
	\makeevenfoot{byu}%
	    {}{}{}
\else
	% one-sided pagination (default)
	% Place page number in top right corner of all pages except chapter title pages
	% Place chapter title in top left corner of all pages except chapter title pages
	\makeoddhead{byu}%
		{\normalfont\small\slshape\leftmark}{}{\normalfont\small\slshape\thepage}	% include chapter title
	%    {}{}{\normalfont\small\slshape\thepage}	% don't include chapter title
	\makeoddfoot{byu}%
	    {}{}{}
\fi

\pagestyle{byu}
\appto\mainmatter{\pagestyle{byu}}
% -----------------------------------------------


% package for making an index at the end. needs to be loaded before hyperref package though, so it's here instead of in customcommands.tex
\usepackage{imakeidx}


% -------------- hyperlinks -------------
%\RequirePackage[colorlinks,bookmarks,bookmarksnumbered,allcolors=royal,breaklinks]{hyperref}  % make hyperlinks royal
%\RequirePackage[colorlinks,bookmarks,bookmarksnumbered,allcolors=black,breaklinks]{hyperref}	   % make hyperlinks black
%\RequirePackage[colorlinks,bookmarks,bookmarksnumbered,allcolors=mediumgray,breaklinks]{hyperref}	   % make hyperlinks mediumgray
\RequirePackage[colorlinks,bookmarks,bookmarksnumbered,hidelinks,breaklinks]{hyperref}	   % make hyperlinks same color as surrounding text

% hack to suppress annoying hyperref warnings related to \listoffigures, \listoftables, and \printnomenclature
\pdfstringdefDisableCommands{%
  \let\leavevmode@ifvmode\@empty
  \let\kern\@empty
}
% ---------------------------------------

% ---------- cross referencing ------------
\RequirePackage[capitalise]{cleveref}
% ---------------------------------------

% ----------- bibliography --------------
\if@fancy
	% do nothing for fancy style here, see fancy style options below
\else
	\usepackage[backend=biber,style=ieee,doi=true]{biblatex}
	\urlstyle{same}
	% Change type of DOI: and DOI URL to match text
	\DeclareFieldFormat{doi}{%
		DOI\addcolon\space
		\ifhyperref
		{\lowercase{\href{https://doi.org/#1}{\nolinkurl{#1}}}}
		{\lowercase{\nolinkurl{#1}}}}
\fi
% ---------------------------------------

% ----------- user input --------------
\def\@customtitle{\relax}
\newcommand{\customtitle}[1]{\gdef\@customtitle{#1}}

\def\@bauthor{\relax}
\newcommand{\bauthor}[1]{\gdef\@bauthor{#1}}

\def\@department{\relax}
\newcommand{\department}[1]{\gdef\@department{#1}}

\def\@institution{\relax}
\def\@shortinstitution{\relax}
\newcommand{\institution}[2][]{\gdef\@shortinstitution{\ifstrempty{#1}{#2}{#1}} \gdef\@institution{#2}}

\institution{Brigham Young University}

\def\@degree{\relax}
\newcommand{\degree}[1]{\gdef\@degree{#1}}

\if@phd
	\degree{Doctor of Philosophy}
\else
	\degree{Master of Science}
\fi

\newcommand*{\@committeelist}{}
\newcommand*{\committeemember}[1]{\listadd{\@committeelist}{#1}}
\newcommand*{\committeechair}[1]{\committeemember{#1, Chair}}

\def\@copyrightyear{\number\the\year}
\newcommand{\copyrightyear}[1]{\gdef\@copyrightyear{#1}}

\def\@keywords{\relax}
\newcommand{\keywords}[1]{\gdef\@keywords{#1}}
% ---------------------------------------

% ----------- title page --------------
\RequirePackage{etoolbox}
\newcommand*{\@printcommitteemember}[1]{#1\\}
\newlength{\titlebottomadjust}  \setlength{\titlebottomadjust}{-4pt}	% defined to tweak the length of the title page so that
																		% last line of text is at 1 inch above page bottom
\newcommand{\titlepage}{%
	% for fancy style, set the page margins to be equal to center title on page
	\if@fancy
		\setlrmarginsandblock{1.75in}{1.75in}{*}
		\checkandfixthelayout
	\fi
	\pdfbookmark{BYU TITLE PAGE}{titlepage}
	\pagestyle{empty}
	\begin{Spacing}{2.0}
		\begin{center}
			\@title \\
			\vspace{\fill}
			\@author \\
			\vspace{\fill}
			\begin{SingleSpace}
				A \if@phd dissertation \else thesis \fi submitted to the faculty of \\
				\@institution \\
				in partial fulfillment of the requirements for the degree of \\
			\end{SingleSpace}
			\@degree \\
			\vspace{\fill}
			\begin{SingleSpace}
				\forlistloop{\@printcommitteemember}{\@committeelist}
			\end{SingleSpace}
			\vspace{\fill}
			\@department \\
			\@institution \\
			\vspace{\fill}
			Copyright \textcopyright{} \@copyrightyear{} \@author \\
			All Rights Reserved
			\par
			\enlargethispage{\titlebottomadjust}
			\clearpage
		\end{center}
	\end{Spacing}
	% for fancy style, reset the page layout back for margin notes
	\if@fancy
		\setlrmarginsandblock{1.0in}{2.75in}{*}
		\checkandfixthelayout
		\setmarginnotes{0.125in}{2.0in}{\onelineskip}
	\fi
}
% ----------------------------------------

% -------- custom title page -------------
% This title page follows immediately after the title page required by BYU Grad Studies.
% This title page can be customized as the student desires, but should include the title,
% student name, department, degree, and the BYU Engineering logo.
%

\newcommand{\customtitlepage}{%
	\pdfbookmark{TITLE PAGE}{customtitlepage}
	\pagestyle{empty}
	\begin{Spacing}{2.2}
		\vspace*{116pt}
		\noindent
		{\color{royal}\Huge\sffamily\itshape\@customtitle} \\
		\vspace{\fill}
	\end{Spacing}
	\begin{Spacing}{1.4}
		{\Large
		\noindent
		\@author \\
		\@department \\
		\@degree} \\
		\vspace{\fill}
		\begin{figure}[htbp]
			\includegraphics[width=3.5in]{college_logo}
		\end{figure}
	\end{Spacing}
}
% ----------------------------------------

% ------------ abstract ------------------
\renewenvironment{abstract}{%
	\clearpage
	\pdfbookmark[0]{ABSTRACT}{abstract}
	\setlength{\parskip}{\baselineskip}
	\vspace*{35pt}	% spacing set to match chapter header spacing
	\begin{flushleft}
		% {\normalfont\LARGE\bfseries Abstract}    						 % black serif boldface
		{\normalfont\LARGE\sffamily\colorchaptitle\itshape{Abstract}}    % royal sans serif italics
		\par
	\end{flushleft}
	%% The following items should be omitted if the optional custom title page is used
	% \@title \\[\baselineskip]
	% \@author \\
	% \@department, \@shortinstitution \\
	% \@degree
	\vspace*{24pt}  % spacing set to match chapter header spacing
}{%
	\par
	\vspace{\fill}
	\noindent
	Keywords: \@keywords
}
% ----------------------------------------

% ---------- dedication ------------
\newenvironment{dedication}{%
	\clearpage
	\pdfbookmark[0]{DEDICATION}{dedication}
	\setlength{\parskip}{\baselineskip}
	\vspace*{35pt}	% spacing set to match chapter header spacing
	\begin{flushleft}
		{\normalfont\LARGE\sffamily\colorchaptitle\itshape{Dedication}}    	% royal sans serif italics
	\end{flushleft}
	\vspace*{24pt}  % spacing set to match chapter header spacing
}
% ----------------------------------------

% ---------- acknowledgments ------------
\newenvironment{acknowledgments}{%
	\clearpage
	\pdfbookmark[0]{ACKNOWLEDGMENTS}{acknowledgments}
	\setlength{\parskip}{\baselineskip}
	\vspace*{35pt}	% spacing set to match chapter header spacing
	\begin{flushleft}
		% {\normalfont\LARGE\bfseries Acknowledgments}    						% black serif boldface
		{\normalfont\LARGE\sffamily\colorchaptitle\itshape{Acknowledgments}}    	% royal sans serif italics
	\end{flushleft}
	\vspace*{24pt}  % spacing set to match chapter header spacing
}
% ----------------------------------------

% ---------- epigraph ------------
\newenvironment{customepigraph}{%
	\clearpage
	\pdfbookmark[0]{EPIGRAPH}{epigraph}
	\setlength{\parskip}{\baselineskip}
	\vspace*{35pt}	% spacing set to match chapter header spacing
	% \begin{flushleft}
	% 	% {\normalfont\LARGE\bfseries Acknowledgments}    						% black serif boldface
	% 	{\normalfont\LARGE\sffamily\colorchaptitle\itshape{Epigraph}}    	% royal sans serif italics
	% \end{flushleft}
	\vspace*{24pt}  % spacing set to match chapter header spacing
}
% ----------------------------------------






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- ADD FANCY OPTIONS --- %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\if@fancy

% ----------- page layout -------------
% change the page layout to allow margin notes
\setlrmarginsandblock{1.0in}{2.75in}{*}
\setulmarginsandblock{1in}{*}{1}
\checkandfixthelayout
% ---------------------------------------

% ------------- marginnotes -------------
% for notes in the margin, define separation from the textblock, max note width and min vertical distance between notes
\setmarginnotes{0.125in}{2.0in}{\onelineskip}
% ---------------------------------------

% % -------------- chapter ----------------
% % Revise the chapter style to create a fancier heading
% \newcommand{\mhsp}{\hspace{10pt}}
% \makechapterstyle{byuchapter}{%
% 	\setlength{\beforechapskip}{54pt}
% 	\setlength{\afterchapskip}{48pt}
% 	\renewcommand*{\chapnamefont}{\normalfont\LARGE\itshape}
% 	\renewcommand*{\chapnumfont}{\normalfont\LARGE\sffamily\itshape\colorchapnum}
% 	\renewcommand*{\chaptitlefont}{\normalfont\LARGE\sffamily\itshape\colorchaptitle}
% 	\renewcommand*{\afterchapternum}{}
% 	\renewcommand*{\printchaptername}{}
% 	\setlength{\midchapskip}{15mm}
% 	\renewcommand*{\chapternamenum}{}
% 	\renewcommand*{\printchapternum}{%
% 		\sidebar{\raisebox{-76pt}[0pt][0pt]{\makebox[0pt][l]{%
% 					\resizebox{!}{\midchapskip}{\mhsp\chapnumfont\thechapter}}}}}
% 	\renewcommand*{\printchaptertitle}[1]{\chaptitlefont ##1}
% 	\renewcommand*{\afterchaptertitle}%
% 	{\vskip\afterchapskip}
% }
% \chapterstyle{byuchapter}

% -------------- chapter ----------------
% Revise the chapter style to create a fancier heading
\newcommand{\mhsp}{\hspace{10pt}}
\makechapterstyle{byuchapter}{%
	\setlength{\beforechapskip}{54pt}
	\setlength{\afterchapskip}{48pt}
	\renewcommand*{\chapnamefont}{\normalfont\LARGE\itshape}
	\renewcommand*{\chapnumfont}{\normalfont\LARGE\sffamily\bfseries\colorchapnum}
	\renewcommand*{\chaptitlefont}{\normalfont\LARGE\sffamily\itshape\colorchaptitle}
	\renewcommand*{\afterchapternum}{}
	\renewcommand*{\printchaptername}{}
	\setlength{\midchapskip}{15mm}
	\renewcommand*{\chapternamenum}{%
		\sidebar{\raisebox{-86.5pt}[0pt][0pt]{\rotatebox[origin=c]{90}{%
					\resizebox{\midchapskip}{12pt}{\sffamily\colorchaptitle Chapter}}}}}
	\renewcommand*{\printchapternum}{%
		\sidebar{\raisebox{-76pt}[0pt][0pt]{\makebox[0pt][l]{%
					\resizebox{!}{\midchapskip}{\mhsp\chapnumfont\thechapter}}}}}
	\renewcommand*{\printchaptertitle}[1]{\chaptitlefont ##1}
	\renewcommand*{\afterchaptertitle}%
	{\vskip\afterchapskip}
}
\chapterstyle{byuchapter}

% -------------- appendix ----------------
\makechapterstyle{byuappendix}{%
	\setlength{\beforechapskip}{54pt}
	\setlength{\afterchapskip}{48pt}
	\renewcommand*{\chapnamefont}{\normalfont\LARGE\itshape}
	\renewcommand*{\chapnumfont}{\normalfont\LARGE\sffamily\bfseries\colorchapnum}
	\renewcommand*{\chaptitlefont}{\normalfont\LARGE\sffamily\itshape\colorchaptitle}
	\renewcommand*{\afterchapternum}{}
	\renewcommand*{\printchaptername}{}
	\setlength{\midchapskip}{15mm}
	\renewcommand*{\chapternamenum}{%
		\sidebar{\raisebox{-86.5pt}[0pt][0pt]{\rotatebox[origin=c]{90}{%
					\resizebox{\midchapskip}{10pt}{\sffamily\colorchaptitle Appendix}}}}}
	\renewcommand*{\printchapternum}{%
		\sidebar{\raisebox{-76pt}[0pt][0pt]{\makebox[0pt][l]{%
					\resizebox{!}{\midchapskip}{\mhsp\chapnumfont\thechapter}}}}}
	\renewcommand*{\printchaptertitle}[1]{\chaptitlefont ##1}
	\renewcommand*{\afterchaptertitle}%
	{\vskip\afterchapskip}
}
\appto\appendix{\chapterstyle{byuappendix}}
% ---------------------------------------

% ------------- footnotes ----------------
\RequirePackage{marginfix}		% This package fixes the issue of margin notes overflowing the page: it makes them break into the next page instead

% we place footnotes in the side margin and typeset them using the parameters associated with \footnotes
\renewcommand*{\@makefntext}[1]{\makefootmark\raggedright #1}  % allow footnotes to be ragged right
%\renewcommand*{\thefootnote}{\fnsymbol{footnote}}	% use symbols to denote footnotes
\renewcommand*{\thefootnote}{\itshape\alph{footnote}}	% use letters to denote footnotes
\newcommand{\shsp}{\hspace{4pt}}
\footmarkstyle{\textsuperscript{\scriptsize\itshape #1}}
\renewcommand{\foottextfont}{\scriptsize}			% slightly smaller fontsize since we are placing footnotes in margin
\setlength{\footmarkwidth}{0em}
\setlength{\footmarksep}{0em}
\setlength{\footparindent}{1em}						% paragraph indent for a multi-paragraph footnote
%\counterwithout{footnote}{chapter}					% do not restart footnote numbering with each chapter
\footnotesinmargin									% put footnotes in the side margin
% ---------------------------------------

% ------------- reference margin notes for backreferencing ---------------
% used by byubib.bbx
\newcounter{rnote}
\newcommand\refnote[1]{\stepcounter{rnote}\marginpar{\MakeLowercase{#1}}}

% ---------------------------------------

% -------- figures ------------
\sidecapmargin{outer} % this puts the sidecaptions in the outer margin (supports both single-sided and two-sided printing)
\setsidecappos{b} % this aligns the side caption with the bottom of the figure
\setsidecaps{0.125in}{2.0in} % this controls the spacing from the edge of the text, and then the width of the box the caption sits in.
% this makes the name of the caption (e.g. Figure 1) bold, small, and blue (the number and word are controlled with separate commands)
\captionnamefont{\bfseries\footnotesize\color{royal}}
\captiontitlefont{\footnotesize}
% \subcaptionlabelfont{\bfseries\footnotesize\color{royal}}
% this redefinition left justifies the captions in both the right and left margins (needed for two-sided printing)
\renewcommand*{\sidecapstyle}{%
%%  \captionnamefont{\bfseries}
  \ifscapmargleft
    \captionstyle{\raggedright}%
  \else
    \captionstyle{\raggedright}%
\fi}

% ------------------------------

% ----------- bibliography --------------
% setup biblatex, specify citation and bibliography styles
\RequirePackage[
	hyperref=true,
	backref=true,		% adds page number where citation occurs to references
	citestyle=byubib,	% citation style (byubib.cbx)
	bibstyle=byubib,	% bibliography style (byubib.bbx)
	giveninits=true,
	doi=true,			% add doi information to bibliography
	url=true,			% adds hyperlinks to citations and bibliography
	maxcitenames=1,		% only one surname in citations
	mincitenames=1,
	dateabbrev=false,
	urldate=long,
	labeldateparts
]{biblatex}

% Change type of DOI: and DOI URL to match text
\urlstyle{same}
\DeclareFieldFormat{doi}{%
	DOI\addcolon\space
	\ifhyperref
	{\lowercase{\href{https://doi.org/#1}{\nolinkurl{#1}}}}
	{\lowercase{\nolinkurl{#1}}}}

\DeclareFieldFormat{labelnumberwidth}{\color{primary}#1}%
\setlength{\biblabelsep}{10pt}%
\renewbibmacro{in:}{}

% ------------  Citations ------------------------
% Some optional citation commands
\DeclareCiteCommand{\mciten} % only prints number and puts reference in margin
	{}
	{\printtext[bibhyperref]{\printfield{labelnumber}}%
	   \usebibmacro{sidecite}}
	{}
	{}
\DeclareCiteCommand{\citen} % only prints number, but not in margin (useful for figures)
	{}
	{{\printtext[bibhyperref]{\printfield{labelnumber}}}}
	{}
	{}
\DeclareCiteCommand{\nobibcite} % only prints number and puts reference in margin
	{}
	{\usebibmacro{sidecite}}
	{}
	{}
\DeclareCiteCommand{\sidecitenumber}
	{}
	{{\printtext[bibhyperref]{[\printfield{labelnumber}]}}}
	{}
	{}
\DeclareCiteCommand{\sidecitename}
	{}
	{\printnames{labelname}\unspace{\printtext[bibhyperref]{[\printfield{labelnumber}]}}}
	{}
	{}

%% make scite the default citation command for fancy format
\let\scite=\mcite		% use scite instead of replacing cite

% sidecite bibliography macro to create citations in the side margin
\DeclareFieldFormat{sidecite}{\marginpar{\linespread{1.0}\scriptsize#1}}
\newbibmacro{sidecite}{%
  \printtext[sidecite]{%
     \raggedright%
     \printtext[bibhyperref]{\textcolor{primary}{\printfield{labelnumber}}\:}
     \textcolor{mediumgray}{%
	    \printnames{labelname}%
	    \unspace, %
		\printfield[citetitle]{labeltitle}%
		\newunit
		\usebibmacro{shortjournal}%
		\unspace, %
		\printfield{year}.%
		\vspace{-0.75\baselineskip}%    % this command decreases the spacing between citations when they come in groups
		}%
     }%
  \setunit{}}
% ---------------------------------------

% ---------- for testing ---------------
\newcommand{\myshorttext}{{\color{mediumgray} Lorem ipsum dolor sit amet, consectetuer adipiscing elit. Etiam lobortis facilisis sem. Nullam nec mi et neque pharetra sollicitudin. Praesent imperdiet mi nec ante.}}

\newcommand{\myshorttexta}{{\color{mediumgray} Lorem ipsum dolor sit amet, consectetuer adipiscing elit. Etiam lobortis facilisis sem. Nullam nec mi et neque pharetra sollicitudin. Praesent imperdiet mi nec ante. Donec ullamcorper, felis non sodales commodo, lectus velit ultrices augue, a dignissim nibh lectus placerat pede. Vivamus nunc nunc, molestie ut, ultricies vel, semper in, velit.}}

% -------------------------------------

\fi

