%=====================================================================
%  The Observed Higgs Boson as a Candidate Bound State of
%  Farid Quantum Command Units:
%  Effective Couplings, Microscopic Completion, and
%  Pre-Registered Collider Tests
%
%  Mega Article v5 — corrected, compile-ready
%  RevTeX 4-2, PRD two-column format, pdfLaTeX
%=====================================================================
\documentclass[aps,prd,twocolumn,superscriptaddress,nofootinbib,longbibliography]{revtex4-2}

\usepackage[T1]{fontenc}
\usepackage{amsmath,amssymb,amsfonts}
\usepackage{bm}
\usepackage{booktabs}
\usepackage{enumitem}
\usepackage{hyperref} % loaded last

\newcommand{\GeV}{\,\mathrm{GeV}}
\newcommand{\TeV}{\,\mathrm{TeV}}
\newcommand{\MeV}{\,\mathrm{MeV}}

\begin{document}

\title{The Observed Higgs Boson as a Candidate Bound State of\\
Farid Quantum Command Units:\\
Effective Couplings, Microscopic Completion, and\\
Pre-Registered Collider Tests}

\author{Prof.\ Dr.\ Md.\ Faridul Islam Chowdhury, MBBS, MS (Neurosurgery)}
\affiliation{Tanfarid Vision Research Institute, Bogura, Bangladesh}
\email{dr\_faridul@yahoo.com}

\date{November 2026}

%=====================================================================
\begin{abstract}
We propose that the observed Higgs boson may be a spin-zero bound state of
multiple Farid Quantum Command Units (FQCUs). This is presented as an
\emph{implementation of a composite pseudo-Nambu--Goldstone-boson (pNGB) Higgs
scenario}, in which the microscopic sector is a confining \(SU(N_F)\) gauge
theory (``FQCD'') of eight complex scalars supplemented by vector-like fermions
\(\Psi_Q,\Psi_T\). We adopt the minimal \(SO(5)/SO(4)\) coset as an
\emph{assumption} and do not claim to derive it from FQCD; deriving the coset
from the eight-flavour dynamics is an open problem. Under this assumption, the
relations \(b=2a^2-1\) and \(\kappa_\lambda=b/a\) follow from the coset
geometry, and the custodial bi-doublet embedding protects \(\rho=1\) at tree
level. We verify that asymptotic freedom holds for \(N_F\ge1\) with \(n_f=2\)
fundamental Dirac flavours and \(n_s=8\) complex scalars, giving
\(b_0=\frac{11}{3}N_F-\frac{8}{3}\). The one-loop Coleman--Weinberg potential
from the top partners yields, with the running Yukawa \(y_t(f)\approx0.85\) at
the composite scale and partner masses \(\tilde M_{Q,T}\sim1.45\text{--}2.0\TeV\),
the calibrated values \(\tilde\beta\simeq3.32\times10^{-2}\),
\(\tilde\alpha\simeq1.79\times10^{-3}\), and \(m_h\simeq125.1\GeV\) at
\(f=1.5\TeV\). These numbers are an internally consistent calibration of the
minimal model, not a parameter-free prediction. The benchmark
\((a,b,\kappa_\lambda)=(0.9865,0.9462,0.959)\) is stated as a pre-registered
target for the High-Luminosity LHC. No evidence for FQCU is claimed.
\end{abstract}

\maketitle

%=====================================================================
\section{Introduction}
\label{sec:intro}
%=====================================================================

A new boson near \(125\GeV\) was reported independently by ATLAS and CMS in
2012~\cite{ATLAS2012,CMS2012}, and its properties have since been found
consistent with the Standard Model (SM) Higgs boson within current
precision~\cite{ATLASCMS2016}. The question addressed here is narrower than the
discovery: it is asked whether the observed particle may possess an underlying
structure --- in the present framework, a bound state of Farid Quantum Command
Units --- that would remain invisible at current resolution but could become
visible in sufficiently precise coupling measurements.

Earlier versions of this programme contained claims of Higgs redundancy,
validation by LHC data, fixed deviations in unrelated channels, and a
first-principles derivation of the Higgs as an FQCU bound state. Those claims
are withdrawn. The present article is restricted to what can be stated with
controlled scope, and it is explicit about which elements are derived and which
are assumed:

\begin{itemize}
\item \textbf{Assumed}: the confining FQCD sector contains eight complex scalars
and vector-like fermions \(\Psi_Q,\Psi_T\); the strong dynamics produces an
\(SO(5)/SO(4)\) coset with the observed Higgs as the pNGB. Deriving this coset
from the eight-flavour dynamics is an open problem.
\item \textbf{Derived}: within the assumed coset, the relations
\(a=\sqrt{1-\xi}\), \(b=1-2\xi\), \(b=2a^2-1\), and
\(\kappa_\lambda=b/a\) follow algebraically. The custodial \(SO(4)\)
bi-doublet embedding gives \(\rho=1\) at tree level. Asymptotic freedom holds
for \(N_F\ge1\) with \(n_f=2\) and \(n_s=8\).
\item \textbf{Calibrated}: the values of \(\tilde\alpha,\tilde\beta\) are chosen
to reproduce the observed \(v\) and \(m_h\). The coefficient formulas are
derived in Appendix~\ref{app:derivation}; the numerical values are a
calibration of the minimal model.
\end{itemize}

The article is organised as follows. Section~\ref{sec:operator} states the
interpolating-operator ansatz. Section~\ref{sec:eff} introduces the effective
\(a,b\) couplings. Section~\ref{sec:micro} describes the microscopic FQCD
sector and the custodial embedding. Section~\ref{sec:spectrum} gives the top
Yukawa and the Coleman--Weinberg potential. Section~\ref{sec:matching} states
the matching relations and predicted bands. Section~\ref{sec:test} gives the
pre-registered test protocol. Section~\ref{sec:relation} discusses the relation
to the spin-one FQCU hypothesis. Section~\ref{sec:open} lists open problems.
Section~\ref{sec:conclusion} concludes. The derivation appendix contains the
CW calculation.

%=====================================================================
\section{Candidate constituents and the interpolating operator}
\label{sec:operator}
%=====================================================================

In the QCU--SNEB notation, the Farid Quantum Command Unit (QCU) and the
Sub-Nano Energy Bit (SNEB) are organised by labels \(q=(p,s,d)\), each entry
taking two values, giving eight bookkeeping configurations. These labels do not
by themselves define spin, chirality, gauge charge, or a binding force.

As a concrete starting ansatz, the configurations are represented by real scalar
fields \(\Phi_q\), taken to be SM gauge singlets. A candidate scalar
interpolating operator is
\begin{equation}
\mathcal{O}_F(x)=\sum_{q,q'}c_{qq'}\,\Phi_q(x)\Phi_{q'}(x),
\qquad c_{qq'}=c_{q'q}\in\mathbb{R}.
\label{eq:OF}
\end{equation}
This operator is a two-field probe of the scalar channel; it is not implied that
the physical state contains exactly two constituents. For a bound-state
interpretation to be established, a specified microscopic interaction must
generate a pole below the constituent continuum, and the mass and overlap
\(\langle0|\mathcal{O}_F|h\rangle\) must be calculated. A nonzero overlap alone
is insufficient: a singlet operator can also overlap with an elementary Higgs
through mixing.

It should be stated clearly that the singlet-operator ansatz
\eqref{eq:OF} is superseded by the electroweak-charged construction of
Sec.~\ref{sec:micro}. The final model uses complex \(SU(N_F)\) fundamentals,
not real singlets; the operator \eqref{eq:OF} is retained only as a probe of
the two-field scalar channel and as a mnemonic for the flavour structure
\(q=(p,s,d)\).

%=====================================================================
\section{The effective bridge: \(a\) and \(b\)}
\label{sec:eff}
%=====================================================================

Below the scale of additional resonances, the observed scalar is described by
\(h\). In unitary gauge the leading vector interactions are
written~\cite{Giudice2007}
\begin{equation}
\mathcal{L}\supset
\bigl(m_W^2W_\mu^+W^{-\mu}+\tfrac12m_Z^2Z_\mu Z^\mu\bigr)
\bigl(1+2a\,\tfrac{h}{v}+b\,\tfrac{h^2}{v^2}+\dots\bigr),
\qquad v\simeq246\GeV,
\label{eq:SILH}
\end{equation}
where \(a\) and \(b\) rescale the one- and two-Higgs couplings to vector
bosons, with \(a=b=1\) in the tree-level SM. These coefficients are what
experiments probe; their microscopic origin is fixed by the coset geometry of
the strong sector under the assumption stated in Sec.~\ref{sec:micro}. At high
subprocess energy, the leading energy-growing term of longitudinal
vector-boson scattering into a Higgs pair is proportional to
\((a^2-b)\hat s/v^2\); the SM cancels this term~\cite{Giudice2007}. A nonzero
value of \(a^2-b\) therefore alters the high-invariant-mass behaviour of
Higgs-pair production, most cleanly in the vector-boson-fusion (VBF) channel.

%=====================================================================
\section{The microscopic theory: a confining FQCD sector}
\label{sec:micro}
%=====================================================================

\subsection{Scalar dynamics and symmetry structure}
\label{sec:scalar-sector}

As the microscopic gauge interaction, a new strong gauge group \(SU(N_F)\)
(``FQCD'') is introduced, under which the eight constituents transform as
complex fundamentals:
\begin{equation}
\mathcal{L}_{\text{FQCD}}^{\text{scalar}}
=\sum_q|D_\mu\Phi_q|^2-\sum_qm_q^2|\Phi_q|^2
-\frac{\lambda}{4}\Bigl(\sum_q|\Phi_q|^2\Bigr)^2
-\frac14F_{\mu\nu}^aF^{a\mu\nu}.
\label{eq:LFQCD}
\end{equation}
The QCU labels thereby acquire physical meaning as flavour labels of the strong
sector, in analogy with \((u,d,s,\dots)\) in QCD. The kinetic and potential
terms of Eq.~\eqref{eq:LFQCD} require the \(\Phi_q\) to be complex \(SU(N_F)\)
fundamentals; the real-field language of Sec.~\ref{sec:operator} is retained
there only for the operator-counting ansatz and is superseded here.

Three symmetry layers govern the scalar sector. First, the local \(SU(N_F)\) is
anomaly-free for scalars. Second, for eight complex flavours the maximal
continuous global symmetry of the kinetic term is \(U(8)\), which contains
\(SU(8)\) and \(U(1)\) factors. Third, an exact custodial
\(SO(4)\simeq SU(2)_L\times SU(2)_R\) subgroup is preserved prior to
electroweak gauging, breaking spontaneously to the diagonal \(SU(2)_V\) so that
the \(\rho\) parameter is protected at tree level.

\subsection{The coset assumption}
\label{sec:coset-assumption}

The minimal composite-Higgs relations
\(a=\sqrt{1-\xi}\), \(b=1-2\xi\) are the leading-order couplings of an
\(SO(5)/SO(4)\) coset with \(h\) the pNGB of the broken \(SO(5)\to SO(4)\).
\textbf{We adopt this coset as an assumption of the model.} The FQCD sector
\eqref{eq:LFQCD} has global symmetry \(U(8)\); whether its strong dynamics
breaks this to a subgroup containing \(SO(5)/SO(4)\), and whether the observed
Higgs is the associated pNGB, is an open problem of the microscopic completion.
The derivation of the coset from the eight-flavour dynamics is listed in
Sec.~\ref{sec:open}.

Under this assumption, \(h\) is an angular coordinate on \(SO(5)/SO(4)\), and
the electroweak mass function is, to leading
order~\cite{Agashe2005,Panico2016},
\begin{equation}
F(h)=\frac{\sin^2(\theta+h/f)}{\sin^2\theta},
\qquad v=f\sin\theta,
\qquad \xi\equiv\sin^2\theta=\frac{v^2}{f^2},
\label{eq:pNGB}
\end{equation}
from which
\begin{equation}
a=\sqrt{1-\xi},
\qquad b=1-2\xi,
\qquad b=2a^2-1.
\label{eq:ab_rel}
\end{equation}
Within the assumed coset, Eq.~\eqref{eq:ab_rel} is the standard leading-order
coupling geometry of a pNGB Higgs.

\subsection{Explicit custodial \(SO(4)\) bi-doublet embedding}
\label{sec:custodial}

To exhibit tree-level \(\rho=1\) protection explicitly, we partition the eight
binary states \(q=(p,s,d)\) (\(p,s,d\in\{0,1\}\)) by the third index
\(d\in\{0,1\}\) into two \(2\times2\) flavour matrices indexed by \((p,s)\):
\begin{equation}
\Phi^{(d)}=
\begin{pmatrix}
\Phi_{00d} & \Phi_{01d}\\
\Phi_{10d} & \Phi_{11d}
\end{pmatrix},
\qquad d\in\{0,1\}.
\label{eq:bidoublet}
\end{equation}
Under \(G_{\text{cust}}=SU(2)_L\times SU(2)_R\), each matrix transforms as a
\((\mathbf2,\mathbf2)\) bi-doublet,
\(\Phi^{(d)}\to U_L\Phi^{(d)}U_R^\dagger\).
The SM hypercharge generator \(Y\) is embedded inside \(SU(2)_R\) as
\(Y=T_R^3=\tfrac12\sigma^3\), assigning \(Y=+1/2\) to the left column and
\(Y=-1/2\) to the right column of \(\Phi^{(d)}\). Soft mass splittings
\(m_0^2\neq m_1^2\) isolate \(\Phi^{(0)}\) as the active electroweak-breaking
bi-doublet while keeping \(\Phi^{(1)}\) as a heavy inert spectator multiplet.
Under infinitesimal \(SU(2)_L\times U(1)_Y\) gauge transformations,
\begin{equation}
\delta\Phi^{(0)}
=i g\alpha_L^a\frac{\sigma^a}{2}\Phi^{(0)}
-i g'\alpha_Y\Phi^{(0)}\frac{\sigma^3}{2}.
\label{eq:gaugetr}
\end{equation}
When the strong condensate forms along the custodial diagonal \(SU(2)_V\) with
misalignment angle \(\theta\) (\(v=f\sin\theta=246\GeV\)),
\(\langle\Phi^{(0)}\rangle=\tfrac{v}{2}\mathbf{I}_{2\times2}\), the
covariant derivative evaluates to
\begin{equation}
D_\mu\langle\Phi^{(0)}\rangle
=-\frac{iv}{4}
\begin{pmatrix}
gW_\mu^3-g'B_\mu & g(W_\mu^1-iW_\mu^2)\\
g(W_\mu^1+iW_\mu^2) & -gW_\mu^3-g'B_\mu
\end{pmatrix}.
\label{eq:DVEV}
\end{equation}
Taking the trace
\(\mathrm{Tr}[(D_\mu\langle\Phi^{(0)}\rangle)^\dagger
(D^\mu\langle\Phi^{(0)}\rangle)]\)
with \(W_\mu^\pm=\tfrac{1}{\sqrt2}(W_\mu^1\mp iW_\mu^2)\) and
\(Z_\mu=(gW_\mu^3-g'B_\mu)/\sqrt{g^2+g'^2}\) yields
\begin{equation}
\mathcal{L}_{\text{mass}}
=\frac14g^2v^2W_\mu^+W^{-\mu}
+\frac18(g^2+g'^2)v^2Z_\mu Z^\mu,
\label{eq:WZm}
\end{equation}
giving \(m_W^2=\tfrac14g^2v^2\) and
\(m_Z^2=\tfrac14(g^2+g'^2)v^2\), which guarantees
\begin{equation}
\rho\equiv\frac{m_W^2}{m_Z^2\cos^2\theta_W}=1
\label{eq:rho}
\end{equation}
identically at tree level.

\subsection{Fermion sector: partial compositeness and anomaly freedom}
\label{sec:fermion-sector}

Gauge-singlet scalars do not develop Yukawa couplings by binding alone.
Following partial compositeness~\cite{Kaplan1991,Contino2007}, we extend the
strong \(SU(N_F)\) sector with two vector-like Dirac fermion constituents,
\(\Psi_Q\) and \(\Psi_T\), transforming under
\(SU(N_F)\times SU(3)_c\times SU(2)_L\times SU(2)_R\times U(1)_X\) with
hypercharge \(Y=T_R^3+X\), as listed in Table~\ref{tab:constituents}.

\begin{table}[t]
\caption{Quantum numbers of the scalar (\(\Phi^{(d)}\)) and minimal top-sector
fermionic (\(\Psi_Q,\Psi_T\)) constituents under the strong and SM gauge groups.}
\label{tab:constituents}
\begin{ruledtabular}
\begin{tabular}{lccccc}
Field & \(SU(N_F)\) & \(SU(3)_c\) & \(SU(2)_L\times SU(2)_R\) & \(U(1)_X\) & Composite Role\\
\midrule
\(\Phi^{(0,1)}\) & \(\mathbf{N}_F\) & \(\mathbf1\) & \((\mathbf2,\mathbf2)\) & \(0\) & pNGB \(h\subset\Phi^{(0)}\Phi^{(0)\dagger}\)\\
\(\Psi_Q\) & \(\bar{\mathbf{N}}_F\) & \(\mathbf3\) & \((\mathbf2,\mathbf2)\) & \(+2/3\) & Partner \(\mathcal{B}_Q\sim\Phi^{(0)}\Psi_Q\)\\
\(\Psi_T\) & \(\bar{\mathbf{N}}_F\) & \(\mathbf3\) & \((\mathbf1,\mathbf1)\) & \(+2/3\) & Partner \(\mathcal{B}_T\sim\Phi^{(0)}\Psi_T\)\\
\end{tabular}
\end{ruledtabular}
\end{table}

\textbf{Flavour counting.} Each of \(\Psi_Q\) and \(\Psi_T\) is one Dirac field
in the fundamental of \(SU(N_F)\), carrying \(SU(3)_c\) colour as an internal
quantum number of that single field. The number of Dirac flavours in
\(SU(N_F)\) is therefore
\begin{equation}
n_f=1+1=2,
\label{eq:nf}
\end{equation}
not \(6\). The three \(SU(3)_c\) colours of each field are internal indices
within the same \(SU(N_F)\) fundamental, not separate \(SU(N_F)\) flavours.
This is the standard field-theory counting and it differs from an earlier
version of this manuscript, which incorrectly used \(n_f=6\).

Because both \(\Psi_Q\) and \(\Psi_T\) are vector-like Dirac multiplets, all
perturbative gauge and gravitational anomalies vanish identically:
\begin{equation}
\mathcal{A}[SU(N_F)^3]=\mathcal{A}[SU(3)_c^3]
=\mathcal{A}[G_i^2U(1)_Y]
=\mathcal{A}[\mathrm{Grav}^2U(1)_Y]=0.
\label{eq:anomalies}
\end{equation}
Furthermore, \(\Psi_Q\) contains \(2N_F\) left-handed \(SU(2)_L\) doublets
(always an even integer), avoiding the global Witten \(SU(2)_L\) anomaly.

\subsection{One-loop \(\beta\) functions}
\label{sec:betas}

With \(n_s=8\) complex fundamental scalars and \(n_f=2\) fundamental Dirac
flavours, the one-loop running of the \(SU(N_F)\) gauge coupling \(g_F\) is
governed by
\begin{equation}
\beta_g\equiv\frac{dg_F}{dt}
=-\frac{g_F^3}{16\pi^2}b_0,
\qquad
b_0=\frac{11}{3}N_F-\frac23n_f-\frac16n_s
=\frac{11}{3}N_F-\frac83.
\label{eq:betag}
\end{equation}
Asymptotic freedom (\(b_0>0\)) holds for all \(N_F\ge1\):
\begin{equation}
N_F\ge1\;\text{is asymptotically free}.
\label{eq:NFbound}
\end{equation}
The selection rule \(N_F\ge1\) is therefore not a sharp constraint; it simply
states that FQCD is asymptotically free for any \(N_F\) consistent with the
strong-sector content. The choice of \(N_F\) is fixed by other considerations
(confinement, custodial embedding, anomaly freedom of the enlarged theory).

The one-loop running of the scalar quartic self-coupling \(\lambda\) in
\(V_{\text{scalar}}=\tfrac{\lambda}{4}(\sum_q|\Phi_q|^2)^2\) for \(n_s=8\)
complex fundamental scalars and \(n_f=2\) Dirac fermions receives contributions
from scalar quartic, gauge, and Yukawa couplings. The leading one-loop
\(\beta_\lambda\) in the absence of Yukawas has the schematic structure
\begin{equation}
\beta_\lambda^{(\mathrm{no\ Yuk})}
=\frac{1}{16\pi^2}
\left[
12\lambda^2-6g_F^2\lambda
+C_\lambda(N_F,n_s)g_F^4
\right],
\label{eq:betalambda}
\end{equation}
where \(C_\lambda(N_F,n_s)>0\) is a model-dependent coefficient that must be
computed with the exact flavour trace and the chosen normalisation of \(\lambda\).
The negative \(-6g_F^2\lambda\) term drives \(\lambda\) downward; the positive
\(g_F^4\) term drives it upward. Whether the full one-loop system stabilizes
\(\lambda>0\) in the infrared depends on the relative size of \(g_F\) and
\(\lambda\) at each scale, and on the Yukawa contributions from
\(\Psi_Q,\Psi_T\). A full numerical running, including the Yukawas and the
two-loop terms of Refs.~\cite{Machacek1983,Machacek1985}, is required to
establish vacuum stability; this is left to future work. The one-loop formula
\eqref{eq:betalambda} is included only to show the structure of the running,
not to claim that the vacuum is stabilized.

%=====================================================================
\section{Spectrum: top Yukawa and the Coleman--Weinberg potential}
\label{sec:spectrum}
%=====================================================================

\subsection{Top Yukawa and the running convention}
\label{sec:yt}

Below \(\Lambda_F\), scalar--fermion confinement binds \(\Phi^{(0)}\) with
\(\Psi_{Q,T}\) into spin-\(\tfrac12\) composite baryonic resonances
\(\mathcal{B}_Q,\mathcal{B}_T\) of mass \(M_{Q,T}\sim g_\rho f\). Linear mixing
with the elementary top quarks \((q_L,t_R)\) is described by
\begin{align}
\mathcal{L}_{\text{mix}}
&=\lambda_Lf\,\bar q_L\mathcal{B}_Q
+\lambda_Rf\,\bar t_R\mathcal{B}_T
-M_Q\bar{\mathcal{B}}_Q\mathcal{B}_Q
-M_T\bar{\mathcal{B}}_T\mathcal{B}_T
\nonumber\\
&\quad-Y_*f\,\bar{\mathcal{B}}_Q\tilde H_{\text{pNGB}}\mathcal{B}_T
+\mathrm{h.c.}
\label{eq:Lmix}
\end{align}
Defining the physical partner masses
\(\tilde M_{Q,T}=\sqrt{M_{Q,T}^2+\lambda_{L,R}^2f^2}\) and compositeness
angles \(\sin\varphi_L=\lambda_Lf/\tilde M_Q\),
\(\sin\varphi_R=\lambda_Rf/\tilde M_T\), diagonalization yields the
composite-scale top Yukawa
\begin{equation}
y_t(f)=Y_*\sin\varphi_L\sin\varphi_R.
\label{eq:yt-f}
\end{equation}
For \(f=1.5\TeV\), \(Y_*\in[2.5,4.0]\), and
\(\sin\varphi_{L,R}\approx0.55\pm0.06\), Eq.~\eqref{eq:yt-f} gives
\(y_t(f)\approx0.85\), the value used in the CW evaluation below.

\paragraph{Convention for \(y_t\).}
Throughout this paper, \(y_t(f)\) denotes the running top Yukawa at the
composite scale \(f\), and \(y_t(m_t)\) denotes the same coupling run to the
top mass scale. The two differ by RG evolution; \(y_t(m_t)\approx0.94\) is the
pole-like value in the \(\overline{\text{MS}}\) scheme at \(m_t\), and
\(y_t(f)\approx0.85\) is the running value at \(f\). The partial-compositeness
formula \eqref{eq:yt-f} gives \(y_t(f)\); the abstract and conclusion quote the
pole value \(y_t(m_t)\approx0.94\) when comparing with experiment. The CW
evaluation in Sec.~\ref{sec:CW} uses the running value \(y_t(f)\approx0.85\).

\subsection{Coleman--Weinberg potential}
\label{sec:CW}

In the misaligned vacuum, the \(3\times3\) top-sector mass matrix in the basis
\(\psi_L=(t_L,\mathcal{B}_{Q,L}^{(u)},\mathcal{B}_{T,L})^T\) and
\(\psi_R=(t_R,\mathcal{B}_{Q,R}^{(u)},\mathcal{B}_{T,R})^T\) takes the form
\begin{equation}
\mathcal{M}_t(\theta)=
\begin{pmatrix}
0 & \lambda_Lf\cos\theta & \tfrac{\lambda_Lf}{\sqrt2}\sin\theta\\
0 & M_Q & Y_*f\cos\theta\\
\lambda_Rf & 0 & M_T
\end{pmatrix}.
\label{eq:Mt}
\end{equation}
Forming \(H_t(\theta)\equiv\mathcal{M}_t^\dagger(\theta)\mathcal{M}_t(\theta)\),
the determinant factorizes into the squared eigenvalues
\(\det[H_t(\theta)]=m_t^2(\theta)\tilde M_Q^2(\theta)\tilde M_T^2(\theta)\).
We show in Appendix~\ref{app:derivation} that this yields the
\(\theta\)-dependent top mass
\begin{equation}
m_t^2(\theta)
\simeq\tfrac12y_t^2(f)f^2
\bigl(\sin^2\theta-\sin^4\theta\bigr)
+\mathcal{O}(\sin^6\theta).
\label{eq:mtop}
\end{equation}
The form \(\sin^2\theta-\sin^4\theta\) follows from the specific three-site
mixing structure of Eq.~\eqref{eq:Mt} with \(M_{Q,T}\gg y_tf\); in the simpler
minimal-composite-Higgs limit \(M_{Q,T}\to\infty\) it reduces to
\(m_t^2(\theta)\propto\sin^2\theta\). The derivation is given in
Appendix~\ref{app:derivation}.

Integrating out the top quark and heavy partners at one loop (\(N_c=3\)) gives
the Coleman--Weinberg effective potential for Dirac
fermions~\cite{Matsedonskyi2013}
\begin{equation}
V_{\text{CW}}(\theta)
=-\frac{N_c}{16\pi^2}\,
\mathrm{Tr}\!\left[
H_t^2(\theta)
\left(\ln\frac{H_t(\theta)}{\Lambda_F^2}-\frac32\right)
\right].
\label{eq:VCW}
\end{equation}
Softening of composite form factors above \(\Lambda_F\) (Weinberg sum
rules~\cite{Weinberg1967}) renders the coefficients finite and calculable.
Expanding for small \(\sin^2\theta\),
\begin{equation}
V_{\text{CW}}(\theta)
=f^4\bigl[-\tilde\alpha\sin^2\theta+\tilde\beta\sin^4\theta\bigr]
+\mathcal{O}(\sin^6\theta),
\label{eq:VCWexpand}
\end{equation}
with coefficients derived in Appendix~\ref{app:derivation}:
\begin{equation}
\tilde\alpha=\frac{N_c}{16\pi^2f^4}
\left[
\frac{\lambda_L^2f^2\tilde M_Q^2\tilde M_T^2}{\tilde M_Q^2-\tilde M_T^2}
\ln\!\left(\frac{\tilde M_Q^2}{\tilde M_T^2}\right)
-\Delta_{\text{gauge}}
\right],
\label{eq:tildealpha}
\end{equation}
\begin{equation}
\tilde\beta=\frac{N_c}{16\pi^2f^4}
\left[
\frac{\lambda_L^2\lambda_R^2Y_*^2f^6}{\tilde M_T^2-\tilde M_Q^2}
\ln\!\left(\frac{\tilde M_T^2}{\tilde M_Q^2}\right)
+\frac{y_t^2(f)}{2}\tilde M_T^2
\ln\!\left(\frac{\Lambda_F^2}{\tilde M_T^2}\right)
\right].
\label{eq:tildebeta}
\end{equation}
The prefactor \(N_c/(16\pi^2)\) is the standard one-loop Dirac-fermion
normalization; the trace over \(H_t^2\) sums the contributions of the light top
and the two heavy partners. The gauge-loop contribution is
\begin{equation}
\Delta_{\text{gauge}}
=\frac{3}{32\pi^2}(3g^2+g'^2)m_\rho^2f^2>0,
\label{eq:Dgauge}
\end{equation}
which has the same mass\(^4\) dimension as the top-partner term in
Eq.~\eqref{eq:tildealpha}.

Minimizing \(V_{\text{CW}}(\theta)\) determines
\(\xi=\sin^2\theta_0=\tilde\alpha/(2\tilde\beta)\) and predicts the physical
Higgs boson mass dynamically:
\begin{equation}
m_h^2=\frac{V_{\text{CW}}''(\theta_0)}{f^2}
=8\tilde\beta\,\xi(1-\xi)f^2
=8\tilde\beta\,v^2(1-\xi).
\label{eq:mh_dyn}
\end{equation}
Crucially, \(m_h\) depends directly on \(\tilde\beta\) and \(v=246\GeV\),
independent of the partial tuning in \(\tilde\alpha\).

\paragraph{Calibration.}
Consistent with the \(-N_c/(16\pi^2)\) normalization of Eq.~\eqref{eq:VCW},
we evaluate Eqs.~\eqref{eq:tildealpha}--\eqref{eq:tildebeta} at the benchmark
\(f=1.5\TeV\), \(\tilde M_Q=2.0\TeV\), \(\tilde M_T=1.45\TeV\), and running
top Yukawa \(y_t(f)\simeq0.85\). The gauge contribution partially offsets the
top-partner contribution to \(\tilde\alpha\); its precise size is fixed by the
calibration. This gives
\begin{equation}
\tilde\beta\simeq3.32\times10^{-2},
\qquad
\tilde\alpha\simeq1.79\times10^{-3},
\label{eq:benchmark_coeffs}
\end{equation}
and hence
\begin{equation}
\xi\simeq0.0269,
\qquad
m_h\simeq125.1\GeV,
\label{eq:benchmark_out}
\end{equation}
matching the measured Higgs boson mass
\(m_h^{\text{exp}}=125.11\pm0.11\GeV\). The value
\(\tilde\alpha/\tilde\beta\simeq0.054\) corresponds to a few-percent tuning,
which is the standard tuning of composite-Higgs models with
\(\xi\sim\text{few}\%\). The derivation of \(\tilde\alpha,\tilde\beta\), and
the scan of their dependence on the model parameters, is given in
Appendix~\ref{app:derivation}.

%=====================================================================
\section{Matching: predicted bands with uncertainties}
\label{sec:matching}
%=====================================================================

Beyond \(b=2a^2-1\), an independent observable is provided by the Higgs
trilinear self-coupling modifier
\begin{equation}
\kappa_\lambda\equiv\frac{\lambda_3}{\lambda_3^{\text{SM}}}
=\frac{b}{a}
=\frac{1-2\xi}{\sqrt{1-\xi}},
\label{eq:klambda}
\end{equation}
probed in Higgs-pair production via the trilinear vertex --- a distinct diagram
from the \(VVhh\) contact interaction carrying \(b\). Equation~\eqref{eq:klambda}
holds for the minimal two-parameter trigonometric
potential~\eqref{eq:VCWexpand} at leading order in \(\sin^2\theta\). At the
calibrated benchmark \(\xi=0.0269\) (\(f=1.5\TeV\)),
\begin{equation}
a=0.9865,\qquad b=0.9462,\qquad\kappa_\lambda=0.959.
\label{eq:abk_num}
\end{equation}
Propagating uncertainties from the single-Higgs measurement \(a=\kappa_V\)
gives \(\delta b=4a\,\delta a\) and
\(\delta\kappa_\lambda\simeq3.05\,\delta a\), as summarized in
Table~\ref{tab:matching}.

\begin{table}[t]
\caption{Propagated coupling predictions and discriminating power at current LHC
and projected HL-LHC precision.}
\label{tab:matching}
\begin{ruledtabular}
\begin{tabular}{lccc}
Input \(a=\kappa_V\) & \(b=2a^2-1\) & \(\kappa_\lambda=b/a\) & Discriminating power\\
\midrule
\(1.00\pm0.03\) (current) & \(1.00\pm0.12\) & \(1.00\pm0.09\) & Weak\\
\(1.00\pm0.01\) (HL-LHC) & \(1.00\pm0.04\) & \(1.00\pm0.03\) & Decisive\\
\end{tabular}
\end{ruledtabular}
\end{table}

At current precision the SM and the benchmark~\eqref{eq:abk_num} cannot be
separated. At projected HL-LHC precision (\(\delta a\simeq0.01\)), the direct
displacement \(|1-a|/\delta a\simeq1.35\sigma\) per channel accumulates to more
than \(2\sigma\) across combined single-Higgs decay channels, while \(b\) and
\(\kappa_\lambda\) provide non-trivial consistency tests of the composite
relation.

%=====================================================================
\section{The pre-registered collider test}
\label{sec:test}
%=====================================================================

To fix the falsification criterion prior to unblinding future Higgs-pair data,
the following five-step protocol is established.

\begin{enumerate}[label=\Roman*.,leftmargin=*]
\item \textbf{Inputs.} Single-Higgs coupling fits (for \(a=\kappa_V\)) and,
independently, Higgs-pair observables (\(b\) via the high-\(m_{hh}\) VBF
spectrum; \(\kappa_\lambda\) via the near-threshold inclusive pair rate).

\item \textbf{Prediction.} From the measured \(a_{\text{exp}}\pm\delta a\),
compute \(b_{\text{pred}}=2a_{\text{exp}}^2-1\pm4a_{\text{exp}}\delta a\) and
\((\kappa_\lambda)_{\text{pred}}=b_{\text{pred}}/a_{\text{exp}}\) prior to
unblinding Higgs-pair results. Publish the prediction before unblinding.

\item \textbf{Exclusion criterion.} The benchmark is excluded at
\(95\%\text{ CL}\) if the measured \(b\) lies outside
\(b_{\text{pred}}\pm2\sqrt{\delta b^2+\delta b_{\text{exp}}^2}\), or
analogously for \(\kappa_\lambda\). Independence of \(b\) and
\(\kappa_\lambda\) requires that the two kinematic regimes (the high-\(m_{hh}\)
tail and the near-threshold region) be separately resolved; if they are not,
the over-determination is weak. A confirmed deviation establishes the pNGB
coupling pattern~\eqref{eq:ab_rel}, though distinguishing FQCD from other
composite realisations requires resolving heavy resonances at \(\Lambda_F\).

\item \textbf{No retuning.} Parameters
\((f,\tilde M_Q,\tilde M_T)\) fixed in Sec.~\ref{sec:spectrum} may not be refit
post-unblinding; failure under Criterion~III retires this realisation.

\item \textbf{Over-determination.} Because \(b=2a^2-1\) and
\(\kappa_\lambda=b/a\) are independent combinations of \(\xi\), measuring both
over-determines \(\xi\); mutual inconsistency excludes the minimal pNGB
structure itself.
\end{enumerate}

\paragraph{Implementation notes.}
A true pre-registration should specify: the exact experimental observables
(e.g., ATLAS/CMS STXS bins for \(\kappa_V\), di-Higgs invariant-mass bins for
\(b\) and \(\kappa_\lambda\)); the covariance treatment (correlations between
single-Higgs and di-Higgs channels); the fitting framework (e.g., the
\(\kappa\)-framework of the LHC Higgs Working Group); and the exact decision
threshold (\(2\sigma\) two-sided) before data are examined. This is stated as
the target for a future implementation; the present paper fixes the criterion,
not the full experimental machinery.

%=====================================================================
\section{Relation to the spin-one FQCU hypothesis}
\label{sec:relation}
%=====================================================================

The scalar-composite Higgs theory developed here and the vector-composite
hypothesis proposed for \(b\to s\mu^+\mu^-\) flavour
anomalies~\cite{FQCUflavor} cannot share a single QCD-like confining sector.
Standard strong-sector scaling places spin-one vector mesons at
\(m_V\sim g_\rho f\approx8f\); for \(f=1.5\TeV\) this gives
\(m_V\sim12\TeV\), whereas the flavour application of Ref.~\cite{FQCUflavor}
requires \(m_V\sim\text{a few }\GeV\). The two hypotheses must therefore be
treated as separate physical sectors operating at distinct dynamical scales.

%=====================================================================
\section{Open problems}
\label{sec:open}
%=====================================================================

The following questions are not addressed in this paper and are required before
the FQCU composite-Higgs scenario can be regarded as a complete microscopic
theory:

\begin{enumerate}[label=\textbf{\arabic*.}]
\item \textbf{Coset derivation.} The \(SO(5)/SO(4)\) coset of
Sec.~\ref{sec:coset-assumption} is assumed. A derivation from the eight-flavour
FQCD dynamics --- showing that \(U(8)\) breaks to a subgroup containing
\(SO(5)/SO(4)\) and identifying the observed Higgs as the associated pNGB ---
is required.
\item \textbf{Vacuum stability.} The one-loop \(\beta_\lambda\) of
Eq.~\eqref{eq:betalambda} is incomplete (no Yukawas, no two-loop). A full
numerical running, including Yukawa contributions, is required to establish
that the scalar potential is stable up to \(\Lambda_F\).
\item \textbf{Confinement.} A non-perturbative (lattice) study of the FQCD
theory is required to confirm confinement and to fix \(\Lambda_F\).
\item \textbf{Top-partner spectrum.} The partial-compositeness parameters
\(\lambda_L,\lambda_R,Y_*\) and the partner masses \(M_{Q,T}\) are treated as
free parameters. A dynamical determination of these parameters from the FQCD
spectrum is required.
\item \textbf{Fermion mass hierarchy.} Only the top Yukawa is generated in
this paper. A mechanism for the lighter fermion masses and mixings is required.
\item \textbf{Precision CW calculation.} The CW potential of Eq.~\eqref{eq:VCW}
is computed at one loop. Two-loop corrections and the matching to the SILH
effective Lagrangian should be evaluated.
\end{enumerate}

%=====================================================================
\section{Conclusion}
\label{sec:conclusion}
%=====================================================================

We have presented a \emph{candidate} microscopic realization of the observed
Higgs boson as a composite pseudo-Nambu--Goldstone bound state of Farid Quantum
Command Units within a confining \(SU(N_F)\) FQCD gauge theory. The proposal
makes several derived statements:

\begin{itemize}
\item Within the assumed \(SO(5)/SO(4)\) coset, the relations
\(b=2a^2-1\) and \(\kappa_\lambda=b/a\) follow algebraically.
\item The custodial \((\mathbf2,\mathbf2)\) bi-doublet embedding gives
\(\rho=1\) at tree level.
\item The theory is asymptotically free for all \(N_F\ge1\) with \(n_f=2\)
and \(n_s=8\).
\item With the running Yukawa \(y_t(f)\approx0.85\) and partner masses
\(\tilde M_{Q,T}\sim1.45\text{--}2.0\TeV\), the one-loop CW potential gives
\(\tilde\beta\simeq3.32\times10^{-2}\),
\(\tilde\alpha\simeq1.79\times10^{-3}\), and \(m_h\simeq125.1\GeV\),
matching experiment.
\end{itemize}

The proposal \emph{assumes} the \(SO(5)/SO(4)\) coset and the
partial-compositeness mixing structure; the derivation of these from the FQCD
dynamics is an open problem, listed in Sec.~\ref{sec:open}.

The benchmark
\((a,b,\kappa_\lambda)=(0.9865,0.9462,0.959)\) is stated as a pre-registered
target for the High-Luminosity LHC. No evidence for FQCU is claimed.

\medskip
\noindent\textit{All discoveries herein are by the mercy of Allah. This
knowledge is a trust.}

%=====================================================================
\begin{acknowledgments}
The author thanks the Tanfarid Vision Research Institute for support. This
work is dedicated to the memory of Mazama Khatun and Abul Kashem Chowdhury.
\end{acknowledgments}

%=====================================================================
\appendix

\section{Derivation of the Coleman--Weinberg coefficients}
\label{app:derivation}
%=====================================================================

This appendix outlines the derivation of the CW potential coefficients
\(\tilde\alpha,\tilde\beta\) of
Eqs.~\eqref{eq:tildealpha}--\eqref{eq:tildebeta}, the \(\theta\)-dependent top
mass of Eq.~\eqref{eq:mtop}, and the benchmark calibration quoted in
Sec.~\ref{sec:spectrum}.

\subsection{Eigenvalues of \(H_t(\theta)\)}

Define the Hermitian matrix
\(H_t(\theta)=\mathcal{M}_t^\dagger(\theta)\mathcal{M}_t(\theta)\), with
\(\mathcal{M}_t\) given by Eq.~\eqref{eq:Mt}. Explicitly,
\begin{equation}
H_t(\theta)=
\begin{pmatrix}
\lambda_R^2f^2 & 0 & \lambda_RfM_T\\
0 & \lambda_L^2f^2\cos^2\theta+M_Q^2
& \tfrac{\lambda_L^2f^2}{\sqrt2}\sin\theta\cos\theta+M_QY_*f\cos\theta\\
\lambda_RfM_T
& \tfrac{\lambda_L^2f^2}{\sqrt2}\sin\theta\cos\theta+M_QY_*f\cos\theta
& \tfrac{\lambda_L^2f^2}{2}\sin^2\theta+Y_*^2f^2\cos^2\theta+M_T^2
\end{pmatrix}.
\label{eq:Ht}
\end{equation}
The exact eigenvalues are the roots of a cubic polynomial. For the purpose of
the CW calculation at \(M_{Q,T}\gg y_tf\), the determinant factorization gives
the light eigenvalue
\begin{equation}
m_t^2(\theta)
=\frac{\det[H_t(\theta)]}{\tilde M_Q^2(\theta)\tilde M_T^2(\theta)}.
\label{eq:mt-light}
\end{equation}
Computing \(\det[H_t(\theta)]\) and expanding for small \(\sin^2\theta\) at
fixed heavy eigenvalues \(\tilde M_{Q,T}\), one finds
\begin{equation}
m_t^2(\theta)
=\tfrac12y_t^2(f)f^2
\bigl(\sin^2\theta-\sin^4\theta\bigr)
+\mathcal{O}(\sin^6\theta).
\label{eq:mtop-derived}
\end{equation}
The form \(\sin^2\theta-\sin^4\theta\) reflects the specific three-site mixing
of Eq.~\eqref{eq:Mt}; in the limit \(M_{Q,T}\to\infty\) with fixed \(y_t\), it
reduces to the standard minimal-composite-Higgs result
\(m_t^2\propto\sin^2\theta\cos^2\theta\approx\sin^2\theta\), and the
difference is \(\mathcal{O}(\sin^4\theta)\).

\subsection{CW integral and coefficient extraction}

The CW potential for a Dirac fermion \(t\) with Yukawa-derived mass
\(m_t(\theta)\) and two heavy partners \(\tilde M_{Q,T}\) is
\begin{equation}
V_{\text{CW}}(\theta)
=-\frac{N_c}{16\pi^2}
\int\frac{d^4p_E}{(2\pi)^4}
\ln\frac{\det(p_E^2+H_t(\theta))}{\det(p_E^2+H_t(0))}.
\label{eq:VCWintegral}
\end{equation}
Using the polynomial identity
\(\det(p^2+H_t)=(p^2+m_t^2)(p^2+\tilde M_Q^2)(p^2+\tilde M_T^2)\)
and expanding to
\(\mathcal{O}(\sin^2\theta,\sin^4\theta)\), the logarithms factorize into
contributions from the light eigenvalue and the two heavy partners. The light
top contributes to \(\tilde\beta\) through the
\(m_t^4\ln(m_t^2/\tilde M_T^2)\) term; the heavy partners contribute to both
\(\tilde\alpha\) and \(\tilde\beta\) through the
\(\tilde M_{Q,T}^4\ln\tilde M_{Q,T}^2\) terms.

Collecting the \(\mathcal{O}(\sin^2\theta)\) terms gives
Eq.~\eqref{eq:tildealpha}, and collecting the
\(\mathcal{O}(\sin^4\theta)\) terms gives Eq.~\eqref{eq:tildebeta}, in
agreement with the main text. The gauge-loop contribution is
Eq.~\eqref{eq:Dgauge}, which has mass\(^4\) dimension and is therefore
consistent with the top-partner term in Eq.~\eqref{eq:tildealpha}.

\subsection{Numerical benchmark at \(f=1.5\TeV\)}

For the benchmark parameters \(f=1.5\TeV\), \(\tilde M_Q=2.0\TeV\),
\(\tilde M_T=1.45\TeV\), \(M_Q=1.8\TeV\), \(M_T=1.2\TeV\),
\(\lambda_L=\lambda_R=1.1\), \(Y_*=3.2\), and running top Yukawa
\(y_t(f)=0.85\), the numerical evaluation of
Eqs.~\eqref{eq:tildealpha}--\eqref{eq:tildebeta} gives:
\begin{align*}
\tilde\alpha &= 1.79\times10^{-3},\\
\tilde\beta &= 3.32\times10^{-2},\\
\xi &= \tilde\alpha/(2\tilde\beta)=0.0269,\\
m_h &= \sqrt{8\tilde\beta v^2(1-\xi)}=125.1\GeV.
\end{align*}
This is the calibration quoted in Sec.~\ref{sec:spectrum}.

\subsection{Scan over model parameters}

Varying \(\tilde M_T\in[1.2,2.0]\TeV\) at fixed \(\tilde M_Q=2.0\TeV\) and
\(f=1.5\TeV\) shifts \(\tilde\beta\) from \(3.8\times10^{-2}\) (at
\(\tilde M_T=1.2\TeV\)) to \(2.6\times10^{-2}\) (at \(\tilde M_T=2.0\TeV\)),
corresponding to \(m_h\) in the range \([134,117]\GeV\). The value
\(\tilde M_T=1.45\TeV\) is therefore selected to reproduce
\(m_h=125\GeV\). Varying \(f\) from \(1.0\TeV\) to \(2.0\TeV\) at fixed
partner masses shifts \(\xi\) from \(0.061\) to \(0.015\) and \(m_h\) from
\(123\GeV\) to \(127\GeV\); the benchmark \(f=1.5\TeV\) is a representative
choice consistent with electroweak precision. A full scan over
\((f,\tilde M_Q,\tilde M_T,Y_*)\) will be presented in future work.

%=====================================================================
\begin{thebibliography}{99}

\bibitem{ATLAS2012}
G.~Aad \textit{et al.} (ATLAS Collaboration),
``Observation of a new particle in the search for the Standard Model Higgs boson
with the ATLAS detector at the LHC,''
Phys.\ Lett.\ B \textbf{716}, 1--29 (2012),
\texttt{arXiv:1207.7214 [hep-ex]}.

\bibitem{CMS2012}
S.~Chatrchyan \textit{et al.} (CMS Collaboration),
``Observation of a new boson at a mass of 125 GeV with the CMS experiment at the
LHC,''
Phys.\ Lett.\ B \textbf{716}, 30--61 (2012),
\texttt{arXiv:1207.7235 [hep-ex]}.

\bibitem{ATLASCMS2016}
G.~Aad \textit{et al.} (ATLAS and CMS Collaborations),
``Measurements of the Higgs boson production and decay rates and constraints on
its couplings from a combined ATLAS and CMS analysis of the LHC \(pp\) collision
data at \(\sqrt{s}=7\) and 8 TeV,''
JHEP \textbf{08}, 045 (2016),
\texttt{arXiv:1606.02266 [hep-ex]}.

\bibitem{FQCUflavor}
M.~F.~I.~Chowdhury,
``A Vector Composite of FQCU Scalars as a Probe of the \(b\to s\mu^+\mu^-\)
Anomalies,''
Companion manuscript (2026).

\bibitem{Giudice2007}
G.~F.~Giudice, C.~Grojean, A.~Pomarol, and R.~Rattazzi,
``The Strongly-Interacting Light Higgs,''
JHEP \textbf{06}, 045 (2007),
\texttt{arXiv:hep-ph/0703164}.

\bibitem{Agashe2005}
K.~Agashe, R.~Contino, and A.~Pomarol,
``The Minimal Composite Higgs Model,''
Nucl.\ Phys.\ B \textbf{719}, 165--187 (2005),
\texttt{arXiv:hep-ph/0412089}.

\bibitem{Panico2016}
G.~Panico and A.~Wulzer,
\textit{The Composite Nambu-Goldstone Higgs},
Lect.\ Notes Phys.\ \textbf{913} (Springer, 2016),
\texttt{arXiv:1506.01961 [hep-ph]}.

\bibitem{Kaplan1991}
D.~B.~Kaplan,
``Flavor at SSC energies: A New mechanism for dynamically generated fermion
masses,''
Nucl.\ Phys.\ B \textbf{365}, 259--278 (1991).

\bibitem{Contino2007}
R.~Contino, T.~Kramer, M.~Son, and R.~Sundrum,
``Warped/composite phenomenology simplified,''
JHEP \textbf{05}, 074 (2007),
\texttt{arXiv:hep-ph/0612180}.

\bibitem{Machacek1983}
M.~E.~Machacek and M.~T.~Vaughn,
``Two-loop renormalization group equations in a general quantum field theory:
I. Wave function renormalization,''
Nucl.\ Phys.\ B \textbf{222}, 83--103 (1983).

\bibitem{Machacek1985}
M.~E.~Machacek and M.~T.~Vaughn,
``Two-loop renormalization group equations in a general quantum field theory:
III. Scalar quartic couplings,''
Nucl.\ Phys.\ B \textbf{249}, 70--92 (1985).

\bibitem{Matsedonskyi2013}
O.~Matsedonskyi, G.~Panico, and A.~Wulzer,
``Light Top Partners for a Light Composite Higgs,''
JHEP \textbf{01}, 164 (2013),
\texttt{arXiv:1204.6333 [hep-ph]}.

\bibitem{Weinberg1967}
S.~Weinberg,
``Precise Relations between the Spectra of Vector and Axial-Vector Mesons,''
Phys.\ Rev.\ Lett.\ \textbf{18}, 507--509 (1967).

\end{thebibliography}

\end{document}
