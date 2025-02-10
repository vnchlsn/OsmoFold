# Welcome to OsmoFold v.0.4.0

Vincent Nicholson<sup>1</sup> and Thomas C. Boothby<sup>1</sup>

<sup>1</sup> - Department of Molecular Biology, Univeristy of Wyoming, Laramie WY

OsmoFold is a tool that uses Auton and Bolen's empherical transfer free energy measurements to calculate the effect of osmolytes on protein folding. If you haven't already done so, we recommend reading the following literature.


1.   [Application of the transfer model to understand how naturally occurring osmolytes affect protein stability](https://pubmed.ncbi.nlm.nih.gov/17875431/)
2.   [Predicting the energetics of osmolyte-induced protein folding/unfolding](https://pubmed.ncbi.nlm.nih.gov/16214887/)
3.   [Additive transfer free energies of the peptide backbone unit that are independent of the model compound and the choice of concentration scale](https://pubmed.ncbi.nlm.nih.gov/14756570/)
4.   [Its Preferential Interactions with Biopolymers Account for Diverse Observed Effects of Trehalose](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4572414/)

Below are some recent papers, preprints, and reviews that cover this method or an adapted version of it. They might be useful for your understanding of the protocol.

1.   [Disordered proteins interact with the chemical environment to tune their
protective function during drying](https://elifesciences.org/reviewed-preprints/97231)
2.   [LEA_4 motifs function alone and in conjunction with synergistic cosolutes to protect a labile enzyme during desiccation](https://www.biorxiv.org/content/10.1101/2024.09.04.611296v1.full.pdf)
3.   [Osmolyte-IDP Interactions During Desiccation](https://www.sciencedirect.com/science/article/pii/S1877117324001765?via%3Dihub)

In short, OsmoFold will quantitatively asses the impact of several omsolytes on a given protein conformational change. Please note that OsmoFold does **not guarentee that said conformational change is physiologically relevant.** To put it simply, OsmoFold requires the existence of a known, binary conformational change in the protein(s) of interest. Additionally, it cannot rule out the presence of specific protein small molecule interactions, such as direct binding. In general, OsmoFold should be used to support experimental data, not replace it.
