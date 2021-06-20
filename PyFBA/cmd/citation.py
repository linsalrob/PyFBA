"""
Provide the citations! Make it easy for people!
"""
import argparse
import sys

def preamble():
    """
    Some text to thank people for citing us!
    """
    print(f"""
Thank you for citing PyFBA and the ModelSEED!
The two papers by Cuevas et al. are for PyFBA, and the Henry et al. paper describes the ModelSEED that PyFBA depend 
upon. 
    """)


def text_citations():
    """
    Provide the citations as plain text
    :return:
    """

    return """
Cuevas DA, Garza D, Sanchez SE, Rostron J, Henry CS, Vonstein V, Overbeek RA, Segall A, Rohwer F,  Dinsdale EA, Edwards RA. 2014.
Elucidating genomic gaps using phenotypic profiles. F1000Research. 3:210
doi: 10.12688/f1000research.5140.2
https://f1000research.com/articles/3-210

Cuevas, Daniel A., Janaka Edirisinghe, Chris S. Henry, Ross Overbeek, Taylor G. O’Connell, and Robert A. Edwards. 2016.
From DNA to FBA: How to Build Your Own Genome-Scale Metabolic Model.
Frontiers in Microbiology 7 (June): 907.
http://journal.frontiersin.org/article/10.3389/fmicb.2016.00907/full

Henry CS, DeJongh M, Best AA, Frybarger PM, Linsay B, Stevens RL. 2010. 
High-throughput generation, optimization and analysis of genome-scale metabolic models. 
Nat Biotechnol 28:977–982.
https://www.nature.com/articles/nbt.1672
    """


def bibtex():
    """
    Provide the citations in bibtex format
    :return:
    """

    return """

@ARTICLE{Cuevas2016-wg,
  title    = "Elucidating genomic gaps using phenotypic profiles",
  author   = "Cuevas, Daniel A and Garza, Daniel and Sanchez, Savannah E and
              Rostron, Jason and Henry, Chris S and Vonstein, Veronika and
              Overbeek, Ross A and Segall, Anca and Rohwer, Forest and
              Dinsdale, Elizabeth A and Edwards, Robert A",
  abstract = "Advances in genomic sequencing provide the ability to model the
              metabolism of organisms from their genome annotation. The
              bioinformatics tools developed to deduce gene function through
              homology-based methods are dependent on public databases; thus,
              novel discoveries are not readily extrapolated from current
              analysis tools with a homology dependence. Multi-phenotype Assay
              Plates (MAPs) provide a high-throughput method to profile
              bacterial phenotypes by growing bacteria in various growth
              conditions, simultaneously. More robust and accurate
              computational models can be constructed by coupling MAPs with
              current genomic annotation methods. PMAnalyzer is an online tool
              that analyzes bacterial growth curves from the MAP system which
              are then used to optimize metabolic models during in silico
              growth simulations. Using Citrobacter sedlakii as a prototype,
              the Rapid Annotation using Subsystem Technology (RAST) tool
              produced a model consisting of 1,367 enzymatic reactions. After
              the optimization, 44 reactions were added to, or modified within,
              the model. The model correctly predicted the outcome on 93\% of
              growth experiments.",
  journal  = "F1000Res.",
  volume   =  3,
  month    =  oct,
  year     =  2016,
  url      = "https://f1000research.com/articles/3-210/v2/pdf",
  language = "en",
  issn     = "2046-1402",
  doi      = "10.12688/f1000research.5140.2"
}

@ARTICLE{Cuevas2016-si,
  title     = "From {DNA} to {FBA}: How to Build Your Own {Genome-Scale}
               Metabolic Model",
  author    = "Cuevas, Daniel A and Edirisinghe, Janaka and Henry, Chris S and
               Overbeek, Ross and O'Connell, Taylor G and Edwards, Robert A",
  abstract  = "Microbiological studies are increasingly relying on in silico
               methods to perform exploration and rapid analysis of genomic
               data, and functional genomics studies are supplemented by the
               new perspectives that genome-scale metabolic models offer. A
               mathematical model consisting of a microbe's entire metabolic
               map can be rapidly determined from whole-genome sequencing and
               annotating the genomic material encoded in its DNA. Flux-balance
               analysis (FBA), a linear programming technique that uses
               metabolic models to predict the phenotypic responses imposed by
               environmental elements and factors, is the leading method to
               simulate and manipulate cellular growth in silico. However, the
               process of creating an accurate model to use in FBA consists of
               a series of steps involving a multitude of connections between
               bioinformatics databases, enzyme resources, and metabolic
               pathways. We present the methodology and procedure to obtain a
               metabolic model using PyFBA, an extensible Python-based
               open-source software package aimed to provide a platform where
               functional annotations are used to build metabolic models
               (http://linsalrob.github.io/PyFBA). Backed by the Model SEED
               biochemistry database, PyFBA contains methods to reconstruct a
               microbe's metabolic map, run FBA upon different media
               conditions, and gap-fill its metabolism. The extensibility of
               PyFBA facilitates novel techniques in creating accurate
               genome-scale metabolic models.",
  journal   = "Front. Microbiol.",
  publisher = "Frontiers Media SA",
  volume    =  7,
  pages     = "907",
  month     =  jun,
  year      =  2016,
  url       = "http://dx.doi.org/10.3389/fmicb.2016.00907",
  keywords  = "flux-balance analysis; genome annotation; in silico modeling;
               metabolic modeling; metabolic reconstruction; model SEED",
  language  = "en",
  issn      = "1664-302X",
  pmid      = "27379044",
  doi       = "10.3389/fmicb.2016.00907",
  pmc       = "PMC4911401"
}

@ARTICLE{Henry2010-zo,
  title    = "High-throughput generation, optimization and analysis of
              genome-scale metabolic models",
  author   = "Henry, Christopher S and DeJongh, Matthew and Best, Aaron A and
              Frybarger, Paul M and Linsay, Ben and Stevens, Rick L",
  abstract = "Genome-scale metabolic models have proven to be valuable for
              predicting organism phenotypes from genotypes. Yet efforts to
              develop new models are failing to keep pace with genome
              sequencing. To address this problem, we introduce the Model SEED,
              a web-based resource for high-throughput generation, optimization
              and analysis of genome-scale metabolic models. The Model SEED
              integrates existing methods and introduces techniques to automate
              nearly every step of this process, taking approximately 48 h to
              reconstruct a metabolic model from an assembled genome sequence.
              We apply this resource to generate 130 genome-scale metabolic
              models representing a taxonomically diverse set of bacteria.
              Twenty-two of the models were validated against available gene
              essentiality and Biolog data, with the average model accuracy
              determined to be 66\% before optimization and 87\% after
              optimization.",
  journal  = "Nat. Biotechnol.",
  volume   =  28,
  number   =  9,
  pages    = "977--982",
  month    =  sep,
  year     =  2010,
  url      = "http://dx.doi.org/10.1038/nbt.1672",
  language = "en",
  issn     = "1087-0156, 1546-1696",
  pmid     = "20802497",
  doi      = "10.1038/nbt.1672"
}


    """


def other_citations():
    """
    Other citations that should be included. A dict of strings
    :return:
    :rtype:
    """

    return {
    'conda-forge': """
conda-forge community. (2015). 
The conda-forge Project: Community-based Software Distribution Built on the conda Package Format and Ecosystem. 
Zenodo. http://doi.org/10.5281/zenodo.4774216
    """
    }

def other_citations_bibtex():
    """
    Other citations that should be included. A dict of strings
    :return:
    :rtype:
    """

    return {
        'conda-forge': """
@misc{conda_forge_community_2015_4774216,
  author       = {conda-forge community},
  title        = {{The conda-forge Project: Community-based Software
                   Distribution Built on the conda Package Format and
                   Ecosystem}},
  month        = jul,
  year         = 2015,
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.4774216},
  url          = {https://doi.org/10.5281/zenodo.4774216}
}
    }
    """
    }

def cite_me_please():
    """
    Parse the arguments and list the citations.
    """

    parser = argparse.ArgumentParser(description='Cite PyFBA Flux Balance Analysis')
    parser.add_argument('-o', '--output', help='write the citations to a file')
    parser.add_argument('-b', '--bibtex', help='write the citations in bibtex format',
                        action='store_true')
    args = parser.parse_args(sys.argv[2:])

    if args.bibtex:
        cits = bibtex()
        other_cits = other_citations_bibtex()
    else:
        cits = text_citations()
        other_cits = other_citations()

    preamble()

    if args.output:
        with open(args.output, 'w') as out:
            out.write(f"{cits}\n")
            for f in other_cits:
                out.write(f"{other_cits[f]}\n")
    else:
        print(cits)
        print("Additionally, we use these resources, so please cite them:")
        for f in other_cits:
            print(f"{f}:\n{other_cits[f]}\n")
