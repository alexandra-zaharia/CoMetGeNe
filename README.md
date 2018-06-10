# CoMetGeNe
A bioinformatics tool for the exploration of **Co**nserved **Met**abolic and **Ge**nomic **Ne**ighborhoods

Project website: [cometgene.lri.fr](https://cometgene.lri.fr)

Publication (in preparation): Zaharia, Labedan, Froidevaux, Denise. CoMetGeNe: mining conserved neighborhood patterns in metabolic and genomic contexts. *BMC Bioinformatics*. 2018 

## Background

Genomic data and chemical reactions embody the dual aspect of metabolism that allows exploring the links between genome evolution and chemical evolution of enzyme-catalyzed reactions. 

It is well-known that neighboring reactions corresponding to neighboring genes underline an evolutionary advantage in keeping the genes involved in succeeding reactions in close proximity. Finding almost identical sequences of reactions being catalyzed by products of neighboring genes in various species suggests that such sequences are made up of key enzymatic steps, best performed when their encoding genes are adjacent and co-transcribed. This type of metabolic and genomic organization strongly suggests the various species have been under strong evolutionary pressure to optimize the expression of enzyme-coding genes involved in successive reactions.

## What is `CoMetGeNe`?

`CoMetGeNe` is an exploration tool for conserved organizational motifs at both the metabolic and genomic levels:

* In a first step named **trail finding**, `CoMetGeNe` detects chains of metabolic reactions for a given species that are catalyzed by products of neighboring genes.
* In a second (and optional) step named **trail grouping**, `CoMetGeNe` compares trail finding results for several species in order to identify conserved metabolic and genomic patterns. This step allows to answer specific questions such as *Which genes involved in the biosynthesis of peptidoglycan are neighbors in the bacterial species in my data set?* or more general questions such as *Which reactions are generally performed by enzymes encoded by neighboring genes?*

`CoMetGeNe` automatically retrieves metabolic pathway maps and genomic information for the species of your choice from the [KEGG](https://www.kegg.jp) (Kyoto Encyclopedia of Genes and Genomes) knowledge base. `CoMetGeNe` usage is as simple as telling it on which query species to run. In order to do so, one simply determines the three- or four-letters KEGG code for the target organism from the KEGG [catalog](https://www.genome.jp/kegg/catalog/org_list.html) of complete genomes.

## Trail finding

`CoMetGeNe` models a metabolic pathway as a directed graph in which vertices are reactions. In this directed graph, an arc from a vertex *r<sub>i</sub>* to a vertex *r<sub>j</sub>* means that reaction *r<sub>i</sub>* produces a metabolite that serves as substrate to reaction *r<sub>j</sub>*.

Roughly speaking, a **trail** is a path in a graph that can include repeated vertices, but not repeated arcs. `CoMetGeNe` identifies trails of metabolic reactions such that the reactions in the trails are catalyzed by enzymes encoded by neighboring genes. More specifically, for every arc (*r<sub>i</sub>*, *r<sub>j</sub>*) in the graph modeling a metabolic pathway, `CoMetGeNe` attempts to determine a trail *T* of reactions passing through the arc (*r<sub>i</sub>*, *r<sub>j</sub>*) such that

1. the genes involved in the reactions of *T* are neighbors, and
1. no other trail *T'* passing trhough (*r<sub>i</sub>*, *r<sub>j</sub>*) such that the genes involved in reactions of *T'* are neighbors contains more unique reactions than *T*.

Property 2 above refers to a concept termed **span**: the span of a trail *T* represents the number of unique reactions in *T*.

Flexibility is allowed in the definition of neighborhood: `CoMetGeNe` is able to skip a few reactions and/or genes.

Trail finding is performed using the `CoMetGeNe.py` Python script.

### Example

Suppose `CoMetGeNe.py` is executed as follows:

```
python2 CoMetGeNe.py eco data/eco/ -dG 2 -dD 1 -o eco.out
```

Metabolic pathway maps for *Escherichia coli* K-12 MG1655 (`eco`) are automatically downloaded from KEGG and stored under `data/eco/`. At most two genes (`-dG 2`) and one reaction (`-dD 1`) can be skipped. Results are saved in the output file `eco.out`.

### Output files

For the above example, `CoMetGeNe` identifies 238 trails of span ranging from 2 to 10 (i.e., the 238 trails contain from 2 to 10 unique metabolic reactions). Below is the output corresponding to a trail of span 3:

```
...
path_eco00564.kgml: Found a trail of span 3 containing skipped vertices
path_eco00564.kgml: 110 -> [139] -> 104 -> 123
path_eco00564.kgml: R02054 -> [R04864] -> R02053 -> R03416
path_eco00564.kgml: eco:b3821 -> [eco:b2836] -> eco:b3821 -> eco:b3825
path_eco00564.kgml: 3.1.1.32 -> [2.3.1.40] -> 3.1.1.4 -> 3.1.1.5
path_eco00564.kgml: Skipped genes: eco:b3823, eco:b3822
...
```

* `path_eco00564.kgml` is the file name for the pathway map `00564` in `eco` (glycerophospholipid metabolism), retrieved automatically from KEGG by `CoMetGeNe`.
* The four lines with entities separated by arrows (`->`) represent the trail `R02054 -> R02053 -> R03416` in four distinct manners, using:
  * The KGML identifiers of the reactions in the trail (`110 -> 104 -> 123`).
  * The KEGG R numbers associated to the reactions in the trail (`R02054 -> R02053 -> R03416`). Span is computed in terms of distinct R numbers in the trail.
  * The names of genes whose products are involved in reactions in the trail (`eco:b3821 -> eco:b3821 -> eco:b3825`).
  * The EC numbers associated to the reactions in the trail (`3.1.1.32 -> 3.1.1.4 -> 3.1.1.5`).
* The reaction `R04864` was skipped (allowed because `CoMetGeNe.py` was executed with the option `-dD 1`): it is shown in square brackets, along with the corresponding KGML identifier (`139`), associated gene (`b2836`), and EC number (`2.3.1.40`).
* Two genes were skipped (allowed because `CoMetGeNe.py` was executed with the option `-dG 2`): `eco:b3823` and `eco:b3822`.

### Parallel execution

If a large number of species needs to be analyzed, an important speedup can be attained if `CoMetGeNe.py` is ran in parallel. This functionality is provided by the script `CoMetGeNe_launcher.py`:

* Metabolic pathway maps are retrieved from KEGG using 3 threads (the maximum permitted as of June 2018).
* Genomic information is retrieved from KEGG using 2 threads (the maximum permitted as of June 2018).
* Trail finding is performed on the maximum number of available threads (e.g., 4 or 8 for a quad-core CPU).

These values can be adjusted in `CoMetGeNe_launcher.py`, where the data directory for metabolic pathways as well as the species to analyze can also be specified.

## Trail grouping

Once `CoMetGeNe` trails are identified for several species, the conservation of metabolic and genomic organizational motifs can be investigated at an interspecific level. For simplicity, trails of metabolic reactions catalyzed by products of neighboring genes for a given species are called **metabolic and genomic patterns**. The role of trail grouping is to identify **conserved** such patterns for several species.

Given a reference species *S* among the ones for which trail grouping has been performed (using `CoMetGeNe.py` or `CoMetGeNe_launcher.py`), the script `grouping.py` can exploit trails of the reference species in either of the following two ways:

* `CoMetGeNe` trails of *S* are analyzed in terms of genomic conservation across the other species in the data set, referred to as *grouping trails by genes*. This consists in determining whether the genes of *S* involved in these trails have neighboring homologues in other species. See *Genomic conservation patterns* below.

* `CoMetGeNe` trails of *S* are analyzed in terms of metabolic conservation across the other species in the data set, referred to as *grouping trails by reactions*. This consists in determining whether reactions in the `CoMetGeNe` trails of *S* are also performed by products of neighboring geens in other species. See *Metabolic conservation patterns* below.

### Note regarding directory structure

In order to perform trail grouping, the metabolic pathway maps of all species in the data set need to be stored in KGML format in a single directory with subdirectories for every species. The subdirectory names need to be the three- or four-letter KEGG codes for the species in question. For example, a correct directory structure can look like this:

```
data/bsu/path_aae00010.kgml, path_aae00020.kgml, ...
data/pae/path_bbn00010.kgml, path_bbn00020.kgml, ...
data/eco/path_eco00010.kgml, path_eco00020.kgml, ...
data/ype/path_mpn00010.kgml, path_mpn00020.kgml, ...
```

It is important to preserve this type of directory structure if `CoMetGeNe.py` is launched directly; in case `CoMetGeNe_launcher.py` was used, this particular directory structure is ensured. 

### Genomic conservation patterns (grouping by genes)

Trail grouping by genes identifies conservation patterns between a reference species and the other species in the data set in terms of genomic organization.

For example, suppose `grouping.py` is executed as follows:

```
python2 grouping.py genes results/ data/ eco -o tsg_eco.csv
```

This results in detecting genomic conservation patterns (trail grouping by `genes`) for species `eco`, using `CoMetGeNe` results stored in `results/` and metabolic pathway maps stored in KGML format in `data/` (see the *Note regarding directory structure* above). The output of trail grouping by genes is stored in CSV format in the file `tsg_eco.csv`.

The CSV file contains a line for every gene of the reference species *S* involved in `CoMetGeNe` trails of *S* that are common to *S* and at least one other species from the dataset. Groups of neighboring genes in *S* involved in `CoMetGeNe` trails of *S* are separated by the line `***`. In this CSV file, the line for a gene *g* of *S* contains:

* The name of gene *g*.
* The name of the chromosome on which *g* is located.
* The strand on the chromosome on which *g* is located (`+` for the positive strand, `-` for the negative strand).
* A column for every other species in the data set that can take either of the following values:
  * A cross (`x`) if *g* has an homologue in the other species that is a neighbor of at least one other gene involved in the trail;
  * A dot (`.`) if *g* has no such homologue.

**Example:** Suppose trail finding was performed for species `aae`, `bbn`, `eco`, and `mpn`. A small part of the CSV obtained when grouping `CoMetGeNe` trails by genes for `eco` as the reference species is reproduced below (slightly re-formatted for readability purposes):

```
eco_gene;chr;str;aae;bbn;mpn
b0114;   chr; + ; . ; . ; x 
b0115;   chr; + ; . ; . ; x
b0116;   chr; + ; x ; . ; x
```
From the table above, it can be seen that:

* Species `aae` has at least two neighboring homologues to the gene `b0116` in `eco`;
* Species `bbn` has no neighboring homologues for the three genes in `eco`;
* Species `mpn` has neighboring homologues for all the three genes in `eco`.

### Metabolic conserveration patterns (grouping by reactions)

Trail grouping by reactions identifies conservation patterns between a reference species and the other species in the data set in terms of metabolic organization.

For example, suppose `grouping.py` is executed as follows:

```
python2 grouping.py reactions results/ data/ eco -o tsr_eco.csv
```

This results in detecting metabolic conservation patterns (trail grouping by `reactions`) for species `eco`, using `CoMetGeNe` results stored in `results/` and metabolic pathway maps stored in KGML format in `data/` (see the *Note regarding directory structure* above). The output of trail grouping by reactions is stored in CSV format in the file `tsr_eco.csv`.

The CSV file contains a line for every reaction of the reference species *S* involved in `CoMetGeNe` trails of *S*. Groups of reactions involved in `CoMetGeNe` trails of *S* are separated by the line `***`. Note that a given reaction may appear several times in the CSV file, if it occurs in several `CoMetGeNe` trails of *S*. In this CSV file, the line for a reaction *r* in a `CoMetGeNe` trail of *S* contains:

* The KEGG R number for reaction *r*.
* The gene name(s) of the gene(s) of *S* involved in reaction *r*.
* The KEGG pathway map ID(s) for the pathway(s) in which the R number associated to *r* occurs.
* A column for every other species *S'* in the data set that can take one of the three following values:
  * A cross (`x`) if *r* is performed in species *S'* by the product of at least one gene neighboring at least one other gene involved in the `CoMetGeNe` trail to which reaction *r* belongs;
  * A dot (`.`) if *r* is performed in species *S'* by the product of a gene that is not a neighbor of at least one other gene involved in the `CoMetGeNe` trail to which reaction *r* belongs.
  * A circle (`o`) if *r* is absent from species *S'*.

**Example:** Suppose trail finding was performed for species `aae`, `bbn`, `eco`, and `mpn`. A small part of the CSV obtained when grouping `CoMetGeNe` trails by reactions for `eco` as the reference species is reproduced below (slightly re-formatted for readability purposes):

```
reaction;eco_gene;pathway;                      aae;bbn;mpn
R07618;  b0116;   00010 00020 00280 00620 00640; . ; o ; x 
R03270;  b0114;   00010 00020 00620;             o ; o ; x
R00014;  b0114;   00010 00020 00620;             o ; o ; x
R02569;  b0115;   00010 00020 00620;             o ; o ; x
```

From the table above, it can be seen that:

* Species `aae` performs only reaction `R07618` but none of the three other reactions;
* Species `bbn` performs none of the four reactions;
* Species `mpn` performs all four reactions using products of neighboring genes.

## Requirements

* Python 2.7
* cPickle (Python library)
* [lxml](http://lxml.de) (Python library)
* [NetworkX](https://networkx.github.io) (Python library)
* Active internet connection (necessary for automatic data retrieval from KEGG)
* Multi-core CPU (recommended)

## License 

The `CoMetGeNe` software is available under the MIT license.
