#  Gapfilling thoughts and comments

New thoughts and comments on gapfilling, our new big problem!

# A gapfilling model

First, we need a model to test whether gapfilling works.

The model models/pickles/Citrobacter.67826.8.p (note this is a binary pickle file) is a raw model built from the
functional roles listed in
[models/Citrobacter/263199/67826.8.assigned_functions](models/Citrobacter/263199/67826.8.assigned_functions). These
roles are directly from the genome annotation in RAST, and of course do not result in growth.

It is relatively easy to read this code and run a model.

First, we need to create a pickle object that we can reuse time and again. (We don't include these in the github repo
because they are binary files and dependent on your python installation.) You can easily create the pickle file like
this. This code will also run the fba for you, so if this works you should have everything installed correctly!

``` 
python  flux_balance_analysis/pickle_fba.py -a models/Citrobacter/263199/67826.8.assigned_functions -g 67826.8 -m media/ArgonneLB.txt -p Citrobacter.67826.8.p 
```

Now we can use that file to load the reactions and enzyme data, and start playing with it:

```python

import os
import sys
import flux_balance_analysis as fba

data = fba.unpickle_fba("models/pickles/Citrobacter.67826.8.p")
# data has the following keys: 
# media : The media object consisting of compound.Compunds 
# enzymes_and_reactions : All the enzymes and reactions in our database 
# roles : The roles we identified from the assigned_functions file 
# reactions_to_run : The reactions we are trying to run 
# drainflux : The drainflux reaction object 
# biomass : The biomass equation object 
# status : The status of the FBA when it was pickled 
# value : The value of the FBA when it was pickled 
# fba_sm : The stoichiometic matrix as a hash 
# fba_cp : The compounds in the sm as an array 
# fba_rc : The reactions in the sm as an array

media = data['media']
enz_r = data['enzymes_and_reactions']
reactions2run = data['reactions_to_run'] 
df = data['drainflux'] 
bme = data['biomass']

status, value, growth = fba.run_fba(enz_r, reactions2run, media, bme, fba)

print("Status: " + str(status) + " Value: " + str(value) + " Growth: " + str(growth))

``` 

Note that one of the key items here is the object _reactions2run_. This is a python
[set](https://docs.python.org/2/tutorial/datastructures.html#sets) which just has reaction IDs (typically in the format
rxn12345) from
[Biochemistry/ModelSEEDDatabase/Biochemistry/reactions.master.tsv](../Biochemistry/ModelSEEDDatabase/Biochemistry/reacti
ons.master.tsv). **To change the reactions that are included in the model, just change the set of reaction IDs in
reactions2run.** (Note: you may need to make a copy of reactions2run using copy.copy if you want to make multiple
changes and see their effects).


When I run this I get the output

``` 
Status: opt Value: 7.11542086152e-14 Growth: False 
```

The status is from the linear program solver and is suggesting this was an optimal solution. The value is the flux which
in this case is not enough for growth, and growth is a boolean about whether this model grew or not. Basically growth is
True if the value is >1, and it is largely a convenience feature to allow you to ask

```python
if growth: 
    something 
else:
    something else 
```

The goal of gapfilling is to add the minimum set of reactions to get growth to be True _under the right conditions_ and
to be False when appropriate too!

## Gapfilling parameters

We want to get growth from the gapfilling, and we have to choose the right set of reactions. Our [reactions
list](../Biochemistry/ModelSEEDDatabase/Biochemistry/reactions.master.tsv) has 34,696 reactions listed in it, and we
need to choose the appropriate ones!

 
## Ways that we can do gapfilling!

Here are some of the ways that we have tried gapfilling this model and some of the problems and issues that we run into!

### Essential reactions

As we describe [elsewhere](../NOTES.md#sbml-files) there are a set of essential reactions that every model has. There
are 110 reactions that are in every model (one of which is biomass), so we may as well start with those and make sure we
have them.

Note that of these reactions, 18 do not have any proteins associated with them (spontaneous reactions).

### Predefined reactions

Note that the Citrobacter requires (at least one of) six predefined reactions to ensure growth. These are additional
reactions that do not have proteins associated with them but are essential to get growth. They are:

Reaction ID | Name | Equation
--- | --- | ---
rxn13477 | (2E,6E)-farnesyl-diphosphate:isopentenyl-diphosphate farnesyltranstransferase (adding 3 isopentenyl units) |  (3) Isopentenyldiphosphate[0] + (1) Farnesyldiphosphate[0] <=> (3) PPi[0] + (3) H+[0] + (1) all-trans-Hexaprenyl  diphosphate[0] 
rxn09345 | Undecaprenyl diphosphate synthase | (8) Isopentenyldiphosphate[0] + (1) Farnesyldiphosphate[0] <=> (8) PPi[0] + (8) H+[0] + (1) Bactoprenyl diphosphate[0]
rxn01207 | 4-methyl-2-oxopentanoate:NAD+ oxidoreductase (CoA-mehtylpropanoylating) | (1) NAD[0] + (1) CoA[0] + (1) 4MOP[0] <=> (1) NADH[0] + (1) CO2[0] + (1) Isovaleryl-CoA[0] 
rxn05039 | R07280 | (1) H2O[0] + (1) 5-Amino-6--5-phosphoribitylaminouracil[0] <=> (1) Phosphate[0] + (1) 4--1-D-Ribitylamino-5-aminouracil[0] 
rxn09225 | rhamnosyltransferase I (LPS core biosynthesis) | (4) H+[0] + (1) dTDP-rhamnose[0] + (1) kdo-phospho-heptosyl-phospho-heptosyl-heptosyl-kdo2-lipidA[0] <=> (1) dTDP[0] + (1) inner core oligosaccharide lipid A[0] 
rxn05255 | folate transport via proton simport | (1) H+[1] + (1) Folate[1] <=> (1) H+[0] + (1) Folate[0]

By adding the `-t` option to the `flux_balance_scripts/gap_fill_multiple_reactions.py` code, we can test the
essentiality of these six reactions to see which one needs to be added. (Note, by default this will test _all_
reactions, but you can change the line `for r in reactions2run:` at about line 571 to be a specific set (e.g. `for r in
predef_reactions:`) to limit to just that set.)


reaction | essential? 
--- | --- 
rxn13477 | not essential
rxn05255 | not essential 
rxn01207 | not essential 
rxn09345 | not essential 
rxn05039 | not essential 
rxn09225 | essential

The only reaction we need to add is

Reaction ID | Name | Equation 
--- | --- | --- 
rxn09225 | rhamnosyltransferase I (LPS core biosynthesis) | (4) H+[0] + (1) dTDP-rhamnose[0] + (1) kdo-phospho-heptosyl-phospho-heptosyl-heptosyl-kdo2-lipidA[0] <=> (1) dTDP[0] + (1) inner core oligosaccharide lipid A[0]

Note that this reaction does not have any proteins associated with it.

Lets take a look at the compounds in this reaction:

Compound | Reactions in this model (before we add rxn09225) | All reactions
--- | --- | ---
H+ (location: c) | 757 | 17160
kdo-phospho-heptosyl-phospho-heptosyl-heptosyl-kdo2-lipidA (location: c) | 1 | 2
dTDP-rhamnose (location: c) | 1 | 31 
dTDP (location: c) | 4 | 110 
inner core oligosaccharide lipid A (location: c) | 1 | 2

The key here is there are several compounds that are only in one reaction (i.e. they do not go anywhere!)

### Reactions that are associated with orphan compounds

Some of our compounds only appear once or twice in the reactions, and so we need to add reactions that would integrate
those to the network. However, we need to ignore external compounds (unless they need to be imported - see below).
Therefore, we only consider intracellular compounds that are only present in a few reactions, and add the rest of the
reactions that they are associated with. It would be great if we used this as one piece of evidence in our probability
approach.

Max number of reactions | Internal compounds | Reactions assoc. with internal compounds | External compounds | Reactions assoc. with external compounds
--- | --- | --- | --- | --- 
1 | 0 | 0 | 0 | 0 
2 | 370 | 3322 | 124 | 16711 
3 | 765 | 6938 | 141 | 17326 
4 | 879 | 8471 | 146 | 17497 
5 | 953 | 10084 | 147 | 17563 
6 | 994 | 11407 | 148 | 17603 
7 | 1012 | 11839 | 149 | 18024 
8 | 1025 | 12469 | 149 | 18024 
9 | 1040 | 12891 | 149 | 18024 
10 | 1048 | 13517 | 149 | 18024

In other words, there are 370 intracellular compounds that are in a single reaction in our network (0 < reactions <2),
and those compounds are in an additional 3,322 reactions that are not currently in our network. Thus, we should consider
adding them.

To add these use
```python
orphan_cpd_reactions = gf.suggest_by_compound(enz_r, reactions2run, 2)
```


### Reactions that will import everything from the outside that we don't have

We need to make sure that all external compounds _in the media_ are involved in at least one import reaction to make
sure we can get it in. I am assuming that transport is not a limitation to growth, which maybe an error, but I am not
sure.

In our base model these compounds are not imported:

    * Cd2+ (location: e) is not imported 
    * chromate (location: e) is not imported 
    * Hg2+ (location: e) is not imported
    * Molybdate (location: e) is not imported
    * Ni2+ (location: e) is not imported
    * Cu2+ (location: e) is not imported
    * L-Cystine (location: e) is not imported

In our test, we limit our reactions to (1) import reactions and (2) reactions that have the external compound on the
left hand side (we don't test whether the compound itself is actually imported!)

The approach proposes an additional 19 reactions that will import all of these compounds. We should probably more
selective about those reactions.

### Reactions from subsystems

Our reaction list that we have generated from the roles in *assigned_functions* connects us to the subsystems that
should be present in the genome. We can figure out which other roles should be there based on incomplete inclusion of
subsystems. We *should* generate this by taking into account things like scenarios and variant codes to ensure that we
have active subsystems, but at the moment we are just adding all roles that are also present in the subsystem and
getting the connection to that. (Note: for convenience,
[flux_balance_analysis](../flux_balance_analysis/roles_and_reactions.py) has two helper functions: *roles_to_reactions*
and *reactions_to_roles*. These both accept a set of reactions or roles and convert to a hash of the other. In the hash,
the key is the reaction (or role) and the value is a set of roles (or reactions) that connects to).


### Reactions from close genomes

RAST generates a list of approximately 30 close genomes (for example, see
[closest.genomes](models/Citrobacter/263199/closest.genomes) for the _Citrobacter_ version), and we can use the [SEED
Servers](SEED_Servers_Python) to pull the list of roles in those genomes. Then we can convert those to reactions, and
add any that have been missed. We generate a [list of the functional roles and probability that they should be included](../models/Citrobacter/263199/closest.genomes.reactions)

### Reactions from other members of the same genera

We can also create a list of all roles in nearby genera and add that to our model, in the same way.

To see the code that made both of these lists, check out [Initial thoughts on
gapfilling](https://github.com/linsalrob/FuzzyMetabolicNetworks/blob/master/gapfilling/README.md#gapfilling-th
e-models)

### Limiting by compound presence 

At each of these steps, we can limit our choice of reactions to only those that  have compounds in our network. Note
that as we have [seen before](https://github.com/linsalrob/FuzzyMetabolicNetworks/blob/master/NOTES.md#how-many-reactions-is-each-compound-involved-in)
there are some compounds that are in a lot of reactions (e.g. H2O, ATP), and so we need to limit the maximum number 
of reactions a compound can be in for us to consider including it.


Max reactions per compound: | 10 | 20 | 30 | 40 | 50 | 60 | 70 | 80 | 90 | 100 |
 --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
 essential reactions | 42 | 42 | 42 | 43 | 44 | 44 | 44 | 44 | 49 | 49 
 subsystem reactions | 607 | 640 | 649 | 663 | 672 | 674 | 674 | 679 | 696 | 723 
 close reactions | 15 | 16 | 17 | 17 | 17 | 17 | 17 | 21 | 21 | 21 
genus reactions | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 | 2 
reactions with proteins | 148 | 166 | 193 | 208 | 208 | 214 | 218 | 218 | 254 | 262 
predefined reactions | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6 | 6

![Reactions suggested with links to compounds](../images/cpd_rcts.png "There seems to be an inflexion between 80 and 90
compounds per reaction")

### Reactions with other proteins

As a nearly last resort we probably want to add in all the reactions that have proteins associated with them that we
have not included. I don't really like this as it throws everything in together, but it is one way to make the model
grow.

### Reactions without proteins

This is a really bad idea, but you could add all 30,000 reactions that don't have proteins! This should ensure growth,
but good luck doing the computation and/or figuring out what it means.

### Reactions based on probability

We can calculate the probability that a reaction should be added based on the proportion of compounds on the left side
of the equation. Presumably, if there are a lot of compounds on the left side of the equation, there is a greater 
likelihood that reaction should be included.

#### Code example

	All of these approaches are put together in a single example piece of code `gap_fill_multiple_reactions.py` which will
	do each of the above (but not limit things by proteins) and try to run the reaction. You can run that like this:

``` 
python  flux_balance_scripts/gap_fill_multiple_reactions.py  -p models/pickles/Citrobacter.67826.8.p -g
models/Citrobacter/263199/citrobacter.reactions -c models/Citrobacter/263199/closest.genomes.reactions

```

### Model construction summary

The current steps to build and gap fill the model can be summarized as:


1. Original reactions: the reactions from the genome 
2. Essential reactions: other reactions that are required 
3. Predefined reactions: one or two reactions that must be included 
4. Orphan reactions: reactions that only connect to a single compound 
5. Subsystems: reactions that complete subsystems found in the genome 
6. Media reactions: reactions that ensure import of media components 
7. Close genomes: reactions from close genomes 
8. Other members of this genera: reactions from other members of this genera 
9. With proteins: other reactions that have been assigned to proteins

### Unique reactions

With these 9 steps, how many of the reactions are unique and how many are the same?

The table shows the number of unique reactions per approach, and the similarity between any two approaches:

Approach | Original | Essential | Predefined | Media | Subsystems | Close Genomes | Other Genera | Orphans | With proteins | Unique
 --- | ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---:| ---: 
Original | **1304** | 56 | 0 | 0 | 1270 | 1228 | 1283 | 0 | 1304 | 0 
Essential | 56 | **109** | 0 | 0 | 87 | 55 | 74 | 37 | 91 | 13 
Predefined | 0 | 0 | **1** | 0 | 0 | 0 | 0 | 1 | 0 | 0 
Media | 0 | 0 | 0 | **19** | 3 | 1 | 1 | 0 | 3 | 16 
Subsystems | 1270 | 87 | 0 | 3 | **2084** | 1500 | 1451 | 309 | 2084 | 0 
Close Genomes | 1228 | 55 | 0 | 1 | 1500 | **1558** | 1362 | 104 | 1558 | 0
Other Genera | 1283 | 74 | 0 | 1 | 1451 | 1362 | **1495** | 81 | 1495 | 0 
Orphans | 0 | 37 | 1 | 0 | 309 | 104 | 81 | **3489** | 374 | 3109 
With proteins | 1304 | 91 | 0 | 3 | 2084 | 1558 | 1495 | 374 | **2539** | 336

Each number is the number of reactions two approaches to suggesting reactions for the model have in common. Thus the
original reactions and the essential reactions have 56 reactions in common.

The number in bold (the diagonal) is the number of reations an approach has in common with itself, thus there are 1,304
reactions in the initial prediction, and 109 reactions in the essential reactions.

The last column is the number of unique reactions in that approach. Mostly they are zero, which means that all the
reactions would be suggested by at least one other approach. But for example, there are 13 essential reactions that are
not suggested by any other approach.


### Permutations

The code flux_balance_scripts/permutate_gapfilling.py` runs all possible permutations of these nine steps and figures out the
minimum number of reactions that we would need to add to create growth of our model:

``` 
python flux_balance_scripts/permutate_gapfilling.py -p models/pickles/Citrobacter.67826.8.p -g
models/Citrobacter/263199/citrobacter.reactions -c models/Citrobacter/263199/closest.genomes.reactions 
```

**Note**: This code runs the FBA in parallel, and will likely freeze up your computer!

### Figuring out required reactions (bisecting reactions)

So what are these reactions that we need to add?

This is derived from comments in [flux_balance_scripts](../flux_balance_scripts) that I wrote on bisection and tested
with the *Acidobacter* model I have been playing with. I amended it here for the *Citrobacter* model that I am playing
with.

# Bisecting out what reactions are interesting

The approach that we take is to create our list of suggested reactions, and then try to figure out the minimal set that
complete our growth.

We basically take the list and split it in half, and then test each half. The rules that I decided upon are:

* If both halves grow we choose one (the left one) and carry on with that. (Note, we should probably store both halves
* as potential solutions). If the left half grows we carry on with that until we get to a single reaction and we call
* that good If the right half grows we carry on with that until we get to a single reaction and we call that good If
* neither half grows we rejoin both halves and shuffle them so that we put the reactions in a different order and try
* again. We allow upto 5 tries


## A wrangle on bisection

The 50:50 split on the bisections is fast (O log time complexity) but it is a bit *heavy handed* for finding the
appropriate reactions - often there are lots of reactions that are needed. Therefore, I wrote a variation on the
bisection code that makes an uneven split and tests both parts. I start with a 40%:60% split, and then iterate through
dividing the lower half by two (20%, 10%, 5%, 2.5% etc) until we get to a list with 1 element. At each split we test
both halves and if one fraction grows (growth returns true) we discard the other fraction and continue. If neither
fraction results in growth we shuffle the order of the list and then start again at 50%. We do this upto 5 times before
giving up and moving on.

This allows us to refine the additional reactions that are necessary for growth of the model in addition to those that
we have from the genome sequence.

The bisection run on the *Citrobacter* test case was pretty successful:


Approach | base reactions | suggested reactions | essential reactions
Base model that doesn't grow | 1304 | | 
Beginning model that grows | 2540 | |
With proteins | 2199 | 341 | 2
probability | 1706 | 495 | 53
genera | 1690 | 69 | 29
close genomes | 1389 | 330 |24
predefined | 1412 | 1 | 1
Final model that grows | 1413 | |

2540 is 1304 + 1 + 330 + 69 + 495 + 341
2199 is 1304 + 1 + 330 + 69 + 495 
1706 is 1304 + 1 + 330 + 69 + 2
1690 is 1304 + 1 + 330 + 2 + 53
1389 is 1304 + 1 + 2 + 53 + 29


We went from 2540 reactions in the beginning gapfilled model and ended up with 1413 reactions in the final model. 
We only added 109 reactions in gapfilling.

### Code example

You can run the bisections code like this. Note that the `gap_fill_bisections.py` has a helper method that handles
multiple bisections for us, although this is not (currently) threaded.

```
python flux_balance_scripts/gap_fill_bisections.py  -p Citrobacter.67826.8.p -g models/Citrobacter/263199/citrobacter.reactions -c models/Citrobacter/263199/closest.genomes.reactions > temp/gf_bisections.out 2> temp/gf_bisections.err
```



## Likelihood based gapfilling

The approach that Benedict *et al.,* took [Likelihood-Based Gene Annotations for Gap Filling and Quality Assessment in
Genome-Scale Metabolic Models](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003882) has four 
distinct steps to calculate the likelihood that a reaction should be included in a model. This is described in the 
section *'Calculating reaction likelihoods'* on pg. 11 of their paper:

1. Conversion of gene annotations to reactions.
2. Linking other genes to annotations by blasting against their database of all members of a reaction and storing any
gene that had >=80% of the maximum score for that organism
3. Computation of the presence of complexes based on the presence of functional roles
4. Computation of the presence of reactions based on the presence of complexes

We actually assume the last two steps already in our model ... although perhaps we should not! When the 
*assigned_functions* file is read we use [fba.roles_to_reactions](../fba/roles_and_reactions.py) to convert from those
functional roles to a set of reaction IDs. There is an assumption in that file at the moment that if any member of the
complex is present, the whole complex is present, and thus the reaction is present. We may need to amend that 
so that the reaction is only present if some proportion of the complex is present.

In that paper, they do not include a discussion of the addition of reactions for which there is no known protein 
involved!


# Running gapfilling on minimal media models

I ran the *Citrobacter* gapfilling on each of the minimal media conditions for which the *Citrobacter* model was 
supposed to grow. In the file [Citrobacter/Citrobacter_sedlakii_growth](../models/Citrobacter/Citrobacter_sedlakii_growth)
 we have each of the 96 media conditions and whether or not *C. sedlakii* grew in those conditions. I took all instances 
where we had growth and built an independent gapfilled model:

```bash
for m in $(grep True models/Citrobacter/Citrobacter_sedlakii_growth  | cut -f 1); do 
python flux_balance_scripts/compare_minimal_media.py  -m media/$m.txt -p Citrobacter.67826.8.p \
 > temp/$m.gf.txt 2> temp/$m.gf.err & echo $m; 
done
```

This took ~9-12 hours per condition (but I ran them all in parallel on rambox), and generated a model that grows. Each
model was different (see below).

I summarized the reactions added to each model:

```bash
(echo "{"; 
grep -h 'To gapfill Citrobacter.67826.8.p for growth on' temp/*err | sed -e \
's/To gapfill Citrobacter.67826.8.p for growth on /"/; s/ we added /" : /; s/$/,/' ; 
echo "}") > ~/public_html/gf.txt
```

and then using python extracted interesting data:

### Number of gapfilled reactions


```python
c=[]
with open('temp/gf.txt', 'r') as f:
    for l in f:
        c.append(l.strip())
h=eval(''.join(c))
for k in sorted(h.keys(), key=lambda k: len(h[k]), reverse=True):
    sys.stderr.write(k.replace('media/', '').replace('.txt', '') + " | " + str(len(h[k])) + "\n")
```

Media conditions | Number of reactions added for growth
---|---
MOPS_NoC_D-Fructose | 48
MOPS_NoC_L-Arabinose | 53
MOPS_NoC_D-Glucosamine | 54
MOPS_NoN_L-Arginine | 54
MOPS_NoC_D-Mannose | 58
MOPS_NoC_Succinate | 58
MOPS_NoC_D-Galactose | 59
MOPS_NoN_Cytosine | 60
MOPS_NoC_L-Rhamnose | 61
MOPS_NoC_Alpha-D-Glucose | 61
MOPS_NoC_Inosine | 62
MOPS_NoC_Dulcitol | 62
MOPS_NoN_Cytidine | 62
MOPS_NoC_Glycerol | 62
MOPS_NoC_L-Fucose | 64
MOPS_NoN_N-Acetyl-D-Glucosamine | 64
MOPS_NoC_Adenosine | 64
MOPS_NoN_Adenine | 66
MOPS_NoC_D-Xylose | 67
MOPS_NoC_Lactate | 67
MOPS_NoC_Melibiose | 68
MOPS_NoC_D-Serine | 70
MOPS_NoC_L-Asparagine | 71
MOPS_NoC_D-Arabinose | 71
MOPS_NoC_Pyruvate | 72
MOPS_NoC_L-Glutamine | 79
MOPS_NoN_Adenosine | 80
MOPS_NoC_D-Glucose | 98
MOPS_NoC_D-Ribose | 106
MOPS_NoN_L-Proline | 107
MOPS_NoC_Trehalose | 108
MOPS_NoC_Cellobiose | 120


I think part of this is the stochasticity of the gapfilling approach, especially with the finding additional reactions

### Similarity between solutions

We can see how similar each reaction is to any other:

```python
def jaccard(s1, s2):
    """
    Calculate the Jaccard distance between two sets
    """

    if len(s1) == 0 or len(s2) == 0:
        return 1
    return 1-1.0*len(s1.intersection(s2))/len(s1.union(s2))

ks = h.keys()
ks.sort()
d = {x:{} for x in ks}
for k in h:
    for l in h:
        j=jaccard(h[k], h[l])
        d[k][l] = d[l][k] = j
print(" | " + " | ".join([x.replace('media/', '').replace('.txt', '') for x in ks]))
print(" --- |" + " | ".join([" --- " for x in ks]))
for k in ks:
    sys.stdout.write(k.replace('media/', '').replace('.txt', ''))
    for l in ks:
        sys.stdout.write(" | %0.2f" % d[k][l])
    sys.stdout.write("\n")
```


Media | MOPS_NoC_Adenosine | MOPS_NoC_Alpha-D-Glucose | MOPS_NoC_Cellobiose | MOPS_NoC_D-Arabinose | MOPS_NoC_D-Fructose | MOPS_NoC_D-Galactose | MOPS_NoC_D-Glucosamine | MOPS_NoC_D-Glucose | MOPS_NoC_D-Mannose | MOPS_NoC_D-Ribose | MOPS_NoC_D-Serine | MOPS_NoC_D-Xylose | MOPS_NoC_Dulcitol | MOPS_NoC_Glycerol | MOPS_NoC_Inosine | MOPS_NoC_L-Arabinose | MOPS_NoC_L-Asparagine | MOPS_NoC_L-Fucose | MOPS_NoC_L-Glutamine | MOPS_NoC_L-Rhamnose | MOPS_NoC_Lactate | MOPS_NoC_Melibiose | MOPS_NoC_Pyruvate | MOPS_NoC_Succinate | MOPS_NoC_Trehalose | MOPS_NoN_Adenine | MOPS_NoN_Adenosine | MOPS_NoN_Cytidine | MOPS_NoN_Cytosine | MOPS_NoN_L-Arginine | MOPS_NoN_L-Proline | MOPS_NoN_N-Acetyl-D-Glucosamine
 --- | ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---  |  ---
MOPS_NoC_Adenosine | 0.00 | 0.74 | 0.85 | 0.71 | 0.77 | 0.72 | 0.74 | 0.89 | 0.80 | 0.86 | 0.81 | 0.81 | 0.80 | 0.81 | 0.75 | 0.77 | 0.84 | 0.72 | 0.82 | 0.80 | 0.72 | 0.73 | 0.70 | 0.82 | 0.84 | 0.74 | 0.75 | 0.80 | 0.71 | 0.76 | 0.83 | 0.81
MOPS_NoC_Alpha-D-Glucose | 0.74 | 0.00 | 0.83 | 0.80 | 0.80 | 0.71 | 0.74 | 0.85 | 0.75 | 0.87 | 0.88 | 0.83 | 0.82 | 0.77 | 0.78 | 0.73 | 0.78 | 0.79 | 0.74 | 0.83 | 0.79 | 0.77 | 0.79 | 0.84 | 0.84 | 0.80 | 0.79 | 0.76 | 0.73 | 0.75 | 0.83 | 0.76
MOPS_NoC_Cellobiose | 0.85 | 0.83 | 0.00 | 0.89 | 0.88 | 0.86 | 0.77 | 0.90 | 0.85 | 0.90 | 0.89 | 0.85 | 0.86 | 0.85 | 0.86 | 0.84 | 0.88 | 0.86 | 0.84 | 0.86 | 0.87 | 0.88 | 0.87 | 0.85 | 0.90 | 0.85 | 0.87 | 0.83 | 0.82 | 0.82 | 0.86 | 0.84
MOPS_NoC_D-Arabinose | 0.71 | 0.80 | 0.89 | 0.00 | 0.73 | 0.80 | 0.75 | 0.83 | 0.81 | 0.85 | 0.82 | 0.79 | 0.76 | 0.81 | 0.77 | 0.83 | 0.83 | 0.73 | 0.83 | 0.80 | 0.73 | 0.74 | 0.76 | 0.76 | 0.87 | 0.75 | 0.78 | 0.80 | 0.73 | 0.79 | 0.83 | 0.79
MOPS_NoC_D-Fructose | 0.77 | 0.80 | 0.88 | 0.73 | 0.00 | 0.70 | 0.79 | 0.83 | 0.77 | 0.82 | 0.74 | 0.75 | 0.69 | 0.74 | 0.66 | 0.71 | 0.72 | 0.76 | 0.74 | 0.67 | 0.76 | 0.77 | 0.75 | 0.55 | 0.79 | 0.76 | 0.80 | 0.66 | 0.68 | 0.74 | 0.82 | 0.76
MOPS_NoC_D-Galactose | 0.72 | 0.71 | 0.86 | 0.80 | 0.70 | 0.00 | 0.81 | 0.85 | 0.77 | 0.85 | 0.78 | 0.83 | 0.75 | 0.77 | 0.71 | 0.70 | 0.77 | 0.78 | 0.80 | 0.74 | 0.78 | 0.75 | 0.76 | 0.76 | 0.78 | 0.82 | 0.76 | 0.75 | 0.75 | 0.65 | 0.88 | 0.77
MOPS_NoC_D-Glucosamine | 0.74 | 0.74 | 0.77 | 0.75 | 0.79 | 0.81 | 0.00 | 0.83 | 0.77 | 0.88 | 0.86 | 0.70 | 0.77 | 0.68 | 0.75 | 0.68 | 0.80 | 0.70 | 0.76 | 0.81 | 0.78 | 0.77 | 0.80 | 0.78 | 0.84 | 0.75 | 0.78 | 0.70 | 0.69 | 0.79 | 0.78 | 0.78
MOPS_NoC_D-Glucose | 0.89 | 0.85 | 0.90 | 0.83 | 0.83 | 0.85 | 0.83 | 0.00 | 0.82 | 0.87 | 0.79 | 0.87 | 0.80 | 0.83 | 0.84 | 0.79 | 0.79 | 0.86 | 0.84 | 0.84 | 0.89 | 0.84 | 0.88 | 0.86 | 0.86 | 0.88 | 0.82 | 0.80 | 0.86 | 0.80 | 0.85 | 0.79
MOPS_NoC_D-Mannose | 0.80 | 0.75 | 0.85 | 0.81 | 0.77 | 0.77 | 0.77 | 0.82 | 0.00 | 0.82 | 0.78 | 0.82 | 0.79 | 0.74 | 0.80 | 0.74 | 0.75 | 0.80 | 0.83 | 0.79 | 0.76 | 0.81 | 0.80 | 0.77 | 0.87 | 0.81 | 0.85 | 0.78 | 0.84 | 0.76 | 0.85 | 0.73
MOPS_NoC_D-Ribose | 0.86 | 0.87 | 0.90 | 0.85 | 0.82 | 0.85 | 0.88 | 0.87 | 0.82 | 0.00 | 0.72 | 0.82 | 0.83 | 0.84 | 0.86 | 0.86 | 0.77 | 0.83 | 0.89 | 0.83 | 0.84 | 0.86 | 0.85 | 0.82 | 0.88 | 0.81 | 0.82 | 0.89 | 0.87 | 0.81 | 0.89 | 0.84
MOPS_NoC_D-Serine | 0.81 | 0.88 | 0.89 | 0.82 | 0.74 | 0.78 | 0.86 | 0.79 | 0.78 | 0.72 | 0.00 | 0.71 | 0.74 | 0.80 | 0.84 | 0.81 | 0.76 | 0.80 | 0.87 | 0.78 | 0.82 | 0.82 | 0.82 | 0.75 | 0.82 | 0.80 | 0.78 | 0.85 | 0.86 | 0.71 | 0.90 | 0.78
MOPS_NoC_D-Xylose | 0.81 | 0.83 | 0.85 | 0.79 | 0.75 | 0.83 | 0.70 | 0.87 | 0.82 | 0.82 | 0.71 | 0.00 | 0.74 | 0.81 | 0.81 | 0.83 | 0.84 | 0.72 | 0.85 | 0.81 | 0.78 | 0.79 | 0.80 | 0.72 | 0.82 | 0.70 | 0.84 | 0.85 | 0.77 | 0.73 | 0.88 | 0.80
MOPS_NoC_Dulcitol | 0.80 | 0.82 | 0.86 | 0.76 | 0.69 | 0.75 | 0.77 | 0.80 | 0.79 | 0.83 | 0.74 | 0.74 | 0.00 | 0.80 | 0.76 | 0.80 | 0.78 | 0.67 | 0.83 | 0.73 | 0.79 | 0.75 | 0.78 | 0.75 | 0.84 | 0.80 | 0.83 | 0.81 | 0.78 | 0.73 | 0.88 | 0.73
MOPS_NoC_Glycerol | 0.81 | 0.77 | 0.85 | 0.81 | 0.74 | 0.77 | 0.68 | 0.83 | 0.74 | 0.84 | 0.80 | 0.81 | 0.80 | 0.00 | 0.75 | 0.63 | 0.75 | 0.81 | 0.84 | 0.79 | 0.83 | 0.85 | 0.83 | 0.78 | 0.84 | 0.79 | 0.83 | 0.68 | 0.77 | 0.74 | 0.83 | 0.76
MOPS_NoC_Inosine | 0.75 | 0.78 | 0.86 | 0.77 | 0.66 | 0.71 | 0.75 | 0.84 | 0.80 | 0.86 | 0.84 | 0.81 | 0.76 | 0.75 | 0.00 | 0.78 | 0.80 | 0.80 | 0.81 | 0.72 | 0.81 | 0.77 | 0.81 | 0.70 | 0.79 | 0.78 | 0.83 | 0.69 | 0.67 | 0.70 | 0.84 | 0.79
MOPS_NoC_L-Arabinose | 0.77 | 0.73 | 0.84 | 0.83 | 0.71 | 0.70 | 0.68 | 0.79 | 0.74 | 0.86 | 0.81 | 0.83 | 0.80 | 0.63 | 0.78 | 0.00 | 0.77 | 0.81 | 0.80 | 0.76 | 0.81 | 0.84 | 0.81 | 0.75 | 0.82 | 0.80 | 0.81 | 0.68 | 0.76 | 0.70 | 0.82 | 0.76
MOPS_NoC_L-Asparagine | 0.84 | 0.78 | 0.88 | 0.83 | 0.72 | 0.77 | 0.80 | 0.79 | 0.75 | 0.77 | 0.76 | 0.84 | 0.78 | 0.75 | 0.80 | 0.77 | 0.00 | 0.83 | 0.78 | 0.80 | 0.80 | 0.85 | 0.83 | 0.78 | 0.86 | 0.84 | 0.78 | 0.78 | 0.82 | 0.83 | 0.84 | 0.74
MOPS_NoC_L-Fucose | 0.72 | 0.79 | 0.86 | 0.73 | 0.76 | 0.78 | 0.70 | 0.86 | 0.80 | 0.83 | 0.80 | 0.72 | 0.67 | 0.81 | 0.80 | 0.81 | 0.83 | 0.00 | 0.79 | 0.74 | 0.73 | 0.69 | 0.68 | 0.80 | 0.85 | 0.75 | 0.76 | 0.82 | 0.73 | 0.81 | 0.84 | 0.80
MOPS_NoC_L-Glutamine | 0.82 | 0.74 | 0.84 | 0.83 | 0.74 | 0.80 | 0.76 | 0.84 | 0.83 | 0.89 | 0.87 | 0.85 | 0.83 | 0.84 | 0.81 | 0.80 | 0.78 | 0.79 | 0.00 | 0.78 | 0.81 | 0.79 | 0.78 | 0.79 | 0.85 | 0.83 | 0.82 | 0.77 | 0.77 | 0.85 | 0.79 | 0.79
MOPS_NoC_L-Rhamnose | 0.80 | 0.83 | 0.86 | 0.80 | 0.67 | 0.74 | 0.81 | 0.84 | 0.79 | 0.83 | 0.78 | 0.81 | 0.73 | 0.79 | 0.72 | 0.76 | 0.80 | 0.74 | 0.78 | 0.00 | 0.81 | 0.77 | 0.77 | 0.66 | 0.83 | 0.74 | 0.82 | 0.74 | 0.71 | 0.75 | 0.83 | 0.80
MOPS_NoC_Lactate | 0.72 | 0.79 | 0.87 | 0.73 | 0.76 | 0.78 | 0.78 | 0.89 | 0.76 | 0.84 | 0.82 | 0.78 | 0.79 | 0.83 | 0.81 | 0.81 | 0.80 | 0.73 | 0.81 | 0.81 | 0.00 | 0.74 | 0.71 | 0.80 | 0.86 | 0.76 | 0.78 | 0.84 | 0.77 | 0.80 | 0.84 | 0.80
MOPS_NoC_Melibiose | 0.73 | 0.77 | 0.88 | 0.74 | 0.77 | 0.75 | 0.77 | 0.84 | 0.81 | 0.86 | 0.82 | 0.79 | 0.75 | 0.85 | 0.77 | 0.84 | 0.85 | 0.69 | 0.79 | 0.77 | 0.74 | 0.00 | 0.70 | 0.80 | 0.84 | 0.78 | 0.78 | 0.80 | 0.77 | 0.79 | 0.83 | 0.80
MOPS_NoC_Pyruvate | 0.70 | 0.79 | 0.87 | 0.76 | 0.75 | 0.76 | 0.80 | 0.88 | 0.80 | 0.85 | 0.82 | 0.80 | 0.78 | 0.83 | 0.81 | 0.81 | 0.83 | 0.68 | 0.78 | 0.77 | 0.71 | 0.70 | 0.00 | 0.79 | 0.83 | 0.76 | 0.77 | 0.84 | 0.75 | 0.82 | 0.85 | 0.83
MOPS_NoC_Succinate | 0.82 | 0.84 | 0.85 | 0.76 | 0.55 | 0.76 | 0.78 | 0.86 | 0.77 | 0.82 | 0.75 | 0.72 | 0.75 | 0.78 | 0.70 | 0.75 | 0.78 | 0.80 | 0.79 | 0.66 | 0.80 | 0.80 | 0.79 | 0.00 | 0.81 | 0.71 | 0.84 | 0.75 | 0.73 | 0.70 | 0.86 | 0.79
MOPS_NoC_Trehalose | 0.84 | 0.84 | 0.90 | 0.87 | 0.79 | 0.78 | 0.84 | 0.86 | 0.87 | 0.88 | 0.82 | 0.82 | 0.84 | 0.84 | 0.79 | 0.82 | 0.86 | 0.85 | 0.85 | 0.83 | 0.86 | 0.84 | 0.83 | 0.81 | 0.00 | 0.83 | 0.85 | 0.81 | 0.83 | 0.76 | 0.91 | 0.85
MOPS_NoN_Adenine | 0.74 | 0.80 | 0.85 | 0.75 | 0.76 | 0.82 | 0.75 | 0.88 | 0.81 | 0.81 | 0.80 | 0.70 | 0.80 | 0.79 | 0.78 | 0.80 | 0.84 | 0.75 | 0.83 | 0.74 | 0.76 | 0.78 | 0.76 | 0.71 | 0.83 | 0.00 | 0.82 | 0.83 | 0.73 | 0.74 | 0.83 | 0.83
MOPS_NoN_Adenosine | 0.75 | 0.79 | 0.87 | 0.78 | 0.80 | 0.76 | 0.78 | 0.82 | 0.85 | 0.82 | 0.78 | 0.84 | 0.83 | 0.83 | 0.83 | 0.81 | 0.78 | 0.76 | 0.82 | 0.82 | 0.78 | 0.78 | 0.77 | 0.84 | 0.85 | 0.82 | 0.00 | 0.81 | 0.77 | 0.84 | 0.87 | 0.84
MOPS_NoN_Cytidine | 0.80 | 0.76 | 0.83 | 0.80 | 0.66 | 0.75 | 0.70 | 0.80 | 0.78 | 0.89 | 0.85 | 0.85 | 0.81 | 0.68 | 0.69 | 0.68 | 0.78 | 0.82 | 0.77 | 0.74 | 0.84 | 0.80 | 0.84 | 0.75 | 0.81 | 0.83 | 0.81 | 0.00 | 0.70 | 0.78 | 0.79 | 0.76
MOPS_NoN_Cytosine | 0.71 | 0.73 | 0.82 | 0.73 | 0.68 | 0.75 | 0.69 | 0.86 | 0.84 | 0.87 | 0.86 | 0.77 | 0.78 | 0.77 | 0.67 | 0.76 | 0.82 | 0.73 | 0.77 | 0.71 | 0.77 | 0.77 | 0.75 | 0.73 | 0.83 | 0.73 | 0.77 | 0.70 | 0.00 | 0.73 | 0.79 | 0.83
MOPS_NoN_L-Arginine | 0.76 | 0.75 | 0.82 | 0.79 | 0.74 | 0.65 | 0.79 | 0.80 | 0.76 | 0.81 | 0.71 | 0.73 | 0.73 | 0.74 | 0.70 | 0.70 | 0.83 | 0.81 | 0.85 | 0.75 | 0.80 | 0.79 | 0.82 | 0.70 | 0.76 | 0.74 | 0.84 | 0.78 | 0.73 | 0.00 | 0.87 | 0.74
MOPS_NoN_L-Proline | 0.83 | 0.83 | 0.86 | 0.83 | 0.82 | 0.88 | 0.78 | 0.85 | 0.85 | 0.89 | 0.90 | 0.88 | 0.88 | 0.83 | 0.84 | 0.82 | 0.84 | 0.84 | 0.79 | 0.83 | 0.84 | 0.83 | 0.85 | 0.86 | 0.91 | 0.83 | 0.87 | 0.79 | 0.79 | 0.87 | 0.00 | 0.83
MOPS_NoN_N-Acetyl-D-Glucosamine | 0.81 | 0.76 | 0.84 | 0.79 | 0.76 | 0.77 | 0.78 | 0.79 | 0.73 | 0.84 | 0.78 | 0.80 | 0.73 | 0.76 | 0.79 | 0.76 | 0.74 | 0.80 | 0.79 | 0.80 | 0.80 | 0.80 | 0.83 | 0.79 | 0.85 | 0.83 | 0.84 | 0.76 | 0.83 | 0.74 | 0.83 | 0.00

If we ignore ourself, the range of distances is generally from 0.7-0.9: i.e. not very similar.

### Reactions added in every gapfill:

There are 11 reactions that are added to every condition. Note that if we check for reactions in all but one model, all
but two models, or all but three models, we don't see any additional reactions, suggesting that there are no alternative
reactions we need to worry about (perhaps).

```python
count={}
for k in h:
    for s in h[k]:
        count[s] = count.get(s, 0)+1
for c in count:
    if count[c] == len(h.keys()):
        print(c)
```

Reaction | Name | Equation
--- | --- | ---
rxn00293 | UTP:N-acetyl-alpha-D-glucosamine-1-phosphate uridylyltransferase | (1) UTP[0] + (1) N-Acetyl-D-glucosamine1-phosphate[0] <=> (1) PPi[0] + (1) UDP-N-acetylglucosamine[0]
rxn03164 | UDP-N-acetylmuramoyl-L-alanyl-D-glutamyl-meso-2,6-diaminoheptanedioate:D-alanyl-D-alanine ligase(ADP-forming) | (1) ATP[0] + (1) Ala-Ala[0] + (1) UDP-N-acetylmuramoyl-L-alanyl-D-gamma-glutamyl-meso-2-6-diaminopimelate[0] <=> (1) ADP[0] + (1) Phosphate[0] + (1) H+[0] + (1) UDP-N-acetylmuramoyl-L-alanyl-D-glutamyl-6-carboxy-L-lysyl-D-alanyl- D-alanine[0]
rxn08333 | 1,4-dihydroxy-2-naphthoate octaprenyltransferase | (1) 1-4-Dihydroxy-2-naphthoate[0] + (1) Farnesylfarnesylgeraniol[0] <=> (1) CO2[0] + (1) PPi[0] + (1) 2-Demethylmenaquinol 8[0]
rxn00392 | ATP:riboflavin 5'-phosphotransferase | (1) ATP[0] + (1) Riboflavin[0] <=> (1) ADP[0] + (1) FMN[0] + (1) H+[0]
rxn03638 | Acetyl-CoA:D-glucosamine-1-phosphate N-acetyltransferase | (1) Acetyl-CoA[0] + (1) D-Glucosamine1-phosphate[0] <=> (1) CoA[0] + (1) H+[0] + (1) N-Acetyl-D-glucosamine1-phosphate[0]
rxn10094 | S-adenosylmethione:2-demethylmenaquinone methyltransferase | (1) S-Adenosyl-L-methionine[0] + (1) 2-Demethylmenaquinone 8[0] <=> (1) S-Adenosyl-homocysteine[0] + (1) H+[0] + (1) Menaquinone 8[0]
rxn03395 | S-adenosyl-L-methionine:3-(all-trans-octaprenyl)benzene-1,2-diol 2-O-methyltransferase | (1) S-Adenosyl-L-methionine[0] + (1) 2-Octaprenyl-6-hydroxyphenol[0] <=> (1) S-Adenosyl-homocysteine[0] + (1) H+[0] + (1) 2-Octaprenyl-6-methoxyphenol[0]
rxn03397 | UDP-L-rhamnose:flavonol-3-O-D-glucoside L-rhamnosyltransferase | (1) S-Adenosyl-L-methionine[0] + (1) 2-Octaprenyl-6-methoxy-1,4-benzoquinone[0] <=> (1) S-Adenosyl-homocysteine[0] + (1) H+[0] + (1) 2-Octaprenyl-3-methyl-6-methoxy-1,4-benzoquinone[0]
rxn11946 | R05614 | (1) S-Adenosyl-L-methionine[0] + (1) 2-Octaprenyl-3-methyl-5-hydroxy-6-methoxy-1,4-benzoquinone[0] <=> (1) S-Adenosyl-homocysteine[0] + (1) H+[0] + (1) Ubiquinone-8[0]
rxn02285 | UDP-N-acetylmuramate:NADP+ oxidoreductase | (1) NADP[0] + (1) UDP-MurNAc[0] <=> (1) NADPH[0] + (1) H+[0] + (1) UDP-N-acetylglucosamine enolpyruvate[0]
rxn00122 | ATP:FMN adenylyltransferase | (1) ATP[0] + (1) FMN[0] <=> (1) PPi[0] + (1) FAD[0]

These are the roles associated with each of these reactions. We need to check whether those are in the genome!


* rxn11946
  * 3-demethylubiquinol 3-O-methyltransferase (EC 2.1.1.64)
* rxn00293
  * N-acetylglucosamine-1-phosphate uridyltransferase eukaryotic (EC 2.7.7.23)
  * N-acetylglucosamine-1-phosphate uridyltransferase (EC 2.7.7.23)
  * Glucosamine-1-phosphate N-acetyltransferase (EC 2.3.1.157)
* rxn10094
  * Ubiquinone/menaquinone biosynthesis methyltransferase UbiE (EC 2.1.1.-)
* rxn03164
  * UDP-N-acetylmuramoylalanyl-D-glutamyl-2,6-diaminopimelate--D-alanyl-D-alanine ligase (EC 6.3.2.10)
* rxn00122
  * Riboflavin kinase (EC 2.7.1.26)
  * FMN adenylyltransferase, type 2 eukaryotic (EC 2.7.7.2)
  * FMN adenylyltransferase (EC 2.7.7.2)
* rxn00392
  * Riboflavin kinase (EC 2.7.1.26)
  * FMN adenylyltransferase (EC 2.7.7.2)
* rxn08333
  * 1,4-dihydroxy-2-naphthoate octaprenyltransferase (EC 2.5.1.74)
* rxn03397
  * Ubiquinone/menaquinone biosynthesis methyltransferase UbiE (EC 2.1.1.-)
* rxn03395
  * 3-demethylubiquinol 3-O-methyltransferase (EC 2.1.1.64)
* rxn03638
  * N-acetylglucosamine-1-phosphate uridyltransferase (EC 2.7.7.23)
  * Glucosamine-1-phosphate N-acetyltransferase (EC 2.3.1.157)
* rxn02285
  * UDP-N-acetylenolpyruvoylglucosamine reductase (EC 1.1.1.158)



To find the frequency of all the reactions added we can use:

```python
for c in sorted(count, key=lambda x: count[x], reverse=True):
    print(c + "| " + str(count[c]))
```    

Here are the additional *top* reactions that we needed to add:

reaction id | number of models (out of 32)
--- | ---
rxn00293 | 32
rxn03164 | 32
rxn08333 | 32
rxn00392 | 32
rxn03638 | 32
rxn10094 | 32
rxn03395 | 32
rxn03397 | 32
rxn11946 | 32
rxn02285 | 32
rxn00122 | 32
rxn03087 | 28
rxn08244 | 27
rxn01513 | 24
rxn05276 | 24
rxn05275 | 23
rxn07434 | 20
rxn07435 | 20
rxn32930 | 20
rxn28987 | 20
rxn06335 | 20


Now we have this set, we want to confirm the minimal set of reactions that grows on all of these conditions.

The most straightforward way to do this is to run a knock out experiment with each of these reactions, but there are 925
reactions in total, so that will take a while to run 925 reactions x 32 conditions!

Instead, I added the check for individual knockouts to the `compare_minimal_media.py` code and then reran that analysis 
(I also added a the 11 absolutely required reactions so we don't try and gapfill for those anymore). Also I only used 
20 media conditions so that there are four processors left idle on rambox!

# Finding essential reactions

In the second iteration I added those 11 reactions that were required in every condition, and tested 20 conditions. 

Media | Total Suggestions | Essential Suggestions
--- | ---: | ---:
MOPS_NoC_D-Glucose-6-phosphate | 69 | 33
MOPS_NoC_L-Asparagine | 64 | 30
MOPS_NoC_D-Ribose | 59 | 36
MOPS_NoC_Melibiose | 58 | 32
MOPS_NoC_Succinate | 56 | 38
MOPS_NoC_Inosine | 52 | 32
MOPS_NoC_L-Arabinose | 50 | 29
MOPS_NoC_D-Galactose | 45 | 31
MOPS_NoC_Dulcitol | 44 | 29
MOPS_NoC_Lactate | 42 | 33
MOPS_NoC_L-Fucose | 42 | 27
MOPS_NoC_D-Fructose | 42 | 32
MOPS_NoC_Glycerol | 40 | 33
MOPS_NoC_D-Mannose | 40 | 34
MOPS_NoC_D-Glucose | 39 | 32
MOPS_NoC_Trehalose | 39 | 30
MOPS_NoC_Pyruvate | 38 | 30
MOPS_NoC_Alpha-D-Glucose | 38 | 33
MOPS_NoC_Cellobiose | 37 | 27
MOPS_NoC_D-Serine | 37 | 31

In most cases, the majority of suggestions are essential, but there are still a few that we could trim away.

There are 129 different reactions that are required to complete growth on all of these media. Lets add these to
the existing model and see how they perform on other media.


We can now compare how these perform on different media.

We create a file called `essential_rxns.py` that has the additional 129 reactions, and then run:

```bash
python flux_balance_scripts/reactions_different_media.py -p Citrobacter.67826.8.p -m media/ -r temp/essential_rxns.py -v
```

We combined these reactions with our base of 1315 reactions for a total model with 1444 reactions, 
and we see (this is just a selection of media, but everything has growth):

Media | FBA Value | Growth
--- | --- | ---
MOPS_NoC_Alpha-D-Lactose | 100.0 | True
MOPS_NoN_Tyramine | 100.0 | True
MOPS_NoC_Oxalate | 100.0 | True
MOPS_NoC_D-Glucose | 100.0 | True
MOPS_NoC_Thymidine | 100.0 | True
MOPS_NoC_L-Cysteate | 100.0 | True
MOPS_NoC_L-Aspartate | 100.0 | True
MOPS_NoC_D-Glucosamine | 100.0 | True
MOPS_NoC_L-Fucose | 100.0 | True
MOPS_NoC_D-Xylose | 100.0 | True
MOPS_NoC_Raffinose | 100.0 | True
MOPS_NoN_D-Methionine | 100.0 | True
MOPS_NoN_Biuret | 100.0 | True
MOPS_NoC_D-Cysteine | 100.0 | True
MOPS_NoC_Malate | 100.0 | True
ArgonneLB | 100.0 | True
MOPS_NoC_D-Glucose-6-phosphate | 100 | True
MOPS_NoC_L-Glutamine | 100.0 | True
MOPS_NoN_Thiourea | 100.0 | True
MOPS_NoN_L-Histidine | 100.0 | True
MOPS_NoN_Beta-Phenylethylamine | 100.0 | True
MOPS_NoC_Sorbate | 100.0 | True
MOPS_NoC_L-Arabinose | 100.0 | True


## Trimming out reactions by media

Now that we have a model that grows on everything, but should not, we should be able to remove individual reactions
to make the model grow only where we need it. This is probably going to be tricky!

We used this code to update our pickle file:

```
python flux_balance_scripts/update_pickle.py -p Citrobacter.67826.8.p -r temp/essential_rxns.py -v
```

Our model now has 1444 reactions in it.

