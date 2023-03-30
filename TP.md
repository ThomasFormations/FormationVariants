## TP détection de variants

L'animal SALSA semble présenter un phénotype atypique, une absence de cornes.
Il y a une forte présomption que cette anomalie phénotypique soit d'origine génétique.
Par ailleurs deux autres animaux, ZOUK et TANGO, semblent présenter un phénotype intermédiaire, une seule corne (à gauche).

Pour en avoir le coeur net, on procède à du séquençage WGS de ces animaux pour identifier les variants.
Afin de pouvoir contraster les résultats avec un animal à phénotype normal (avec cornes) nous inégrons l'animal MAMBO qui possède des cornes et que l'on considère comme phénotype "sauvage" (wild type), c'est à dire l'animal contrôle.

Les différentes étapes ci-dessous vont nous guider dans ce travail.

#### Création de l'environement conda.
```
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/environment.yml
mamba env create -f environment.yml
conda activate polled
```

#### Création de l'arborescence
```
mkdir -p data/ref
mkdir -p data/fastq
mkdir -p mapping
mkdir -p variants
```
### Récupération des données
```
cd data/ref
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/reference.fa.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/reference.gff3.bgz
cd ../fastq
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/MAMBO_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/MAMBO_R2.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/SALSA_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/SALSA_R2.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/TANGO_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/TANGO_R2.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/ZOUK_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/ZOUK_R2.fastq.gz
cd ../..
```

### Qualité des lectures
1. Vérifier que le format des fichiers fastq correspond bien au format attendu ?
3. Quelle est la signification de la mesure de qualité d'une base lue ? Comment est-elle codée dans le fichier fastq ?
4. Quelles sont les différents scores observés pour les 4 premières lectures du fichier SALSA_R1.fastq.gz ?
5. A quelle probabilité d'erreur cela correspond-il ?  
(Vous pourrez jeter un coup d'oeil à cette page https://people.duke.edu/~ccc14/duke-hts-2018/bioinformatics/quality_scores.html)

A l'aide du logiciel FASTQ (https://www.bioinformatics.babraham.ac.uk/projects/fastqc) déjà dans votre environnement, vérifiez la qualité des lectures
```
fastqc data/fastq/SALSA_R1.fastq.gz data/fastq/SALSA_R2.fastq.gz
```
Pour l'ensemble des données
```bash
for sample in SALSA MAMBO TANGO ZOUK
do
  fastqc data/fastq/${sample}_R1.fastq.gz data/fastq/${sample}_R2.fastq.gz 
done

```
Les fichiers de synthèse de ces analyses se trouvent dans le répertoire ```data/fastq```

5. Est-ce que cela correspond à ce qui avait été observé pour les 6 premières lectures de Salsa ?


### Alignement des lectures sur la référence

#### Indexation de la séquence de référence
L'outil ```samtools faidx``` 
```
samtools faidx data/ref/reference.fa.gz
```
#### Indexation de la séquence de référence pour bwa
```
bwa index data/ref/reference.fa.gz
```
6. Quelle est la différence entre ces deux types d'indexation ?

#### Alignement des lectures
```
bwa mem  -R "@RG\tID:MAMBO\tLB:MAMBO\tPL:ILLUMINA\tSM:MAMBO" \
                 data/ref/reference.fa.gz \
                 data/fastq/MAMBO_R1.fastq.gz \
                 data/fastq/MAMBO_R2.fastq.gz > mapping/MAMBO.sam
```
Jeter un coup d'oeil aux fichiers sam. 

7. Jeter un coup d'oeil au format SAM ?
8. Quelle est la signification de la mesure de qualité d'un alignement ?
9. Comment procéder pour utiliser 4 processeurs au lieu d'un seul ? Et en utilisant le nombre maximum de processeurs de votre machine ?
10. Comment, selon vous, cette parallèlisation se fait-elle ?

Afin de pouvoir travailler avec le fichier sortie sam, on procède à une transformation en fichier .bam, à son tri et à son indexation.

```
samtools sort mapping/MAMBO.sam -OBAM -o mapping/MAMBO.bam
samtools index mapping/MAMBO.bam
```
11. Quelle est la différence entre le format bam et le format sam ?
12. Que permet ici l'indexation ?

13. Proposer à l'aide d'un pipe, une solution permettant d'éviter le passage par un fichier intermédiaire .sam.  

Les utilitaires samtools stats et samtools flagstats permettent d'obtenir des statistiques sur les alignemnts
```
samtools stats mapping/MAMBO.bam > mapping/MAMBO.bam.stats
samtools flagstats mapping/MAMBO.bam > mapping/MAMBO.bam.flagtstats
```
14. En vous inspirant de la boucle for ci-dessus, écrire quelques lignes de codes permettant de réaliser tous les alignements, et les fichiers d'index associés et les statistiques.

15. En utilisant l'outil ```multiqc``` (dans votre environnement) produire les fichiers de synthèse des différentes statistiques (cf la documentation de multiqc, https://multiqc.info/docs).

16. A l'aide de l'outil IGV visualiser les alignements de SALSA dans la région 28:6,050,303-6,053,006. Plusieurs alignements sont atypiques (de couleurs différentes), quelle en est la raison ?

### Détection de variants

Pour la détection de variants nous utiliserons deux approches, en utilisant d'une par l'outil ```bcftools``` (commandes ```mpileup``` et ```calls```) et d'autre part l'outil ```freebayes```.

Les deux approches permettent de détecter simultanément les variants sur l'ensemble des individus.

#### bcftools

```
ls mapping/*.bam > bamlist.txt
bcftools mpileup -a AD -Ou -f data/ref/reference.fa.gz --bam-list bamlist.txt | bcftools call -mv -Ov -o variants/bcftools_calls.vcf
bcftools stats variants/bcftools_calls.vcf > variants/bcftools_calls.vcf.stats
```
17. Quel est le format d'un fichier vcf ?
18. Que contient le champ de FORMAT PL ?
19. En utilisant la commande ```query``` de bcftools afficher dans le terminal, pour chaque variant les informations suivantes
```
CHROM POS QUAL AC G1 G2 G3 G4
```
chromosome, position, mesure de qualité du variant, AC le nombre d'allèles alternatifs et les génotypes pour les 4 individus.

#### freebayes

```
zcat data/ref/reference.fa.gz > data/ref/reference.fa
samtools faidx  data/ref/reference.fa
freebayes -f data/ref/reference.fa mapping/MAMBO.bam mapping/SALSA.bam mapping/TANGO.bam mapping/ZOUK.bam > variants/freebayes_calls.vcf
bcftools stats variants/freebayes_calls.vcf > variants/freebayes_calls.vcf.stats
```
20. Que peut-on dire de la différence du nombre de variants entre les deux approches ?

On propose de filter les variants sur la qualité
```
cat variants/freebayes_calls.vcf | vcffilter  -f "QUAL > 20" > variants/freebayes_calls.q20.vcf
bcftools stats variants/freebayes_calls.q20.vcf > variants/freebayes_calls.q20.vcf.stats
```
21. Quelle est la signification de la mesure de qualité d'un variant ?

### Annotation des variants

#### Récupération de la base de données du génome bovin

```
# snpEff download -v ARS-UCD1.2.105
```
ou
```
mkdir -p data/snpEffdatabases
cd data/snpEffdatabases
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/ARS-UCD1.2.105.bz2
tar xjf ARS-UCD1.2.105.bz2
cd ../..
```

#### Annotation

```
snpEff -v ARS-UCD1.2.105 -datadir $PWD/data/snpEffdatabases variants/bcftools_calls.vcf > variants/bcftools_calls.annotated.vcf
```

```
snpEff -v ARS-UCD1.2.105 -datadir $PWD/data/snpEffdatabases variants/freebayes_calls.q20.vcf > variants/freebayes_calls.q20.annotated.vcf
```

Le logiciel snpEff fournit en sortie un fichier de synthèse.

22. Quel est le variant le plus sévère détecté par les deux logiciels ? (indice rechercher une mutation non-sens) 
A l'aide du logiciel igv, regader au voisinage de ce variant, les lectures (fichiers bam) et les génotypes (fichier vcf).

23. Que peut-on dire de ces génotypes ?

#### Supplementary materials

24. Le programme de détection de SNP bcftools se base sur une sortie dite mpileup

```
samtools mpileup mapping/SALSA.bam --reference data/ref/reference.fa > SALSA.mpileup
```

25. Extraire les lignes du fichier de sortie qui indiquent des différences par rapport à la référence (on pourra utiliser awk avec une expression régulière).

26. En vous aidant de la documentation de [samtools-mpileup](http://www.htslib.org/doc/samtools-mpileup.html), expliquer ce que renvoie la fonction python suivante
27. Comme peut-on utiliser cette fonction pour calculer la log-vraisemblance des génotypes ?


```python
def parse_mpileup_line(sequence, phred, info)):
    """
    parse a line of an mpileup output 

    Parameters:
    mpileup_line -- a line of an mipileup output as a string
    """
    mpileup = mpileup_line.split("\t")
    chrom = mpileup[0]
    pos = mpileup[1]
    ref = mpileup[2]
    coverage = mpileup[4]
    sequence = mpileup[5]
    phred = mpileup[6]
    
    # dictionary containing all (Phred score) tuple for each
    # base of the sequence
    letters = defaultdict(list)
    nucleotides = ['A', 'C', 'G', 'T']

    sequence = sequence.upper()

    # Offset used to obtain the Phred score
    asciiOffset = 33

    si = 0  # position in the sequence
    pi = 0  # position in the sequence of phred scores

    while(si < len(sequence)):
        currentchar = sequence[si]
        if currentchar == "," or currentchar == ".":
            phredVal = ord(phred[pi]) - asciiOffset
            letters[ref].append(phredVal)
            si += 1
            pi += 1
        elif currentchar in nucleotides:
            phredVal = ord(phred[pi]) - asciiOffset
            letters[currentchar].append(phredVal)
            si += 1
            pi += 1
        elif currentchar == "+" or currentchar == "-":
            res = match(r"[+-](\d+)", sequence[si:len(sequence)])
            indelLength = res.groups()[0]
            si = si + 1 + len(indelLength) + int(indelLength)
            # No change in positionPhred because there is not Phred value for an indel
        elif currentchar == ">" or currentchar == "<":
            # Reference skip (CIGAR "N")
            # Phred value available but skipped
            si += 1
            pi += 1
        elif currentchar == "$":
            si += 1
            # No change in positionPhred because there is not Phred value for
            # the end of a sequence
        elif currentchar == "*":
            si += 1
            # No change in positionPhred because there is not Phred value for
            # the shadow of a deletion
        elif currentchar == "^":
            # The beginning is composed of 2 characters :
            #     "^" followed by the Mapq value
            si += 2
            # No change in positionPhred because there is not Phred value
            # for the beginning of a sequence
        else:
            print("Problem with letter: " + currentchar)
            print("in this cigar string: " + sequence)
            sys.exit(2)

    return(letters)
```

