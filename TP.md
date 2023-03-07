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
cd ../data/fastq
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/MAMBO_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/MAMBO_R2.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/SALSA_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/SALSA_R2.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/TANGO_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/TANGO_R2.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/ZOUK_R1.fastq.gz
wget https://genoweb.toulouse.inra.fr/~faraut/FormationM12023/TP/ZOUK_R2.fastq.gz
```
### Qualité des lectures

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
Quelle est la différence entre ces deux types d'inexation ?

#### Alignement des lectures
```
bwa mem  -R "@RG\tID:MAMBO\tLB:MAMBO\tPL:ILLUMINA\tSM:MAMBO" \
                 data/ref/reference.fa.gz \
                 data/fastq/MAMBO_R1.fastq.gz \
                 data/fastq/MAMBO_R2.fastq.gz > mapping/MAMBO.sam
```
Comment procéder pour utiliser 4 processeurs au lieu d'un seul ? Et en utilisant le nombre maximum de processeurs de votre machine ?
Comme, selon vous, cette parallèlisation se fait-elle ?

Afin de pouvoir travailler avec le fichier sortie sam, on procède à son tri et à une indexation après transformation en fichier bam
```
samtools sort mapping/MAMBO.sam -OBAM -o mapping/MAMBO.bam
samtools index mapping/MAMBO.bam
```
Proposer à l'aide d'un pipe, une solution permettant d'éviter le passage par un fichier intermédiaire .sam.  

Les utilitaires samtools stats et samtools flagstats permettent d'obtenir des statistiques sur les alignemnts
```
samtools stats mapping/MAMBO.bam > mapping/MAMBO.bam.stats
samtools flagstats mapping/MAMBO.bam > mapping/MAMBO.bam.flagtstats
```
En vous inspirant de la boucle for ci-dessus, écrire quelques lignes de codes permettant de réaliser tous les alignements, es les fichiers d'index associés et les statistiques.

En utilisant l'outil ```multiqc``` (dans votre environnement) produire les fichiers de synthèse des différentes statistiques (cf la documentation de multiqc, https://multiqc.info/docs).

### Détection de variants

Pour la détection de variants nous utiliserons deux approches, en utilisant d'une par l'outil ```bcftools``` (commandes ```mpileup``` et ```calls``` et d'autre part l'outil ```freebayes```.

Les deux approches permettent de détecter simultanément les variants sur l'ensemble des individus.

#### bcftools

```
ls mapping/*.bam > bamlist.txt
bcftools mpileup -Ou -f data/ref/reference.fa.gz --bam-list bamlist.txt | bcftools call -mv -Ov -o variants/bcftools_calls.vcf
```

#### freebayes

```
freebayes -f data/ref/reference.fa.gz mapping/MAMBO.bam mapping/SALSA.bam mapping/TANGO.bam mapping/ZOUK.bam > variants/freebayes_calls.vcf
```

