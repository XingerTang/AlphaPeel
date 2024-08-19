-----
Usage
-----

===============
Program options
===============

|Software| takes in several command line arguments to control the program's behaviour. To view a list of arguments, run |Software| without any command line arguments, i.e. ``AlphaPeel`` or ``AlphaPeel -h``. 

Input Arguments
---------------

::

    Input Options:
      -ped_file [PEDIGREE ...]
                          Pedigree file(s) (see format below).
      -geno_file [GENOTYPES ...]
                          Genotype file(s) (see format below).
      -seq_file [SEQFILE ...]
                          Sequence allele read count file(s) (see format below).
      -plink_file [BFILE ...]
                          Plink (binary) file(s).
      -start_snp START_SNP
                          The first marker to consider. The first marker is "1". Default: 1.
      -stop_snp STOP_SNP  The last marker to consider. Default: all markers considered.

|Software| requires a pedigree file (``-ped_file``) and one or more genomic data files to run the analysis.

|Software| supports the following genomic data files: genotype files in the AlphaGenes format (``-geno_file``), sequence allele read in the AlphaGenes format (``-seq_file``), and binary Plink files (``-plink_file``). Use of binary Plink files requires the package ``alphaplinkpython``, which  can be installed via ``pip``, but is only stable for Linux. There are known issues with this package, so we do not advocate its use at the moment.

Use the ``-start_snp`` and ``-stop_snp`` to run the analysis only on a subset of markers.

The input options in the form of ``[xxx ...]`` can take in more than one input file seperated by space.

Output Arguments 
----------------

::

    Output options:
      -out_file PREFIX      The output file prefix. All file outputs will be stored
                            as "PREFIX.dosage.txt" and so on.
      -out_id_order OUT_ID_ORDER
                            Determines the order in which individuals are ordered
                            in the output file based on their order in the
                            corresponding input file. Individuals not in the input
                            file are placed at the end of the file and sorted in
                            alphanumeric order. These individuals can be suppressed
                            with the "-out_id_only" option. Options: id, pedigree,
                            genotypes, sequence, segregation. Default: id.
      -out_id_only          Flag to suppress the individuals not present in
                            the file used with "-out_id_order". It also suppresses "dummy"
                            individuals.
      -n_io_thread N_IO_THREAD
                            Number of threads to use for input/output. Default: 1.


    Peeling output options:
      -no_dosage            Flag to suppress the dosage files.
      -no_param             Flag to suppress writing the model parameter files.
      -seg_prob             Flag to enable writing out the segregation probabilities.
      -phased_geno_prob     Flag to enable writing out the phased genotype probabilities.
      -geno_prob            Flag to enable writing out the genotype probabilities.
      -hap                  Flag to call and write out the haplotypes.
      -geno                 Flag to call and write out the genotypes.
      -geno_threshold [GENO_THRESHOLD ...]
                            Genotype calling threshold(s). Multiple space separated values allowed.
                            Value less than 1/3 will be replaced by 1/3.
      -hap_threshold [HAP_THRESHOLD ...]
                            Haplotype calling threshold(s). Multiple space separated values allowed.
                            Value less than 1/2 will be replaced by 1/2.
      -binary_call_file    Flag to write out the called genotype files as a
                            binary plink output [Not yet implemented].

By default |Software| produces a dosage file and two model parameter files (genotype error rate and recombination rate). Creation of these files can be suppressed with the ``-no_dosage``, and ``-no_param`` options. |Software| can also write out the phased genotype probability file (.phased_geno_prob.txt) with the `-phased_geno_prob` argument and the segregation probability file (.seg_prob.txt) with the `-seg_prob` argument.

The ``-geno_threshold`` and ``-hap_threshold`` arguments respectively control control which genotypes and haplotypes are called. A threshold of 0.9 will give calls only if the probability mass for one genotype (or haplotype) is higher than 0.9. Using a higher-value will increase the accuracy of called genotypes (or haplotypes), but will result in fewer called genotypes (or haplotypes). Since there are three genotypes states and two haplotype states, "best-guess" genotypes and haplotypes are respectively called with a threshold less than ``1/3`` and ``1/2``.

``-binary_call_file`` option can be used to change the output to a plink binary format.

The order in which individuals are output can be changed by using the ``out_id_order`` option. This option changes the order in which individuals are written out to the order in which they were observed in the corresponding file. The ```-out_id_only`` option suppresses the output of dummy individuals (not recommended for hybrid peeling).

The argument ``-n_io_thread`` controls the number of threads/processes used by |Software|. |Software| uses additional threads to parse and format input and output files. Setting this option to a value greater than 1 is only recommended for very large files (i.e. >10,000 individuals).

Peeling arguments 
------------------

::

    Mandatory peeling arguments:
      -method METHOD        Program run type. Either "single" or "multi".
    
    Optional peeling arguments:
      -n_cycle N_CYCLE    Number of peeling cycles. Default: 5.
      -n_thread N_THREAD
                            Number of threads to use. Default: 1.
      -rec_length REC_LENGTH
                            Estimated recombination length of the chromosome in Morgans.
                            [Default 1.00]

    Peeling control arguments:
      -est_geno_error_prob  Flag to re-estimate the genotyping error rates after
                            each peeling cycle.
      -est_seq_error_prob   Flag to re-estimate the sequencing error rates after
                            each peeling cycle.
      -est_rec_prob         Flag to re-estimate the recombination rates after
                            each peeling cycle.
      -est_alt_allele_prob  Flag to re-estimate the alternative allele probabilities after
                            each peeling cycle.
      -no_phase_founder    A flag phase a heterozygous allele in one of the
                            founders (if such an allele can be found).
      -sex_chrom            A flag to indicate that input data is for a sex chromosome. Sex needs to
                            be given in the pedigree file. This is currently an
                            experimental option.

    Genotype probability arguments:
      -geno_error_prob GENO_ERROR_PROB
                            Genotyping error rate. [Default 0.0001]
      -seq_error_prob SEQ_ERROR_PROB
                            Sequencing error rate. [Default 0.001]

``-method`` controls whether the program is run in "single-locus" or "multi-locus" model. Single locus mode does not use linkage information to perform imputation. It is fast, but not very accurate. Multi-locus mode runs multi-locus iterative peeling which uses linkage information to increase accuracy and calculate segregation values.

For hybrid peeling, where a large amount (millions of segregating sites) of sequence allele read counts needs to be imputed, first run the program in multi-locus mode to generate a segregation file, and then run the program in single-locus mode with a known segregation file.

The ``-geno_error_prob``, ``-seq_error_prob`` and ``-rec_length`` arguments control some of the model parameters used in the model. ``-seq_error_prob`` must not be zero. |Software| is robust to deviations in genotyping error rate and sequencing error rate so it is not recommended to use these options unless large deviations from the default are known. Changing the ``-length`` argument to match the genetic map length can increase accuracy in some situations.

The ``-est_geno_error_prob`` and ``-est_seq_error_prob`` options estimate the genotyping error rate and the sequencing error rate based on miss-match between observed and inferred states. This option is generally not necessary and can increase runtime. ``-est_alt_allele_prob`` estimates the alternative allele probabilities after each peeling cycle. This option can be useful if there are a large number of non-genotyped founders.

Hybrid peeling arguments 
------------------------

::

    Single locus arguments:
      -seg_file SEG_FILE    A segregation probabilities file for hybrid peeling.
      -seg_map_file SEG_MAP_FILE
                            A map file for loci in the segregation probabilities file.
      -map_file MAP_FILE    A map file for all loci in hybrid peeling.

In order to run hybrid peeling the user needs to supply a ``-map_file`` which gives the genetic positions for the SNPs in the sequence allele read counts data supplied, a ``-seg_map_file`` which gives the genetic position for the SNPs in the segregation file, and a ``-seg_file`` which gives the segregation values generated via multi-locus iterative peeling. These arguments are not required for running in multi-locus mode.

============
File formats
============

Input file formats
------------------

Pedigree file
=============

Each line of a pedigree file has three values, the individual's id, their father's id, and their mother's id. "0" represents an unknown id.

Example:

::

  id1 0 0
  id2 0 0
  id3 id1 id2
  id4 id1 id2

Genotype file 
=============

Genotype files contain the input genotypes for each individual. The first value in each line is the individual's id. The remaining values are the genotypes of the individual at each locus, either 0, 1, or 2 (or 9 if missing). The following examples gives the genotypes for four individuals genotyped on four markers each.

Example:

::

  id1 0 2 9 0 
  id2 1 1 1 1 
  id3 2 0 2 0 
  id4 0 2 1 0

Sequence allele read counts file
================================

The sequence allele read counts file has two lines for each individual. The first line gives the individual's id and read counts for the reference allele. The second line gives the individual's id and allele read counts for the alternative allele.

Example:

::

  id1 4 0 0 7 # Reference allele for id1
  id1 0 3 0 0 # Alternative allele for id1
  id2 1 3 4 3
  id2 1 1 6 2
  id3 0 3 0 1
  id3 5 0 2 0
  id4 2 0 6 7
  id4 0 7 7 0

Binary plink file
=================

Binary Plink files are supported using the package ``AlphaPlinkPython``. The pedigree supplied by the ``.fam`` file will be used if a pedigree file is not supplied. Otherwise, the pedigree file will be used and the ``.fam`` file will be ignored. 

Map file 
========

The map file gives the chromosome number, the marker name, and the base pair position for each marker in two columns. Only markers on one chromosome should be provided! 

Example:

::

  1 snp_a 12483939
  1 snp_b 192152913
  1 snp_c 65429279
  1 snp_d 107421759


Output file formats
-------------------

Phase file
==========

The phase file gives the phased haplotypes (either 0 or 1) for each individual in two lines. For individuals where we can determine the haplotype of origin, the first line will provide information on the paternal haplotype, and the second line will provide information on the maternal haplotype.

Example:

::

  id1 0 1 9 0 # Paternal haplotype
  id1 0 1 9 0 # Maternal haplotype
  id2 1 1 1 0
  id2 0 0 0 1
  id3 1 0 1 0
  id3 1 0 1 0 
  id4 0 1 0 0
  id4 0 1 1 0

Genotype probability file
=========================

The haplotype file (*.phased_geno_prob.txt*) provides the (phased) allele probabilities for each locus. There are four lines per individual containing the allele probability for the (aa, aA, Aa, AA) alleles where the paternal allele is listed first, and where *a* is the reference (or major) allele and *A* is the alternative (or minor) allele.

Example:

::

  id1    0.9998    0.0001    0.0001    1.0000
  id1    0.0000    0.4999    0.4999    0.0000
  id1    0.0000    0.4999    0.4999    0.0000
  id1    0.0001    0.0001    0.0001    0.0000
  id2    0.0000    1.0000    0.0000    1.0000
  id2    0.9601    0.0000    0.0455    0.0000
  id2    0.0399    0.0000    0.9545    0.0000
  id2    0.0000    0.0000    0.0000    0.0000
  id3    0.9998    0.0001    0.0001    1.0000
  id3    0.0000    0.4999    0.4999    0.0000
  id3    0.0000    0.4999    0.4999    0.0000
  id3    0.0001    0.0001    0.0001    0.0000
  id4    1.0000    1.0000    0.0000    1.0000
  id4    0.0000    0.0000    0.0000    0.0000
  id4    0.0000    0.0000    0.0000    0.0000
  id4    0.0000    0.0000    1.0000    0.0000

Dosage file
===========

The dosage file gives the expected allele dosage for the alternative (or minor) allele for each individual. The first value in each line is the individual ID. The remaining values are the allele dosages at each loci. These values will be between 0 and 2.

Example:

::

  1    0.0003    1.0000    1.0000    0.0001
  2    1.0000    0.0000    1.0000    0.0000
  3    0.0003    1.0000    1.0000    0.0001
  4    0.0000    0.0000    2.0000    0.0000

Segregation file
================

The segregation file gives the joint probability of each pattern of inheritance. There are four lines for each individual representing the probability of inheriting: 

  1. the grand **paternal** allele from the father and the grand **paternal** allele from the mother
  2. the grand **paternal** allele from the father and the grand **maternal** allele from the mother
  3. the grand **maternal** allele from the father and the grand **paternal** allele from the mother
  4. the grand **maternal** allele from the father and the grand **maternal** allele from the mother

Example:

::

  id1    1.0000    0.9288    0.9583    0.9834
  id1    0.0000    0.0149    0.0000    0.0000
  id1    0.0000    0.0554    0.0417    0.0166
  id1    0.0000    0.0009    0.0000    0.0000
  id2    0.9810    0.9842    1.0000    0.9971
  id2    0.0174    0.0158    0.0000    0.0013
  id2    0.0016    0.0000    0.0000    0.0016
  id2    0.0000    0.0000    0.0000    0.0000
  id3    0.0164    0.0149    0.0000    0.0065
  id3    0.9259    0.9288    0.9582    0.9769
  id3    0.0010    0.0009    0.0000    0.0001
  id3    0.0567    0.0554    0.0417    0.0165
  id4    0.0002    0.0000    0.0002    0.0004
  id4    0.0015    0.0000    0.0019    0.0041
  id4    0.1189    0.1179    0.1052    0.0834
  id4    0.8794    0.8821    0.8927    0.9122

Model parameter files
=====================

|Software| outputs three model parameter files, ``.alt_allele_prob.txt``, ``.seq_error_prob.txt``, ``.geno_error_prob.txt``, ``.rec_prob.txt``. These give the minor allele frequency, sequencing error rates, genotyping error rates and recombination rates used. All three files contain a single column with an entry for each marker.

Example ``.alt_allele_prob.txt`` file for four loci:

::

  0.468005
  0.195520
  0.733061
  0.145847

----------------
Code explanation
----------------

The main peeling funciton of |Software| is given by ``tinypeel.Peeling.Peeling.peel()`` function:

.. autofunction:: tinypeel.Peeling.Peeling.peel

The peeling process consists of two parts:

1. The first part performs a Baum-Welch-like algorithm over each family of each generation of the pedigree.
   The details of the HMM are the following:
   
   - Hidden states: phased genotype.
   - Time dimension: generation (parent -> child).
   - Observed states: input genotype and sequence data.
   - Emission function: combined with the observed states, are introduced via the ``penetrance`` of the peeling information container, which is generated from ``tinypeel.tinyhouse.ProbMath.getGenotypeProbabilities()``.
     
     - For genotype input: suppose genotype error is :math:`e`, the following table represents how the input genotype data is encoded to be the probabilities of the phased genotype data, with rows representing the phased genotypes and columns representing the input genotype data:

       .. list-table::
          :header-rows: 1

          * - 
            - aa
            - aA
            - Aa
            - AA
          * - 0
            - 1 - e
            - e/2
            - e/2
            - e/2
          * - 1
            - e/2
            - 1 - e
            - 1 - e
            - e/2
          * - 2
            - e/2
            - e/2
            - e/2
            - 1 - e

     - For sequence input: suppose sequence error is :math:`e`, the following table represents how the input sequence data is encoded to be the phased genotype data, with :math:`r` representing the number of reference alleles and :math:`a` representing the number of alternative alleles.

       .. list-table::
          :header-rows: 1

          * - phased genotype
            - corresponding probability
          * - aa
            - :math:`(1 - e)^r e^a`
          * - aA
            - :math:`\left(\frac{1}{2}\right)^{r + a}`
          * - Aa
            - :math:`\left(\frac{1}{2}\right)^{r + a}`
          * - AA
            - :math:`e^r (1 - e)^a`

      All the values are normalised before use. In the case that both sequence and genotype inputs are provided, the two probabilities are multiplied together.

   - Transmission function: implemented via ``segregationTensor`` defined in ``tinypeel.tinyhouse.ProbMath()``. The ``segregationTensor`` is a 4D numpy array of float32 of size :math:`4 \times 4 \times 4 \times 4`, with 

     * the first dimension represents the paternal genotype (``p``),

     * the second dimension represents the maternal genotype (``m``),
    
     * the third dimension represents the child's genotype (``allele``), and
    
     * the last dimension represents the child's segregation (``seg``).

     Mathematically, ``segregationTensor`` represents the probability of each combination of the sire's genotype, the dam's genotype,
     child's genotype and segregation without any other information (:math:`P(p, m, allele, seg)`). The following are two examples:

     * Example 1: Suppose ``allele = 0`` (genotype = aa) and ``seg = 0`` (segregation = pp), the values of ``segregationTensor[:, :, 0, 0]`` are the following

     .. list-table::
          :header-rows: 1

          * - 
            - m = aa
            - m = aA
            - m = Aa
            - m = AA
          * - p = aa
            - 1
            - 1
            - 0
            - 0
          * - p = aA
            - 1
            - 1
            - 0
            - 0
          * - p = Aa
            - 0
            - 0
            - 0
            - 0
          * - p = AA
            - 0
            - 0
            - 0
            - 0

     * Example 2: Suppose ``allele = 2`` (genotype = Aa) and ``seg = 1`` (segregation = pm), the values of ``segregationTensor[:, :, 2, 1]`` are the following

     .. list-table::
          :header-rows: 1

          * - 
            - m = aa
            - m = aA
            - m = Aa
            - m = AA
          * - p = aa
            - 0
            - 0
            - 0
            - 0
          * - p = aA
            - 0
            - 0
            - 0
            - 0
          * - p = Aa
            - 1
            - 0
            - 1
            - 0
          * - p = AA
            - 1
            - 0
            - 1
            - 0

     A pre-defined error :math:`e` is used here with ``segregationTensor = segregationTensor*(1-e) + e/4``. (mutation error)

     This matrix can be used to generate ``ChildSegs``, which is the matrix controls how the information is passed across generations. By ``tinypeel.Peeling.Peeling.createChildSegs()``, the probabilities of each combination of sire's genotype, dam's genotype and child's genotype can be calculated via summing over the child's segregation states. 

     .. autofunction:: tinypeel.Peeling.Peeling.createChildSegs

     The information are passed with functions ``tinypeel.Peeling.Peeling.projectChildGenotypes()`` and ``tinypeel.Peeling.Peeling.projectParentGenotypes()``.

     .. autofunction:: tinypeel.Peeling.Peeling.projectChildGenotypes

     .. autofunction:: tinypeel.Peeling.Peeling.projectParentGenotypes

     The pre-defined error :math:`e` is used in a similar way as the ``segregationTensor`` on the following variables:
       
     - ``JointParents``,
     - ``probSire`` and ``probDam``,
     - ``childValues``, and
     - ``sirePosterior`` and ``damPosterior``,
     
     which are all intermediate values or the results of the transmission across generation.

   - Probabilities update: The phased genotype probabilities are calculatd via ``anterior * penetrance * posterior``, which

     * ``anterior``: is updated when peeling down, every time a new value is calculated.
     * ``penetrance``: depends only on input.
     * ``posterior``: is updated when peeling up, but only when a peeling cycle is finished.
  
2. The second part performs a Baum-Welch-like algorithm over each locus of each individual. This part is implemented only when the multi-locus peeling mode is used.
   The details of the HMM is the following:

   - Hidden states: segregation states.
   - Time dimension: locus (locus_i -> locus_i+1).
   - Observed states: phased genotype probabilities from part 1.
   - Emission function: A segregation estimate ``pointSeg`` is obtained by ``estimateSegregationWithNorm()``,

     .. autofunction:: tinypeel.Peeling.Peeling.estimateSegregationWithNorm
   
     then same usage of :math:`e` as the ``segregationTensor`` is applied, but now on the resulting output of function for the Baum-Welch algorithm implementation:

     * matched segregation: :math:`1 - \frac{3}{4}e`,
     * unmatched segregation: :math:`\frac{1}{4}e`.
  
   - Transmission function: first generate equal recombination rates `r` across all loci, which the values are calculated by assuming there is exactly 1 recombination happened per snippet of input and the distances between each locus are equal. The transmission function is setting up as the following:

     * if the segregation states are same at locus i and locus i + 1: :math:`(1 - r)^2`,
     * if the segregation states are different by one parent at locus i and locus i + 1: :math:`r \times (1 - r)`,
     * if the segregation states are different by both parents at locus i and locus i + 1: :math:`r^2`.

=====================
Possible improvements
=====================

1. Improvement in the emission function of the part 1:

   - genotype input, more realistic matrix can be used.
   - the updates of errors may can make use of the Baum-Welch algorithm.
2. Improvement in the error term of the transmission funciton of the part 1:

   - single application of ``e`` can be used
   - the value of ``e`` can also be updated based on the Baum-Welch algorithm.
3. Improvement in the posterior update of part 1: update more frequently. (x)
4. Improvement in the transmission funciton of part 2: the updates of recombination rates may can make use of the Baum-Welch algorithm.
5. Improvement in the emission function of part 2: 

   - the error term should be applied before the implementation of the Baum-Welch algorithm.
   - more realistic values can be used.

New: 

6. Implementation of W-shape cycles update to fully use the information of generations with dense information. 
7. Implementation of random update.
7. Implementation of a logic filter of phasing before the computation.


.. |Software| replace:: AlphaPeel
