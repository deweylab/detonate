// This file is autogenerated. Edit the template instead.
std::string get_help_string() { return
"  REF-EVAL: A toolkit of reference-based scores for de novo transcriptome\n"
"                       sequence assembly evaluation\n"
"\n"
"Overview\n"
"\n"
"   REF-EVAL computes a number of reference-based scores. These scores\n"
"   measure the quality of a transcriptome assembly relative to a\n"
"   collection of reference sequences. For information about how to run\n"
"   REF-EVAL, see \"Usage\" and the sections following it below. For\n"
"   information about the score definitions, see \"Score definitions\" and\n"
"   the sections following it below.\n"
"\n"
"Usage\n"
"\n"
"   As an optional first step, estimate the \"true\" assembly, using\n"
"   [1]REF-EVAL-ESTIMATE-TRUE-ASSEMBLY. Alternatively, you can use the\n"
"   full-length reference transcript sequences directly as a reference.\n"
"\n"
"   From now on, we will call the estimated \"true\" assembly or the\n"
"   collection of full-length reference sequences (whichever you choose\n"
"   to use) the reference. Let's assume that the assembly of interest is\n"
"   in A.fa, and the reference is in B.fa.\n"
"\n"
"   If you want to compute alignment-based scores (see --scores below for\n"
"   more info), align the assembly to the reference and vice versa using\n"
"   [2]Blat. We recommend fairly unrestrictive settings, in order to\n"
"   generate many candidate alignments.\n"
"\n"
" $ blat -minIdentity=80 B.fa A.fa A_to_B.psl\n"
" $ blat -minIdentity=80 A.fa B.fa B_to_A.psl\n"
"\n"
"   If you want to compute weighted variants of scores, use [3]RSEM (or\n"
"   [4]RSEM-EVAL) to compute the expression of the assembly and reference\n"
"   relative to the given reads. Let's assume that the reads are in\n"
"   reads.fq.\n"
"\n"
" $ rsem-prepare-reference --no-polyA A.fa A_ref\n"
" $ rsem-prepare-reference --no-polyA B.fa B_ref\n"
" $ rsem-calculate-expression -p 24 --no-bam-output reads.fq A_ref A_expr\n"
" $ rsem-calculate-expression -p 24 --no-bam-output reads.fq B_ref B_expr\n"
"\n"
"   Finally, run REF-EVAL. To compute everything, run:\n"
"\n"
" $ ./ref-eval --scores=nucl,pair,contig,kmer,kc \\\n"
"              --weighted=both \\\n"
"              --A-seqs A.fa \\\n"
"              --B-seqs B.fa \\\n"
"              --A-expr A_expr.isoforms.results \\\n"
"              --B-expr B_expr.isoforms.results \\\n"
"              --A-to-B A_to_B.psl \\\n"
"              --B-to-A B_to_A.psl \\\n"
"              --num-reads 5000000 \\\n"
"              --readlen 76 \\\n"
"              --kmerlen 76 \\\n"
"              | tee scores.txt\n"
"\n"
"   To only compute the kmer compression score (and its dependencies),\n"
"   run:\n"
"\n"
" $ ./ref-eval --scores=kc \\\n"
"              --A-seqs A.fa \\\n"
"              --B-seqs B.fa \\\n"
"              --B-expr B_expr.isoforms.results \\\n"
"              --num-reads 5000000 \\\n"
"              --readlen 76 \\\n"
"              --kmerlen 76 \\\n"
"              | tee scores.txt\n"
"\n"
"   To only compute the unweighted reference-based scores, run:\n"
"\n"
" $ ./ref-eval --scores=nucl,pair,contig \\\n"
"              --weighted=no \\\n"
"              --A-seqs A.fa \\\n"
"              --B-seqs B.fa \\\n"
"              --A-to-B A_to_B.psl \\\n"
"              --B-to-A B_to_A.psl \\\n"
"              | tee scores.txt\n"
"\n"
"   To only compute the scores discussed in the DETONATE paper, run:\n"
"\n"
" $ ./ref-eval --paper \\\n"
"              --A-seqs A.fa \\\n"
"              --B-seqs B.fa \\\n"
"              --B-expr B_expr.isoforms.results \\\n"
"              --A-to-B A_to_B.psl \\\n"
"              --B-to-A B_to_A.psl \\\n"
"              --num-reads 5000000 \\\n"
"              --readlen 76 \\\n"
"              --kmerlen 76 \\\n"
"              | tee scores.txt\n"
"\n"
"   The scores will be written to standard output (hence, above, to\n"
"   scores.txt). Progress information is written to standard error.\n"
"   Further details about the arguments to REF-EVAL are described below.\n"
"   Further details about the scores themselves are given under \"Score\n"
"   definitions\" below.\n"
"\n"
"Usage: Score specification\n"
"\n"
"   --scores arg\n"
"\n"
"           The groups of scores to compute, separated by commas (e.g.,\n"
"           --scores=nucl,contig,kc). It is more efficient to compute all\n"
"           the scores you are interested in using one invocation of\n"
"           REF-EVAL instead of using multiple invocations that each\n"
"           compute one score. The available score groups are as follows:\n"
"\n"
"           Alignment-based score groups:\n"
"\n"
"              * nucl: nucleotide precision, recall, and F1.\n"
"              * contig: contig precision, recall, and F1.\n"
"              * pair: pair precision, recall, and F1.\n"
"\n"
"           Alignment-free score groups:\n"
"\n"
"              * kmer: kmer Kullback-Leibler divergence, Jensen-Shannon\n"
"                divergence, and Hellinger distance.\n"
"              * kc: kmer recall, number of nucleotides, and kmer\n"
"                compression score.\n"
"\n"
"           Required unless --paper is given.\n"
"\n"
"   --weighted arg\n"
"\n"
"           A string indicating whether to compute weighted or unweighted\n"
"           variants of scores, or both (e.g., --weighted=yes):\n"
"\n"
"              * yes: compute weighted variants of scores.\n"
"              * no: compute unweighted variants of scores.\n"
"              * both: compute both weighted and unweighted variants of\n"
"                scores.\n"
"\n"
"           In weighted variants, the expression levels (TPM) of the\n"
"           assembly and reference sequences are taken into account, and\n"
"           hence need to be specified using --A-expr and --B-expr.\n"
"           Unweighted variants are equivalent to weighted variants with\n"
"           uniform expression.\n"
"\n"
"           The distinction between weighted and unweighted variants\n"
"           doesn't make sense for the KC score, so this option is\n"
"           ignored by the KC score.\n"
"\n"
"           Required unless --paper or only --score=kc is given.\n"
"\n"
"   --paper\n"
"\n"
"           As an alternative to the above, if you are only interested in\n"
"           computing the scores described in the main text of our paper\n"
"           [1], you can pass the --paper flag instead of the --scores\n"
"           and --weighted options. In that case, the following scores\n"
"           will be computed:\n"
"\n"
"           Alignment-based scores:\n"
"\n"
"              * unweighted nucleotide F1\n"
"              * unweighted contig F1\n"
"\n"
"           Alignment-free score groups:\n"
"\n"
"              * weighted kmer compression score\n"
"\n"
"           For obvious reasons, the --scores and --weighted options are\n"
"           incompatible with this flag.\n"
"\n"
"           [1] Bo Li*, Nathanael Fillmore*, Yongsheng Bai, Mike Collins,\n"
"           James A. Thompson, Ron Stewart, Colin N. Dewey. Evaluation of\n"
"           de novo transcriptome assemblies from RNA-Seq data.\n"
"\n"
"Usage: Input and output specification\n"
"\n"
"   --A-seqs arg\n"
"\n"
"           The assembly sequences, in FASTA format. Required.\n"
"\n"
"   --B-seqs arg\n"
"\n"
"           The reference sequences, in FASTA format. Required.\n"
"\n"
"   --A-expr arg\n"
"\n"
"           The assembly expression, for use in weighted scores, as\n"
"           produced by RSEM in a file called *.isoforms.results.\n"
"           Required for weighted variants of scores.\n"
"\n"
"   --B-expr arg\n"
"\n"
"           The reference expression, for use in weighted scores, as\n"
"           produced by RSEM in a file called *.isoforms.results.\n"
"           Required for weighted variants of scores.\n"
"\n"
"   --A-to-B arg\n"
"\n"
"           The alignments of the assembly to the reference. The file\n"
"           format is specified by --alignment-type. Required for\n"
"           alignment-based scores.\n"
"\n"
"   --B-to-A arg\n"
"\n"
"           The alignments of the reference to the assembly. The file\n"
"           format is specified by --alignment-type. Required for\n"
"           alignment-based scores.\n"
"\n"
"   --alignment-type arg\n"
"\n"
"           The type of alignments used, either blast or psl. Default:\n"
"           psl. Currently BLAST support is experimental, not well\n"
"           tested, and not recommended.\n"
"\n"
"Usage: Options that modify the score definitions (and hence output)\n"
"\n"
"   --strand-specific\n"
"\n"
"           If this flag is present, it is assumed that all the assembly\n"
"           and reference sequences have the same orientation. Thus,\n"
"           alignments or kmer matches that are to the reverse strand are\n"
"           ignored.\n"
"\n"
"   --readlen arg\n"
"\n"
"           This option only applies to the KC scores. The read length of\n"
"           the reads used to build the assembly, used in the denominator\n"
"           of the ICR. Required for KC scores.\n"
"\n"
"   --num-reads arg\n"
"\n"
"           This option only applies to the KC scores. The number of\n"
"           reads used to build the assembly, used in the denominator of\n"
"           the ICR. Required for KC scores.\n"
"\n"
"   --kmerlen arg\n"
"\n"
"           This option only applies to the kmer and KC scores. This is\n"
"           the length (\"k\") of the kmers used in the definition of the\n"
"           KC and kmer scores. Required for KC and kmer scores.\n"
"\n"
"   --min-frac-identity arg\n"
"\n"
"           This option only applies to contig scores. Alignments with\n"
"           fraction identity less than this threshold are ignored. The\n"
"           fraction identity of an alignment is min(x/y, x/z), where\n"
"\n"
"              * $x$ is the number of bases in the assembly sequence that\n"
"                are aligned to an identical base in the reference\n"
"                sequence, according to the alignment,\n"
"              * $y$ is the number of bases in the assembly sequence, and\n"
"              * $z$ is the number of bases in the reference sequence.\n"
"\n"
"           Default: 0.99.\n"
"\n"
"   --max-frac-indel arg\n"
"\n"
"           This option only applies to contig scores. Alignments with\n"
"           fraction indel greater than this threshold are ignored. For\n"
"           psl alignments, the fraction indel of an alignment is\n"
"           $\\max(w/y, x/z)$, where\n"
"\n"
"              * $w$ is the number of bases that are inserted in the\n"
"                assembly sequence, according to the alignment (\"Q gap\n"
"                bases\"),\n"
"              * $x$ is the number of bases that are inserted in the\n"
"                reference sequence, according to the alignment (\"T gap\n"
"                bases\"),\n"
"              * $y$ is the number of bases in the assembly sequence, and\n"
"              * $z$ is the number of bases in the reference sequence.\n"
"\n"
"           For blast alignments, the fraction indel of an alignment is\n"
"           $\\max(x/y, x/z)$, where\n"
"\n"
"              * $x$ is the number of gaps bases that are inserted in the\n"
"                reference sequence, according to the alignment (\"gaps\"),\n"
"              * $y$ is the number of bases in the assembly sequence, and\n"
"              * $z$ is the number of bases in the reference sequence.\n"
"\n"
"           Default: 0.01.\n"
"\n"
"   --min-segment-len arg\n"
"\n"
"           This option only applies to nucleotide and pair scores.\n"
"           Alignment segments that contain fewer than this number of\n"
"           bases will be discarded. Default: 100. In the DETONATE paper,\n"
"           this was set to the read length.\n"
"\n"
"Usage: Options that modify the algorithm, but not the score definitions\n"
"\n"
"   --hash-table-type arg\n"
"\n"
"           The type of hash table to use, either \"sparse\" or \"dense\".\n"
"           This is only relevant for KC and kmer scores. The sparse\n"
"           table is slower but uses less memory. The dense table is\n"
"           faster but uses more memory. Default: \"sparse\".\n"
"\n"
"   --hash-table-numeric-type arg\n"
"\n"
"           The numeric type to use to store values in the hash table,\n"
"           either \"double\" or \"float\". This is only relevant for KC and\n"
"           kmer scores. Using single-precision floating point numbers\n"
"           (\"float\") requires less memory than using double-precision\n"
"           (\"double\"), but may also result in more numerical error. Note\n"
"           that we use double-precision numbers throughout our\n"
"           calculations even if single-precision numbers are stored in\n"
"           the table, so the additional error should be minimal.\n"
"           Default: \"double\".\n"
"\n"
"   --hash-table-fudge-factor arg\n"
"\n"
"           This is only relevant for KC and kmer scores. When the hash\n"
"           table is created, its initial capacity is set as the total\n"
"           worst-case number of possible kmers in the assembly and\n"
"           reference, based on each sequence's length, divided by the\n"
"           fudge factor. The default, 2.0, is often reasonable because\n"
"           (1) most kmers should be shared by the assembly and the\n"
"           reference, and (2) many kmers will be repeated several times.\n"
"           However, if you have a lot of memory or a really bad\n"
"           assembly, you could try a smaller number. Default: 2.0.\n"
"\n"
"Usage: Options to include additional output\n"
"\n"
"   --trace arg\n"
"\n"
"           If given, the prefix for additional output that provides\n"
"           details about the REF-EVAL scores; if not given, no such\n"
"           output is produced. Currently, the only such output is as\n"
"           follows.\n"
"\n"
"              * (--trace).{weighted,unweighted}_contig_{precision,recall}_matching\n"
"                is a TSV file that describes the matching used to\n"
"                compute the weighted or unweighted contig precision or\n"
"                recall. (Details about the matching are given in the\n"
"                section on score definitions below.) For recall, each\n"
"                row corresponds to a reference sequence $b$. Column 1\n"
"                contains $b$'s name. If $b$ is matched to a contig $a$,\n"
"                then the remaining columns are as follows:\n"
"\n"
"                   * Column 2 contains $a$'s name.\n"
"                   * Column 3 contains the weight of the edge between\n"
"                     $b$ and $a$. (This is set to the uniform weights\n"
"                     $1/|B|$ in the unweighted case, although the\n"
"                     maximum cardinality matching algorithm does not\n"
"                     actually use these weights.)\n"
"                   * Column 4 contains the names of all the contigs $a'$\n"
"                     that are adjacent to $b$ in the bipartite graph\n"
"                     that the matching is based on, separated by commas.\n"
"                     Thus, this column lists all the contigs $a'$ that\n"
"                     have a \"good enough\" match with the reference\n"
"                     sequence $b$, according to the criteria used to\n"
"                     build the bipartite graph. (See the section below\n"
"                     on score definitions for details.)\n"
"\n"
"                Otherwise, if $b$ is not matched to any contig, columns\n"
"                2 and 3 contain \"NA\". For precision, the file has the\n"
"                same format, but with the reference and the assembly\n"
"                interchanged. In other words, each row corresponds to a\n"
"                contig $a$ and contains information about its matching\n"
"                to a reference sequence $b$, or all \"NA\" values if $a$\n"
"                was not matched.\n"
"\n"
"Usage: General options\n"
"\n"
"   -? [ --help ]\n"
"\n"
"           Display this information.\n"
"\n"
"Score definitions\n"
"\n"
"   In the next few sections, we define the scores computed by REF-EVAL.\n"
"   Throughout, $A$ denotes the assembly, and $B$ denotes the reference.\n"
"   (As discussed under \"Usage\" above, the reference can be either an\n"
"   estimate of the \"true\" assembly or a collection of full-length\n"
"   reference transcripts.) Both $A$ and $B$ are thought of as sets of\n"
"   sequences. $A$ is a set of contigs, and $B$ is a set of reference\n"
"   sequences.\n"
"\n"
"Score definitions: contig precision, recall, and F1\n"
"\n"
"   The contig recall is defined as follows:\n"
"\n"
"     * Align the assembly $A$ to the reference $B$. Notation: each\n"
"       alignment $l$ is between a contig $a$ in $A$ and an reference\n"
"       sequence $b$ in $B$.\n"
"     * Throw out alignments that are to the reverse strand, if\n"
"       --strand-specific is present.\n"
"     * Throw out alignments whose fraction identity is less than\n"
"       --min-frac-identity (q.v.\\ for the definition of \"fraction\n"
"       identity\").\n"
"     * Throw out alignments whose fraction indel is greater than\n"
"       --max-frac-indel (q.v.\\ for the definition of \"fraction indel\").\n"
"     * Construct a bipartite graph from the remaining alignments, in\n"
"       which there is an edge between $a$ and $b$ iff there is a\n"
"       remaining alignment $l$ of $a$ to $b$.\n"
"     * If --weighted=yes, specify a weight for each edge between $a$ and\n"
"       $b$, namely $\\tau(b)$, the relative abundance of $b$ within the\n"
"       reference, as specified in --B-expr.\n"
"     * The unweighted contig recall is the number of edges in the\n"
"       maximum cardinality matching of this graph, divided by the number\n"
"       of sequences in the reference $B$.\n"
"     * The weighted contig recall is the weight of the maximum weight\n"
"       matching of this graph.\n"
"\n"
"   The contig precision is defined as follows: Interchange the assembly\n"
"   and the reference, and compute the contig recall.\n"
"\n"
"   The contig F1 is the harmonic mean of the precision and recall.\n"
"\n"
"Score definitions: nucleotide precision, recall, and F1\n"
"\n"
"   The nucleotide recall is defined as follows:\n"
"\n"
"     * Align the assembly $A$ to the reference $B$. Notation: each\n"
"       alignment $l$ is between a contig $a \\in A$ and an reference\n"
"       element $b \\in B$.\n"
"     * Throw out alignments that are to the reverse strand, if\n"
"       --strand-specific is present.\n"
"     * Throw out alignments that are shorter than --min-fragment-length.\n"
"     * Add each remaining alignment to a priority queue, with priority\n"
"       equal to the number of identical bases in the alignment.\n"
"     * Let numer = 0.\n"
"     * While the priority queue is not empty:\n"
"\n"
"          * Pop the alignment $l$ with highest priority.\n"
"          * Add the number of identical bases in the alignment to numer.\n"
"          * Subtract $l$ from all the other alignments in the queue and\n"
"            update their priorities (see below).\n"
"\n"
"     * Let denom be the total number of bases in the reference $B$.\n"
"     * The unweighted nucleotide recall is numer/denom.\n"
"\n"
"   The actual implementation uses a more complicated and efficient\n"
"   algorithm than the one above.\n"
"\n"
"   If --weighted=yes, then (i) \"the number of identical bases\" above is\n"
"   replaced by \"the number of identical bases, times $\\tau(b)$\", in the\n"
"   definition of the priority and the numer, and (ii) \"total number of\n"
"   bases in the reference $B$\" is replaced by \"$\\sum_{b \\in B} \\tau(b)\n"
"   length(b)$\". In other words, each base (of a reference sequence),\n"
"   throughout the computation, is weighted by the expression level of\n"
"   its parent sequence.\n"
"\n"
"   The nucleotide precision is defined as follows: Interchange the\n"
"   assembly and the reference, and compute the nucleotide recall.\n"
"\n"
"   The nucleotide F1 is the harmonic mean of the precision and recall.\n"
"\n"
"   Alignment subtraction is defined as follows.\n"
"\n"
"     * An alignment $l$ from $a$ to $b$ can be thought of as a set of\n"
"       pairs of disjoint intervals\n"
"       $$ \\{ ([s_1(a), e_1(a)], [s_1(b), e_1(b)]), \\dots, ([s_n(a),\n"
"       e_n(a)], [s_n(b), e_n(b)]) \\}, $$\n"
"       where each pair $([s_i(a), e_i(a)], [s_i(b), e_i(b)])$\n"
"       corresponds to an ungapped segment of the alignment: $s_i(a)$ and\n"
"       $e_i(a)$ are the segment's start and end positions within a, and\n"
"       $s_i(b)$ and $e_i(b)$ are the segment's start and end positions\n"
"       within b. In the case of non-strand-specific alignments, $s_i(b)$\n"
"       might be greater than $e_i(b)$.\n"
"     * If $l$ is an alignment from $a$ to $b$, $l'$ is an alignment from\n"
"       $a'$ to $b'$, $a \\neq a'$, and $b \\neq b'$, then the difference\n"
"       $l - l' = l$.\n"
"     * If $l$ is an alignment from $a$ to $b$, $l'$ is an alignment from\n"
"       $a'$ to $b'$, $a = a'$, and $b \\neq b'$, then the difference $l -\n"
"       l' = l''$, defined as follows. Each alignment segment of $l$ is\n"
"       compared to the alignment segments of $l'$. If a segment of $l$\n"
"       overlaps one of the segments of $l'$ wrt $a$, it is truncated so\n"
"       as to avoid the overlap. This truncation may result in zero, one,\n"
"       or two replacement alignment segments. (If the overlapping\n"
"       alignment segment of $l'$ is contained strictly within the\n"
"       segment of $l$, wrt $a$, two segments will result.)\n"
"     * If $l$ is an alignment from $a$ to $b$, $l'$ is an alignment from\n"
"       $a'$ to $b'$, $a \\neq a'$, and $b = b'$, then the difference $l -\n"
"       l' = l''$, defined similarly as in the previous item, except\n"
"       overlaps are examined and resolved wrt $b$.\n"
"\n"
"   A couple of examples of the above are as follows. The comments in\n"
"   [5]test_re_matched.cpp contain even more examples.\n"
"\n"
"   As a first example, consider alignments of an assembly $A = \\{a_0,\n"
"   a_1, a_2\\}$ to a reference $B = \\{b_0, b_1\\}$. In the pictures below,\n"
"   each alignment segment is indicated by a pair diagonal or vertical\n"
"   lines, with its name (initially $x$, $y$, $z$, $w$) in between the\n"
"   two lines.\n"
"\n"
"       b0                    b1\n"
"    B  -----------------     -------------\n"
"            /    \\   /  \\   /  | / |\n"
"           /      \\ /    \\ /   |/  |\n"
"          /  x     /  y   \\  z / w |\n"
"         /        / \\    / \\  /|   |\n"
"    A    ---------   --------- ---------\n"
"         a0          a1        a2\n"
"\n"
"   Assume:\n"
"\n"
"     * $x > y > z > w$, where $>$ compares alignment size measured by\n"
"       the number of identical bases.\n"
"     * $y - x < z$.\n"
"     * $y - x$ is contained in $z$, wrt $A$.\n"
"\n"
"   Step 1: Process alignment $x$, resulting in\n"
"\n"
"       b0                    b1\n"
"    B  ------------------   --------------\n"
"            /        /\\ \\   /  | / |\n"
"           /        /  \\ \\ /   |/  |\n"
"          /  x     /    \\*\\  z / w |        * = y - x\n"
"         /        /      / \\  /|   |\n"
"    A    ---------   --------- ---------\n"
"         a0          a1        a2\n"
"\n"
"   Step 2: Process alignment $z$, resulting in\n"
"\n"
"       b0                    b1\n"
"    B  -----------------    --------------\n"
"            /        /      /    /||\n"
"           /        /      /    / ||\n"
"          /  x     /      /  z /  *|        * = w - z\n"
"         /        /      /    /   ||\n"
"    A    ---------   --------- ---------\n"
"         a0          a1        a2\n"
"\n"
"   Now we have a 1-1 mathing. The intervals of $B$ used to compute\n"
"   recall are as follows:\n"
"\n"
"       b0                    b1\n"
"    B  -----[--------]--    [----][]------\n"
"            /        /      /    /||\n"
"           ...\n"
"\n"
"   As a second example, we start with the same initial set of\n"
"   alignments:\n"
"\n"
"       b0                    b1\n"
"    B  -----------------     -------------\n"
"            /    \\   /  \\   /  | / |\n"
"           /      \\ /    \\ /   |/  |\n"
"          /  x     /  y   \\  z / w |\n"
"         /        / \\    / \\  /|   |\n"
"    A    ---------   --------- ---------\n"
"         a0          a1        a2\n"
"\n"
"   But we make a slightly different set of assumptions:\n"
"\n"
"     * $x > y > w > z$ ($w$ and $z$ are interchanged, compared to the\n"
"       first example).\n"
"     * $y - x < w$.\n"
"     * $y - x > z - w$.\n"
"     * $y - x$ is contained in $z$, wrt $A$.\n"
"\n"
"   Step 1: Process alignment $x$, resulting in\n"
"\n"
"       b0                    b1\n"
"    B  ------------------   --------------\n"
"            /        /\\ \\   /  | / |\n"
"           /        /  \\ \\ /   |/  |\n"
"          /  x     /    \\*\\  z / w |        * = y - x\n"
"         /        /      / \\  /|   |\n"
"    A    ---------   --------- ---------\n"
"         a0          a1        a2\n"
"\n"
"   Step 2: Process alignment w, resulting in\n"
"\n"
"       b0                    b1\n"
"    B  ------------------   --------------\n"
"            /        /\\ \\   /  |   |\n"
"           /        /  \\ \\ /  /|   |\n"
"          /  x     /    \\*\\ +/ | w |        * = y - x\n"
"         /        /      / \\/  |   |        + = z - w\n"
"    A    ---------   --------- ---------\n"
"         a0          a1        a2\n"
"\n"
"   Step 3: Process alignment $y - x$, resulting in\n"
"\n"
"       b0                    b1\n"
"    B  ------------------   --------------\n"
"            /        /\\ \\  +/  |   |\n"
"           /        /  \\*\\//   |   |\n"
"          /  x     /    \\ \\    | w |        * = y - x\n"
"         /        /     // \\   |   |        + = (z - w) - (y - x)\n"
"    A    ---------   --------- ---------\n"
"         a0          a1        a2\n"
"\n"
"   Now we have a 1-1 mathing. The intervals of $B$ used to compute\n"
"   recall are:\n"
"\n"
"       b0                    b1\n"
"    B  -----[--------][-]   []-[---]------\n"
"            /        /\\ \\  //  |   |        * = y - x\n"
"             x         *   +     w          + = (z - w) - (y - x)\n"
"           ...\n"
"\n"
"Score definitions: pair precision, recall, and F1\n"
"\n"
"   The definitions for pair precision, recall, and F1 are exactly the\n"
"   same as for nucleotide precision, recall, and F1, except that instead\n"
"   of bases, we operate on pairs of bases.\n"
"\n"
"   For example, consider the reference $B$ (with one transcript) and\n"
"   assembly $A$ (with two contigs), where horizontal position indicates\n"
"   alignment.\n"
"\n"
"       B: b   = AGCTCGACGT\n"
"       A: a_1 = AGCT\n"
"          a_2 =     CGACGT\n"
"\n"
"   Here, the transcript recall is 0 (because neither $a_1$ nor $a_2$\n"
"   covers $b$ to $\\geq$ 99 percent), but the nucleotide recall is 1\n"
"   (because $a_1$ and $a_2$ jointly cover $b$ completely). The pair\n"
"   recall is somewhere in between, because the following pairs of $b$\n"
"   are correctly predicted (represented by an upper triangular indicator\n"
"   matrix):\n"
"\n"
"           First base\n"
"       S   AGCTCGACGT\n"
"       e A 1111\n"
"       c G  111\n"
"       o C   11\n"
"       n T    1\n"
"       d C     111111\n"
"         G      11111\n"
"       b A       1111\n"
"       a C        111\n"
"       s G         11\n"
"       e T          1\n"
"\n"
"Score definitions: KC and related scores\n"
"\n"
"   The kmer compression score (KC score) is a combination of two\n"
"   measures, weighted kmer recall (WKR) and inverse compression rate\n"
"   (ICR), and is simply\n"
"\n"
"   $$ KC = WKR - ICR. $$\n"
"\n"
"   The WKR measures the fidelity with which a particular assembly\n"
"   represents the kmer content of the reference sequences. Balancing the\n"
"   WKR, the ICR measures the degree to which the assembly compresses the\n"
"   RNA-Seq data. The details of the WKR and ICR measures are provided\n"
"   below.\n"
"\n"
"   To compute the WKR, the relative abundances of the reference elements\n"
"   are required, as specified by --B-expr. Given the reference sequences\n"
"   and their abundances, a kmer occurrence frequency profile, $p$, is\n"
"   computed, with individual kmer occurrences weighted by their parent\n"
"   sequences' abundances: for each kmer $r$, we define\n"
"\n"
"   $$ p(r) = \\frac{ \\sum_{b \\in B} n(r,b) \\tau(b) } { \\sum_{b \\in B}\n"
"   n(b) \\tau(b) } $$\n"
"\n"
"   where $B$ is the set of reference sequences, and for each reference\n"
"   sequence $b \\in B$:\n"
"\n"
"     * $n(r,b)$ is the number of times the kmer $r$ occurs in $b$,\n"
"     * $n(b) $ is the total number of kmers in $b$, and\n"
"     * $\\tau(b)$ is the relative abundance of $b$.\n"
"\n"
"   Letting $R(A)$ be the set of all kmers in the assembly $A$, the\n"
"   weighted kmer recall (WKR) is defined as\n"
"\n"
"   $$ WKR = \\sum_{r \\in R(A)} p(r). $$\n"
"\n"
"   REF-EVAL currently uses --readlen as the kmer length.\n"
"\n"
"   Since recall measures only tell half of the story regarding accuracy,\n"
"   the KC score includes a second term, the ICR, which serves to\n"
"   penalize large assemblies. We define the inverse compression rate\n"
"   (ICR) of an assembly as\n"
"\n"
"   $$ ICR = n_A/(N L), $$\n"
"\n"
"   where\n"
"\n"
"     * $n_A$ is the total number of bases in the assembly $A$,\n"
"     * $N$ is the total number of reads, as specified by --num-reads,\n"
"       and\n"
"     * $L$ is the read length, as specified by --readlen.\n"
"\n"
"Score definitions: kmer scores\n"
"\n"
"   If --weighted=yes, we construct a kmer occurrence frequency profile\n"
"   $p_B$ for $B$ exactly as described in the previous section (about the\n"
"   KC score). We construct a kmer occurrence frequency profile $p_A$ for\n"
"   $A$ similarly. The relative abundances are specified by --A-expr and\n"
"   --B-expr.\n"
"\n"
"   If --weighted=no, we construct the kmer occurrence frequency profiles\n"
"   $p_A$ and $p_B$ in the same way, except that uniform relative\n"
"   abundances are used, i.e., $\\tau(a) = 1/|A|$ for all $a$ in $A$, and\n"
"   $\\tau(b) = 1/|B|$ for all $b$ in $B$, where $|A|$ is the number of\n"
"   contigs in $A$, and $|B|$ is the number of reference sequences in\n"
"   $B$.\n"
"\n"
"   Let $m$ be the \"mean\" profile of $p_A$ and $p_B$:\n"
"\n"
"   $$ m(r) = (1/2) (p_A(r) + p_B(r)) \\qquad\\hbox{for every kmer $r$}. $$\n"
"\n"
"   The Jensen-Shannon divergence between $p_A$ and $p_B$ is defined in\n"
"   terms of the KL divergence between $p_A$ and the mean, and $p_B$ and\n"
"   the mean, as follows:\n"
"\n"
"     * Let $KL(p_A || m) = \\sum_r p_A(r) (\\log_2(p_A(r)) -\n"
"       \\log_2(m(r)))$.\n"
"     * Let $KL(p_B || m) = \\sum_r p_B(r) (\\log_2(p_B(r)) -\n"
"       \\log_2(m(r)))$.\n"
"     * Let $JS(p_A || p_B) = (1/2) (KL(p_A || m) + KL(p_B || m))$.\n"
"\n"
"   In the output file, these three scores are denoted\n"
"   (un)weighted_kmer_KL_A_to_M, (un)weighted_kmer_KL_B_to_M, and\n"
"   (un)weighted_kmer_jensen_shannon, respectively.\n"
"\n"
"   The Hellinger distance between $p_A$ and $p_B$ is defined as\n"
"\n"
"   $$ \\sqrt{ (1/2) \\sum_r (\\sqrt{p_A(r)} - \\sqrt{p_B(r)})^2 } $$\n"
"\n"
"   The total variation distance between $p_A$ and $p_B$ is defined as\n"
"\n"
"   $$ (1/2) \\sum_r |p_A(r) - p_B(r)|, $$\n"
"\n"
"   where $|\\cdot|$ denotes absolute value. Above, $\\sum_r$ denotes a sum\n"
"   over all possible kmers $r$ (most of which will have $p_A(r) = p_B(r)\n"
"   = 0$).\n"
"\n"
"References\n"
"\n"
"   Visible links\n"
"   1. http://deweylab.biostat.wisc.edu/detonate/ref-eval-estimate-true-assembly.html\n"
"   2. http://genome.ucsc.edu/FAQ/FAQblat.html\n"
"   3. http://deweylab.biostat.wisc.edu/rsem/\n"
"   4. http://deweylab.biostat.wisc.edu/detonate/rsem-eval.html\n"
"   5. https://github.com/nfillmore/ref-eval/blob/master/test_re_matched.cpp\n"
; }
