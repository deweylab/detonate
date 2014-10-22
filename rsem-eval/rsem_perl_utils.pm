#!/usr/bin/env perl

package rsem_perl_utils;

use strict;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT = qw(runCommand);
our @EXPORT_OK = qw(runCommand collectResults showVersionInfo nbinom_calc_moments nbinom_load_params nbinom_convert_params);

# command, {err_msg}
sub runCommand {
    print $_[0]."\n";
    my $status = system($_[0]);

    if ($? == -1) {
	my @arr = split(/[ \t]+/, $_[0]);
	print "$arr[0] : $!!\n";
	print "Please check if you have compiled the associated codes by typing related \"make\" commands and/or made related executables ready to use.\n";
	exit(-1);
    }

    if ($status != 0) {
        my $errmsg = "";
        if (scalar(@_) > 1) { $errmsg .= $_[1]."\n"; }
	$errmsg .= "\"$_[0]\" failed! Plase check if you provide correct parameters/options for the pipeline!\n";
	print $errmsg;
        exit(-1);
    }
    print "\n";
}

my @allele_title = ("allele_id", "transcript_id", "gene_id", "length", "effective_length", "expected_count", "TPM", "FPKM", "AlleleIsoPct", "AlleleGenePct", "posterior_mean_count", "posterior_standard_deviation_of_count", "pme_TPM", "pme_FPKM", "AlleleIsoPct_from_pme_TPM", "AlleleGenePct_from_pme_TPM", "TPM_ci_lower_bound", "TPM_ci_upper_bound", "FPKM_ci_lower_bound", "FPKM_ci_upper_bound");

my @transcript_title = ("transcript_id", "gene_id", "length", "effective_length", "expected_count", "TPM", "FPKM", "IsoPct", "posterior_mean_count", "posterior_standard_deviation_of_count", "pme_TPM", "pme_FPKM", "IsoPct_from_pme_TPM", "TPM_ci_lower_bound", "TPM_ci_upper_bound", "FPKM_ci_lower_bound", "FPKM_ci_upper_bound");

my @gene_title = ("gene_id", "transcript_id(s)", "length", "effective_length", "expected_count", "TPM", "FPKM", "posterior_mean_count", "posterior_standard_deviation_of_count", "pme_TPM", "pme_FPKM", "TPM_ci_lower_bound", "TPM_ci_upper_bound", "FPKM_ci_lower_bound", "FPKM_ci_upper_bound");

# type, inpF, outF
sub collectResults {
    my $local_status;
    my ($inpF, $outF);
    my @results = ();
    my $line;

    $inpF = $_[1];
    $outF = $_[2];

    $local_status = open(INPUT, $inpF);
    if ($local_status == 0) { print "Fail to open file $inpF!\n"; exit(-1); }
    
    @results = ();
    
    while ($line = <INPUT>) {
	chomp($line);
	my @local_arr = split(/\t/, $line);
	push(@results, \@local_arr); 
    }

    close(INPUT);

    $local_status = open(OUTPUT, ">$outF");
    if ($local_status == 0) { print "Fail to create file $outF!\n"; exit(-1); }

    my $n = scalar(@results);
    my $m = scalar(@{$results[0]});

    $" = "\t";

    my @out_arr = ();
    for (my $i = 0; $i < $n; $i++) {
	if ($_[0] eq "allele") { push(@out_arr, $allele_title[$i]); }
	elsif ($_[0] eq "isoform") { push(@out_arr, $transcript_title[$i]); }
	elsif ($_[0] eq "gene") { push(@out_arr, $gene_title[$i]); }
	else { print "A bug on 'collectResults' is detected!\n"; exit(-1); }
    }
    print OUTPUT "@out_arr\n";

    for (my $i = 0; $i < $m; $i++) {
	@out_arr = ();
	for (my $j = 0; $j < $n; $j++) { push(@out_arr, $results[$j][$i]); }
	print OUTPUT "@out_arr\n"; 
    }

    close(OUTPUT);
}

# 0, dir
sub showVersionInfo {
    open(INPUT, "$_[0]/WHAT_IS_NEW");
    my $line = <INPUT>;
    chomp($line);
    close(INPUT);
    print "Current version is $line\n";
    exit(0);
}

# 0, input multi-FASTA file
# RETURN, mean and sd
sub nbinom_calc_moments {
    my $doesOpen = open(INPUT, $_[0]);
    if ($doesOpen == 0) { print "Cannot open $_[0]!\n"; exit(-1); }
    my ($line, $seq) = ();
    my @lens = ();
    $line = <INPUT>;
    chomp($line);
    while (substr($line, 0, 1) eq '>') {
	$seq = "";
	my $len = 0;
	while (($seq = <INPUT>) && (substr($seq, 0, 1) ne '>')) {
	    chomp($seq);
	    $len += length($seq);
	}
	push(@lens, $len);
	$line = $seq;
    }
    close(INPUT);

    my ($mu, $mu2, $sigma2) = (0, 0, 0);

    my $n = @lens;
    for (my $i = 0; $i < $n; $i++) {
	$mu += $lens[$i];
	$mu2 += $lens[$i] * $lens[$i];
    }
    $mu /= $n;
    $mu2 /= $n;
    $sigma2 = $mu2 - $mu * $mu;

    return ($mu, $sigma2 ** 0.5);
}

# 0, input_file
# RETURN: mean and sd
sub nbinom_load_params {
    my $doesOpen = open(INPUT, $_[0]);
    if ($doesOpen == 0) { print "Cannot open $_[0]!\n"; exit(-1); }
    my $line = <INPUT>;
    chomp($line);
    my @arr = split(/\t/, $line);
    close(INPUT);

    return @arr;
}

# 0, mean; 1, sd
# RETURN: r and p
sub nbinom_convert_params {
    my ($mu, $var) = ($_[0], $_[1] ** 2);
    my @res = ($mu ** 2 / ($var - $mu), $mu / $var);
    return @res;
}

1;
