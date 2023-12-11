package Bio::Adventure::SeqMisc;
use Modern::Perl;
use autodie qw":all";
use diagnostics;
use warnings qw"all";

use File::Temp qw / tmpnam /;
use parent 'Exporter';
our @EXPORT = qw"$references";
our @EXPORT_OK = qw"$references";
## Every function here should take as input an array reference
## containing the sequence to be shuffled.

=head1 NAME

Bio::Adventure::SeqMisc - Given a primary sequence, collect some information.

=head1 SYNOPSIS

The information collected includes the following:

amino acids in all reading frames, nucleotide/aa/dinucleotide/codon
frequencies.  Weights (in Daltons), pI of amino acids, volume of
amino acids, polarity of amino acids, hydrophobicity, solubility,
pYrimidine/puRine ratio, # hydrogen bonds by codon, charge by
amino acid, pkCOOH, pkNHHH, pkR, occurrence value, accessibility
value, %amino acid buried, isoelectric point, amino acid sequences
in all reading frames, dinucleotide sequences of the 12, 23, 31,
and 13 transitions.  Reverse sequence, complementary sequence,
revcomp.  GC/CT content.

Conversely, perform randomizations using a few algorithms
including via the external squid library.

=head1 DATA STRUCTURES

=over

=item C<$references>

The namespace global variable '$references' contains some generic
information.  Note to self: I've been meaning for years to include
either the codon tables from the codon database or a way to read
in a text file.  Do that, damnit!

=back

=cut
our $references = {
    aminos => {
        A => [['GCA', 'GCC', 'GCG', 'GCU'], ['Ala', 'Alanine']],
        C => [['UGC', 'UGU'], ['Cys', 'Cysteine']],
        D => [['GAC', 'GAU'], ['Asp', 'Aspartic Acid']],
        E => [['GAA', 'GAG'], ['Glu', 'Glutamic Acid']],
        F => [['UUC', 'UUU'], ['Phe', 'Phenylalanine']],
        G => [['GGA', 'GGC', 'GGG', 'GGU'], ['Gly', 'Glycine']],
        H => [['CAC', 'CAU'], ['His', 'Histidine']],
        I => [['AUA', 'AUC', 'AUU'], ['Ile', 'Isoleucine']],
        K => [['AAA', 'AAG'], ['Lys', 'Lysine']],
        L => [['UUA', 'UUG', 'CUA', 'CUC', 'CUG', 'CUU'], ['Leu', 'Leucine']],
        M => [['AUG'], ['Met', 'Methionine']],
        N => [['AAC', 'AAU'], ['Asn', 'Asparagine']],
        P => [['CCA', 'CCC', 'CCG', 'CCU'], ['Pro', 'Proline']],
        Q => [['CAA', 'CAG'], ['Gln', 'Glutamine']],
        R => [['CGA', 'CGC', 'CGG', 'CGU', 'AGG', 'AGA'], ['Arg', 'Arginine']],
        S => [['UCA', 'UCC', 'UCG', 'UCU', 'ACG', 'AGU', 'AGC'], ['Ser', 'Serine']],
        T => [['ACA', 'ACC', 'ACG', 'ACU'], ['Thr', 'Threonine']],
        V => [['GUA', 'GUC', 'GUG', 'GUU'], ['Val', 'Valine']],
        W => [['UGG'], ['Trp', 'Tryptophan']],
        Y => [[ 'UAC', 'UAU'], ['Tyr', 'Tyrosine']],
        '*' => [['UAA', 'UGA', 'UAG'], ['Ter', 'Termination']],
    },
    codons => {
        GCA => 'A', GCC => 'A', GCG => 'A', GCU => 'A',
        UGC => 'C', UGU => 'C',
        GAC => 'D', GAU => 'D',
        GAA => 'E', GAG => 'E',
        UUC => 'F', UUU => 'F',
        GGA => 'G', GGC => 'G', GGG => 'G', GGU => 'G',
        CAC => 'H', CAU => 'H',
        AUA => 'I', AUC => 'I', AUU => 'I',
        AAA => 'K', AAG => 'K',
        UUA => 'L', UUG => 'L', CUA => 'L', CUC => 'L', CUG => 'L', CUU => 'L',
        AUG => 'M',
        AAC => 'N', AAU => 'N',
        CCA => 'P', CCC => 'P', CCG => 'P', CCU => 'P',
        CAA => 'Q', CAG => 'Q',
        CGA => 'R', CGC => 'R', CGG => 'R', CGU => 'R', AGG => 'R', AGA => 'R',
        UCA => 'S', UCC => 'S', UCG => 'S', UCU => 'S', ACG => 'S', AGU => 'S', AGC => 'S',
        ACA => 'T', ACC => 'T', ACG => 'T', ACU => 'T',
        GUA => 'V', GUC => 'V', GUG => 'V', GUU => 'V',
        UGG => 'W',
        UAC => 'Y', UAU => 'Y',
        UAA => '*', UGA => '*', UAG => '*',
        ## A catchall to pick up Ns
        '?' => [['NNN',],['Hmm', 'Dunno']],
        NNN => '?',
        NNA => '?', NNU => '?', NNG => '?', NNC => '?',
        NAA => '?', NAU => '?', NAG => '?', NAC => '?',
        NUA => '?', NUU => '?', NUG => '?', NUC => '?',
        NCA => '?', NCU => '?', NCG => '?', NCC => '?',
        NGA => '?', NGU => '?', NGG => '?', NGC => '?',
        ANN => '?', UNN => '?', GNN => '?', CNN => '?',
        AAN => '?', AUN => '?', AGN => '?', ACN => '?',
        UAN => '?', UUN => '?', UGN => '?', UCN => '?',
        GAN => '?', GUN => '?', GCN => '?', GGN => '?',
        CAN => '?', CUN => '?', CCN => '?', CCN => '?',
        ANA => '?', ANU => '?', ANC => '?', ANG => '?',
        UNA => '?', UNU => '?', UNC => '?', UNG => '?',
        CNA => '?', CNU => '?', CNC => '?', CNG => '?',
        GNA => '?', GNU => '?', GNC => '?', GNG => '?',
    },
    ## Weights of amino acids taken from sigma's catalog
    weights => {
        A => 89.1,
        R => 174.2,
        N => 132.12,
        D => 133.11,
        C => 121.16,
        E => 147.13,
        Q => 146.15,
        G => 75.07,
        H => 155.16,
        O => 131.13,
        I => 131.18,
        L => 131.18,
        K => 146.19,
        M => 149.19,
        F => 165.19,
        P => 115.13,
        U => 139.11,
        S => 105.09,
        T => 119.12,
        W => 204.23,
        Y => 181.19,
        V => 117.15,
    },
    pI => {
        A => 6.0,
        R => 10.76,
        N => 5.41,
        D => 2.77,
        C => 5.07,
        E => 3.22,
        Q => 5.65,
        G => 5.97,
        H => 7.59,
        O => undef,
        I => 6.02,
        L => 5.98,
        K => 9.74,
        M => 5.74,
        F => 5.48,
        P => 6.30,
        U => 5.68,
        S => 5.68,
        T => 5.60,
        W => 5.89,
        Y => 5.66,
        V => 5.96,
    },
    ## Volume metrics taken from Darby and Creighton 1993
    volume => {
        A => 67,
        R => 148,
        N => 96,
        D => 91,
        C => 86,
        E => 109,
        Q => 114,
        G => 48,
        H => 118,
        I => 124,
        L => 124,
        K => 135,
        M => 124,
        F => 135,
        P => 90,
        S => 73,
        T => 163,
        Y => 141,
        V => 105,
        W => 163,
    },
    polarity_rank => {
        A => 9,
        R => 15,
        N => 16,
        D => 19,
        C => 7,
        E => 18,
        Q => 17,
        G => 11,
        H => 10,
        I => 1,
        L => 3,
        K => 20,
        M => 5,
        F => 2,
        P => 13,
        S => 14,
        T => 12,
        Y => 8,
        V => 4,
        W => 6,
    },
    ## Taken from sigma
    hydrophobicity => {
        G => 0,
        F => 100,
        I => 99,
        W => 97,
        L => 97,
        V => 76,
        M => 74,
        Y => 63,
        C => 49,
        A => 41,
        T => 13,
        H => 8,
        S => -5,
        Q => -10,
        R => -14,
        K => -23,
        N => -28,
        E => -31,
        P => -46,
        D => -55,
    },
    solubility => {
        A => 16.65,
        R => 15.0,
        D => 0.778,
        N => 3.53,
        C => 200,  ## "very high"
        E => 0.864,
        Q => 2.5,
        G => 24.99,
        H => 4.19,
        I => 4.117,
        L => 2.426,
        K => 200, ## "Very high"
        M => 3.381,
        F => 2.965,
        P => 162.3,
        S => 5.023,
        T => 200, ## Very high
        W => 1.136,
        Y => 0.0453,
        V => 8.85,
    },
    YR => {
        CCC => 000, CCU => 000, ## P
        UCC => 000, UCU => 000, ## S
        CUC => 000, CUU => 000, ## L
        UUC => 000, UUU => 000, ## F
        CCA => 001, CCG => 001, ## P
        UCA => 001, UCG => 001, ## S
        CUA => 001, CUG => 001, ## L
        UUA => 001, UUG => 001, ## L
        GCC => 100, GCU => 100, ## A
        ACC => 100, ACU => 100, ## T
        GUC => 100, GUU => 100, ## V
        AUC => 100, AUU => 100, ## I
        GCA => 101, GCG => 101, ## A
        ACA => 101, ACG => 101, ## T
        GUA => 101, GUG => 101, ## V
        AUA => 101, AUG => 101, ## I/M
        CGC => 010, CGU => 010, ## R
        UGC => 010, UGU => 010, ## C
        CAC => 010, CAU => 010, ## H
        UAC => 010, UAU => 010, ## Y
        CGA => 011, CGG => 011, ## R
        UGA => 011, UGG => 011, ## */W
        CAA => 011, CAG => 011, ## Q
        UAA => 011, UAG => 011, ## *
        GGC => 110, GGU => 110, ## G
        AGC => 110, AGU => 110, ## S
        GAC => 110, GAU => 110, ## D
        AAC => 110, AAU => 110, ## N
        GGA => 111, GGG => 111, ## G
        AGA => 111, AGG => 111, ## R
        GAA => 111, GAG => 111, ## E
        AAA => 111, AAG => 111, ## K
    },
    hbonds => {
        CCC => 6, CCU => 6, ## P
        UCC => 5, UCU => 5, ## S
        CUC => 5, CUU => 5, ## L
        UUC => 4, UUU => 4, ## F
        CCA => 6, CCG => 6, ## P
        UCA => 5, UCG => 5, ## S
        CUA => 5, CUG => 5, ## L
        UUA => 4, UUG => 4, ## L
        GCC => 6, GCU => 6, ## A
        ACC => 5, ACU => 5, ## T
        GUC => 5, GUU => 5, ## V
        AUC => 4, AUU => 4, ## I
        GCA => 6, GCG => 6, ## A
        ACA => 5, ACG => 5, ## T
        GUA => 5, GUG => 5, ## V
        AUA => 4, AUG => 4, ## I/M
        CGC => 6, CGU => 6, ## R
        UGC => 5, UGU => 5, ## C
        CAC => 5, CAU => 5, ## H
        UAC => 4, UAU => 4, ## Y
        CGA => 6, CGG => 6, ## R
        UGA => 5, UGG => 5, ## */W
        CAA => 5, CAG => 5, ## Q
        UAA => 4, UAG => 4, ## *
        GGC => 6, GGU => 6, ## G
        AGC => 5, AGU => 5, ## S
        GAC => 5, GAU => 5, ## D
        AAC => 4, AAU => 4, ## N
        GGA => 6, GGG => 6, ## G
        AGA => 5, AGG => 5, ## R
        GAA => 5, GAG => 5, ## E
        AAA => 4, AAG => 4, ## K
    },
    charge => {
        D => -1.0, E => -1.0,
        K => 1.0, R => 1.0, H => 1.0,
        G => 0, A => 0, V => 0, L => 0, I => 0, M => 0, F => 0,
        W => 0, P => 0, S => 0, T => 0, C => 0, Y => 0, N => 0, Q => 0,
    },
    pKCOOH => {
        A => 2.34,
        R => 2.17,
        N => 2.02,
        D => 1.88,
        C => 1.96,
        E => 2.19,
        Q => 2.17,
        G => 2.34,
        H => 1.82,
        I => 2.36,
        L => 2.36,
        K => 2.18,
        M => 2.28,
        F => 1.83,
        P => 1.99,
        S => 2.21,
        T => 2.09,
        W => 2.83,
        Y => 2.20,
        V => 2.32,
    },
    pKNHHH => {
        A => 9.69,
        R => 9.04,
        N => 8.80,
        D => 9.60,
        C => 10.128,
        E => 9.67,
        Q => 9.13,
        G => 9.60,
        H => 9.17,
        I => 9.60,
        L => 9.60,
        K => 8.95,
        M => 9.21,
        F => 9.13,
        P => 10.60,
        S => 9.15,
        T => 9.10,
        W => 9.39,
        Y => 9.11,
        V => 9.62,
    },
    pKR => {
        A => undef,
        R => 12.48,
        N => undef,
        D => 3.65,
        C => 8.18,
        E => 4.25,
        Q => undef,
        G => undef,
        H => 6.0,
        I => undef,
        L => undef,
        K => 10.53,
        M => undef,
        F => undef,
        P => undef,
        S => undef,
        T => undef,
        W => undef,
        Y => 10.07,
        V => undef,
    },
    occurrence => {
        A => 7.5,
        R => 5.2,
        N => 4.6,
        D => 5.2,
        C => 1.8,
        E => 6.3,
        Q => 4.1,
        G => 7.1,
        H => 2.2,
        I => 5.5,
        L => 9.1,
        K => 5.8,
        M => 2.8,
        F => 3.9,
        P => 5.1,
        S => 7.4,
        T => 6.0,
        W => 1.3,
        Y => 3.3,
        V => 6.5,
    },
    accessible => {
        A => 67,
        R => 196,
        N => 113,
        D => 106,
        C => 104,
        E => 138,
        Q => 144,
        G => 0,
        H => 151,
        I => 140,
        L => 137,
        K => 167,
        M => 160,
        F => 175,
        P => 105,
        S => 80,
        T => 102,
        W => 217,
        Y => 187,
        V => 117,
    },
    buried => {
        A => 38,
        R => 0,
        N => 10,
        D => 14.5,
        C => 47,
        E => 20,
        Q => 6.3,
        G => 37,
        H => 19,
        I => 65,
        L => 14,
        K => 4.2,
        M => 50,
        F => 48,
        P => 24,
        S => 24,
        T => 25,
        W => 23,
        Y => 13,
        V => 56,
    },
    isoelectric => {
        A => 6.01,
        R => 10.76, ## +1 charge
        N => 5.41,
        D => 2.77,  ## -1 charge
        C => 5.07,
        E => 3.22,  ## -1 charge
        Q => 5.65,
        G => 5.97,
        H => 7.59,  ## +1 charge on 10%
        I => 6.02,
        L => 5.98,
        K => 9.74,  ## +1 charge
        M => 5.74,
        F => 5.48,
        P => 6.48,
        S => 5.68,
        T => 5.87,
        W => 5.88,
        Y => 5.66,
        V => 5.97,
    },
};

=head1 METHODS

=head2 C<new>

Construct a new SeqMisc object.  In so doing, perform the various queries suggested by the
global data structure above...

=cut
sub new {
    my ($class, %arg) = @_;
#    $arg{sequence} = [] if (!defined($arg{sequence}));
    srand();
    my $seqstring = '';
    my $seq_type = ref($arg{sequence});
    if ($seq_type eq 'ARRAY') {
	foreach my $char (@{$arg{sequence}}) {
	    $char = 'U' if ($char eq 'T');
	    $char = 'u' if ($char eq 't');
	    $seqstring .= $char;
	}
    } else {
	$seqstring = $arg{sequence};
	$seqstring =~ tr/Tt/Uu/;
	my @seq = split(//, $seqstring);
	$arg{sequence} = \@seq;
    }
    my $me = bless {
	nt => {0 => 'A', 1 => 'U', 2 => 'C', 3 => 'G',},
	sequence => $arg{sequence},
	seqstring => $seqstring,
	revcomp => $seqstring,
	reverse => $seqstring,
	complement => $seqstring,
	length => scalar(@{$arg{sequence}}),
	aaseq => [],
	aaminusone => [],
	aaminustwo => [],
	dint12seq => [],
	dint23seq => [],
	dint31seq => [],
	dint13seq => [],
	codonseq => [],
	readop => defined($arg{readop}) ? $arg{readop} : 'straight',
	ntfreq => {},
	dintfreq => {},
	codonfreq => {},
    };
    $me->{revcomp} = scalar reverse $me->{revcomp};
    $me->{reverse} = scalar reverse $me->{reverse};
    $me->{revcomp} =~ tr/ACGTUacgtu/UGCAAugcaa/;
    $me->{complement} =~ tr/ACGTUacgtu/UGCAAugcaa/;
    my @ntseq = @{$me->{sequence}};
    $me->{gc_content} = $me->Get_GC(\@ntseq);
    $me->{pyrimidine_content} = $me->Get_CT(\@ntseq);
    my (@aaseq, @aaminusone, @aaminustwo, @codonseq, @dint12seq, @dint23seq, @dint31seq, @dint13seq);
    if ($me->{readop} eq 'orf') {
	my @init = ('A', 'T', 'G');
	next until (splice(@ntseq, 0, 3) eq @init);
    }
    my $finished = 0;
    my $count = 0;
    my ($one, $two, $three, $last_one, $last_two, $last_three);
    while (scalar(@ntseq) > 2 and $finished == 0) {
	$count = $count + 3;
	$me->{dintfreq}{12}{join('', $ntseq[0], $ntseq[1])}++;
	$me->{dintfreq}{23}{join('', $ntseq[1], $ntseq[2])}++;
	$me->{dintfreq}{13}{join('', $ntseq[0], $ntseq[2])}++;
	if ($count == $me->{length}) {
	    $me->{dintfreq}{31}{join('', $ntseq[2], '*')}++;
	} else {
	    $me->{dintfreq}{31}{join('', $ntseq[2], $ntseq[3])}++;
	}
	$last_one = $one;
	$last_two = $two;
	$last_three = $three;
	($one, $two, $three) = splice(@ntseq, 0, 3);
	my $codon_string = join('', $one, $two, $three);
	$codon_string = uc($codon_string);
	$codon_string =~ tr/T/U/;
	my ($minus_one_codon_string, $minus_two_codon_string);
	if (defined($last_three)) {
	    $minus_one_codon_string = join('', $last_three, $one, $two);
	    push(@aaminusone, $references->{codons}{$minus_one_codon_string});
	}
	if (defined($last_two)) {
	    $minus_two_codon_string = join('', $last_two, $last_three, $one);
	    push(@aaminustwo, $references->{codons}{$minus_two_codon_string});
	}

	if (!defined($references->{codons}{$codon_string})) {
	    print "$codon_string is not defined!\n";
	}
	push(@aaseq, $references->{codons}{$codon_string});
	push(@codonseq, $codon_string);
	push(@dint12seq, join('', $one, $two));
	push(@dint23seq, join('', $two, $three));

	#	push(@dint13seq, join('', $one, $ntseq[0]));
	($count == $me->{length}) ? push(@dint31seq, join('', $three, '*')) :
	    push(@dint31seq, join('', $three, $ntseq[0]));
	$me->{codonfreq}{$codon_string}++;
	$me->{ntfreq}{$one}++;
	$me->{ntfreq}{$two}++;
	$me->{ntfreq}{$three}++;

	$finished++ if ($me->{readop} eq 'orf' and ($references->{codons}{$codon_string} eq '*'));
    }    ## End while the length of sequence > 2
    $me->{aaseq} = \@aaseq;
    $me->{aaminusone} = \@aaminusone;
    $me->{aaminustwo} = \@aaminustwo;
    $me->{codonseq} = \@codonseq;
    $me->{dint12seq} = \@dint12seq;
    $me->{dint23seq} = \@dint23seq;
    $me->{dint31seq} = \@dint31seq;
    $me->{dint13seq} = \@dint13seq;

    while (scalar(@ntseq) > 0) {
	my $nt = shift @ntseq;
	$me->{ntfreq}{$nt}++;
	$count++;
    }
    return ($me);
}

=head1 ATTRIBUTION

I think Jonathan Jacobs wrote some of these, but the history of what has been
lost to time I am afraid.

=head2 C<Random>

Randomize a sequence.

=cut
sub Random {
    my $me = shift;
    my @seq = @{$me->{sequence}};
    my @return = ();
    my $count = 0;
    my $codon_count = 0;
    while (scalar(@seq) > 0) {
	shift @seq;
	$return[$count] = $me->{nt}->{int(rand(4))};

	## Avoid premature termination codons using the codon table of interest.
	if ($codon_count == 2) {
	    my $test = $return[$count - 2] . $return[$count - 1] . $return[$count];
	    my $fun = $me->{codons}->{$test};

	    while ($me->{codons}->{$test} eq '*') {
		$return[$count - 2] = $me->{nt}->{int(rand(4))};
		$return[$count - 1] = $me->{nt}->{int(rand(4))};
		$return[$count] = $me->{nt}->{int(rand(4))};
		$test = join('', $return[$count - 2], $return[$count - 1], $return[$count]);
	    }    ## End while we have a ptc
	}    ## End codon check.
	($codon_count == 2) ? $codon_count = 0 : $codon_count++;
	$count++;
    }    ## End top level while
    return (\@return);
}

=head2 C<SameCodons>

Randomize a sequence but maintain codons.

=cut
sub SameCodons {
    my $me = shift;
    my @codons = @{$me->{codonseq}};
    my @return;
    while (scalar(@codons) > 0) {
	my $pull = int(rand(scalar(@codons)));
	my ($first, $second, $third) = split(//, $codons[$pull]);
	push(@return, $first, $second, $third);
	splice(@codons, $pull, 1);
    }
    return (\@return);
}

=head2 C<Get_CT>

Gather the pyrimidine percentage of a sequence.  I don't remember why.

=cut
sub Get_CT {
    my $me = shift;
    my $arr = shift;
    my $par = shift;
    my $ct;
    if (!defined($par)) {
	my $len = scalar(@{$arr});
	my $num_y = 0;
	foreach my $char (@{$arr}) {
	    $num_y++ if ($char eq 'c' or $char eq 't' or $char eq 't' or $char eq 'C' or $char eq 'U' or $char eq 'u');
	}
        $len = 1 if (!defined($len) or ($len == 0));
	$ct = $num_y * 100.0 / $len;
    } else {
	my $num_y = 0;
	my $total = 0;
	for my $c (0 .. $#$par) {
	    next if ($par->[$c] eq '.');
	    $num_y++ if ($arr->[$c] eq 'c' or $arr->[$c] eq 't' or $arr->[$c] eq 'C' or $arr->[$c] eq 'T' or $arr->[$c] eq 'U' or $arr->[$c] eq 'u');
	    $total++;
	}
	$ct = (($total == 0) ? 0 : $num_y * 100.0 / $total);
    }
    $ct = sprintf('%.1f', $ct);
    return($ct);
}

=head2 C<Get_GC>

This gets the more common GC content.

=cut
sub Get_GC {
    my $me = shift;
    my $arr = shift;
    my $par = shift;
    my $gc;
    if (!defined($par)) {
	my $len = scalar(@{$arr});
	my $num_strong = 0;
	foreach my $char (@{$arr}) {
	    $num_strong++ if ($char eq 'c' or $char eq 'g' or $char eq 'G' or $char eq 'C');
	}
        $len = 1 if (!defined($len) or ($len == 0));
	$gc = $num_strong * 100.0 / $len;
    } else {
	my $num_strong = 0;
	my $total = 0;
	for my $c (0 .. $#$par) {
	    next if ($par->[$c] eq '.');
	    $num_strong++ if ($arr->[$c] eq 'c' or $arr->[$c] eq 'g' or $arr->[$c] eq 'C' or $arr->[$c] eq 'G');
	    $total++;
	}
	$gc = (($total == 0) ? 0 : $num_strong * 100.0 / $total);
    }
    $gc = sprintf('%.1f', $gc);
    return($gc);
}

=head2 C<Same31Random>

Maintain the percentages of first, third codon positions when randomizing a sequence.

=cut
sub Same31Random {
    my $me = shift;
    my @dint31seq = @{$me->{dint31seq}};
    my @return;

    my ($third, $first) = undef;
    while (scalar(@dint31seq) > 0) {

	#  The first nucleotide in the codon, checking if we are starting or ending.
	if (defined($first)) {    ## Not at the beginning of the sequence?
	    if ($first eq '*') {    ## Maybe at the end?
		return (\@return);    ## If so, drop out
	    } else {                  ## If not at the end push the first
		push(@return, $first);
	    }
	}         ## End the if not at the beginning of the sequence
	else {    ## at the beginning add a random nucleotide at position 1.
	    push(@return, $me->{nt}->{ int( rand(4))});
	}

	#	!defined($first)    ?   ## 1st nucleotide in the sequence?
	#	  push(@return, $me->{nt}->{int(rand(4))})    :   ## yes
	#		($first eq '*')    ?  ## End of sequence?
	#		  return(\@return)    :  ## yes
	#			push(@return, $first);  ## no

	#  The second nucleotide in the codon
	push(@return, $me->{nt}->{int(rand(4))});

	#  The third nucleotide in the codon
	####  Why do I do the shift here?  Think about the recreation of the nucleotide sequence
	####  for a moment in terms of the first codon and you will see :)
	($third, $first) = split(//, shift(@dint31seq));
	push(@return, $third);
    }    ## End while.
    return (\@return);
}

=head2 C<Sameaa>

Maintain amino acid frequencies but choose a random codon from those available for each amino acid.

=cut
sub Sameaa {
    my $me = shift;
    my $aaref = shift;
    my @aaseq = @{$aaref};
    my @return = ();
    while (scalar(@aaseq) > 0) {
	my $aa = shift @aaseq;
	my ($first, $second, $third) = split(//, $me->{codons}->{$aa}->[0]->[int(rand(scalar(@{$me->{codons}->{$aa}->[0]})))]);
	push(@return, $first, $second, $third);
    }
    return (\@return);
}

=head2 C<Same31Composite>

Step 1:  Acquire a sequence to randomize
Step 2:  pick an amino acid to consider out of pool of amino acids not yet considered.
Step 3:  Choose two random codons in the sequence which code for this amino acid.
  "Consider the consequences of swapping the 1st and jth codons."
 Make the swap if the 1,3 dinucleotides do not change after the swap.
 If they do change, make a list of all reciprocal swaps remaining, pick one at random, make both swaps.
 Mark both pairs as 'swapped'
 Recurse for all other codons for all amino acids.

=cut
sub Same31Composite {
    my $me = shift;
    my @dint31seq = @{$me->{dint31seq}};
    my @aminos = @{$me->{aaseq}};
    my @codons;

    my $first = undef;
    my ($third) = split(//, $dint31seq[0]);
    while (scalar(@dint31seq) > 0) {
	my $aa = shift @aminos;
	#  The first nucleotide in the codon, checking if we are starting or ending.
	!defined($first)
	    ?    ## 1st nucleotide in the sequence?
	    push(@codons, $me->Find_Alternates($aa, "..$third"))
	    :    ## yes
	    ($first eq '*')
	    ?    ## End of sequence?
	    return (\@codons)
	    :    ## yes
	    push(@codons, $me->Find_Alternates($aa, "$first.$third"));

	#			push(@codons, $first);  ## no
	#  The second nucleotide in the codon
	#	push(@codons, $me->Find_Alternates($aa, "$first.$third"));
	#  The third nucleotide in the codon
	####  Why do I do the shift here?  Think about the recreation of the nucleotide sequence
	####  for a moment in terms of the first codon and you will see :)
	($third, $first) = split(//, shift(@dint31seq));
    }    ## End while.
    return (\@codons);
}

=head2 C<Find_Alternates>

Find alternate codons for a given codon.

=cut
sub Find_Alternates {
    my $me = shift;
    my $amino = shift;
    my $spec = shift;
    my @return = ();
    my $starters = $me->{codons}->{$amino}->[0];    ### An array reference like ['GCA','GCC','GCG','GCU']
    return ($starters) if ($spec eq '...');
    my ($dumbone, $dumbtwo, $dumbthree) = split(//, $spec);
    if ($dumbtwo eq '.' and $dumbone eq '.') {
        foreach my $codon (@{$starters}) {
            my ($one, $two, $three) = split(//, $codon);
            push(@return, $codon) if ($dumbthree eq $three);
        }
    } elsif ($dumbtwo eq '.') {
        foreach my $codon (@{$starters}) {
            my ($one, $two, $three) = split(//, $codon);
            push(@return, $codon) if ($dumbthree eq $three and $dumbone eq $one);
        }
    } else {
        print "What!?\n";
    }
    return (\@return);
}

=head2 C<Translate>

Translate some sequence.  There are so many versions of this function, I bet
this one is amont the worst.

=cut
sub Translate {
    my $me = shift;
    my $sequence = $me->{seqstring};
    $sequence =~ tr/atgcT/AUGCU/;
    if ((!$me->{sequence}) and (!defined($sequence))) {
        die("Nothing to work with.");
    } elsif (defined($sequence)) {
        my @seqtmp = split(//, $sequence);
        my $t = new SeqMisc(sequence => \@seqtmp);
        my $aaseq = join("", @{$t->{aaseq}});
        undef $t;
        return ($aaseq);
    } else {
        my $aaseq = join("", @{$me->{aaseq}});
        return ($aaseq);
    }
}

=head2 C<Coin_Random>

Expect an array reference, return same reference, but array changed.
Just do a pseudo random number randomization of each incoming nucleotide
Code shamelessly taken from Jonathan Jacobs' <jlj email> work from 2003.

=cut
sub Coin_Random {
  my $me = shift;
  my $sequence = shift;
  my @nucleotides = ('a','u','g','c');
  for my $c ($#$sequence) {
    $sequence->[$c] = $nucleotides[int(rand(4))];
  }
  return ($sequence);
}

=head2 C<Nucleotide_Montecarlo>

Expect an array reference, Perform a shuffle by substitution
Code shamelessly taken from Jonathan Jacobs' <jlj email> work from 2003.

=cut
sub Nucleotide_Montecarlo {
  my $start_sequence = shift;
  my @st = @{$start_sequence};
  my $new_sequence = [];
  my $new_seq_count = 0;
  while (@st) {
    my $position = rand(@st);
    $new_sequence->[$new_seq_count] = $position;
    $new_seq_count++;
    splice(@st, $position, 1);
  }
  return ($new_sequence);
}

=head2 C<Dinucleotide>

Expect an array reference
Perform a shuffle keeping the same dinucleotide frequencies as original sequence.
Code shamelessly taken from Jonathan Jacobs' <jlj email> work from 2003.
I don't know why, but I love the way Jonathan wrote this.

=cut
sub Dinucletode {
    my $start_sequence = shift;
    my @st = @{$start_sequence};
    my @new_sequence = ();

    #  my @startA = my @startT = my @startG = my @startC = ();
    #  ## Create 4 arrays containing
    ## Create an array containing all the dinucleotides
    my @doublets = ();
    push(@doublets, $start_sequence =~ /(?=(\w\w))/g);
    while (scalar(@{$start_sequence}) > scalar(@new_sequence)) {
	my $chosen_doublet = int(rand(@doublets));
	my ($base1, $base2) = split(//, $doublets[$chosen_doublet]);
	push(@new_sequence, $base1, $base2);
    }
    return (\@new_sequence);
}

=head2 C<Codon_montecarlo>

Expect an array reference
Perform a shuffle keeping the same triplet frequencies as original sequence.
Code shamelessly taken from Jonathan Jacobs' <jlj email> work from 2003.

=cut
sub Codon_montecarlo {
    my $start_sequence = shift;
    my @new_sequence = ();
    my $new_seq_count = 0;
    my @codons = $$start_sequence =~ /(\w\w\w)/g;
    while (@codons) {
	my $position = rand(@codons);
	my @nts = split(//, $codons[$position]);
	push(@new_sequence, @nts);
	splice(@codons, $position, 1);
    }
    return (\@new_sequence);
}

=head2 C<Related_Codon>

Expect an array reference
Perform a shuffle keeping the same triplet frequencies as original sequence.
Code shamelessly taken from Jonathan Jacobs' <jlj email> work from 2003.
This will clearly not work anymore.

=cut
sub Related_Codon {
    my $start_sequence = shift;
    my @codons = $start_sequence =~ /(\w\w\w)/g;
    for my $c (0 .. $#codons) {
	my @potential = split(/\s/, $Randomize::amino_acids{$codons[$c]});
	$codons[$c] = $potential[ rand(@potential)];
    }
    return (\@codons);
}

=head2 C<Array_Shuffle>

This shuffles a referenced array like a deck of cards.
From the Perl cookbook. Uses the Fischer-Yates Shuffle.

=cut
sub ArrayShuffle {
    my $seqArrayREF = shift;
    my @arrayREF = @$seqArrayREF;
    for (my $i = @arrayREF ; --$i ;) {
	my $j = int rand($i + 1);
	next if $i == $j;
	@arrayREF[$i, $j] = @arrayREF[$j, $i];
    }
    return (\@arrayREF);
}

=head2 C<Squid>

Invoke squid and use its better methods for shuffling.

=cut
sub Squid {
    my $inarray = shift;
    my $shuffle = shift;
    my $inseq = '';
    my $shuffle_exe;
    foreach my $char (@{$inarray}) { $inseq = join('', $inseq, $char); }
    if (defined($shuffle)) { $shuffle_exe = $shuffle; }
    else { $shuffle_exe = 'shuffle'; }
    my $out_text;
    {    ## Begin a File::Temp Block
        my $fh = new File::Temp(DIR => qq"$ENV{HOME}/tmp", UNLINK => 0,);
        ## OPEN $fh in Squid
        print $fh ">tmpsquid
$inseq
";
        my $infile = $fh->filename;
        my $command = "$shuffle_exe $infile";
        open(CMD, "$command |");
        ## OPEN CMD in Squid
        while (my $line = <CMD>) {
            chomp $line;
            next if ($line =~ /^\>/);
            $out_text = join('', $out_text, $line);
        }    ## End while
        close(CMD);
        unlink($infile);
        ## CLOSE CMD in Squid
    }    ## End a File::Temp Block -- the tempfile should now no longer exist.
    my @out_array = split(//, $out_text);
    return (\@out_array);
}

=head2 C<Squid_Dinuc>

Use Squid to dinucleotide shuffle some sequence.

=cut
sub Squid_Dinuc {
    my $inarray = shift;
    my $shuffle = shift;
    my $inseq = '';
    my $shuffle_exe;
    foreach my $char (@{$inarray}) { $inseq = join('', $inseq, $char); }
    if (defined($shuffle)) { $shuffle_exe = $shuffle; }
    else { $shuffle_exe = qq($ENV{HOME}/bin/shuffle); }
    my $out_text = '';
    {    ## Begin a File::Temp Block
        my $fh = new File::Temp(DIR => qq"$ENV{HOME}/tmp", UNLINK => 0,);
        print $fh ">tmp
$inseq
";
        my $infile = $fh->filename;
        my $command = qq($shuffle_exe -d $infile);
        open(CMD, "$command |");
        ## OPEN CMD in Squid_Dinuc
        while (my $line = <CMD>) {
            chomp $line;
            next if ($line =~ /^\>/);
            $out_text = join('', $out_text, $line);
        }
        close(CMD);
        unlink($infile);
    }    ## End a ifile::temp block
    my @out_array = split(//, $out_text);
    return (\@out_array);
}

=head2 C<Codon_Distribution>

Gather the distribution of codons from a sequence.

=cut
sub Codon_Distribution {
    my $me = shift;
    my %args = @_;
    my $residues = ['*'] unless ($args{residues});
#    my $sequence = $me->{aaminusone} unless ($args{sequence});
    my $sequence = ($args{sequence} ? $args{sequence} : $me->{aaminusone});
    my $dist_array = ($args{dist_array} ? $args{dist_array} : []);
    my @seq = @{$sequence};
    my $dist = 0;
    LOOP: while (@seq) {
	my $current = shift(@seq);
	foreach my $res (@{$residues}) {
	    if ($current eq $res) {
		$dist_array = $me->Add_One($dist, $dist_array);
		$dist = 0;
		next LOOP;
	    }
	}
	$dist++;
    } ## End of the while loop
    return($dist_array);
}

=head2 C<Add_One>

Count!

=cut
sub Add_One {
    my $me = shift;
    my $pos = shift;
    my $arr = shift;
    my @array = @{$arr};
    my $len = $#array;
    if ($len >= $pos) {
        $arr->[$pos]++;
        return($arr);
    } else {
        my $c = 0;
        while ($c < $pos) {
            $array[$c] = 0 unless ($array[$c]);
            $c++;
        }
        return(\@array);
    }
}

=head1 AUTHOR - atb

Email  <abelew@gmail.com>

=head1 SEE ALSO

    L<Bio::Tools::Run::StandAloneBlast>

=cut

1;
