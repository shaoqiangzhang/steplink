#!/usr/bin/perl 
use warnings;
use strict;
#-----------------------------------------------------------
#Produce hash_of_hash of k-mers and link kmers into contigs
#By Shaoqiang Zhang(zhangshaoqiang@mail.tjnu.edu.cn)
#Apr16, 2016 
#-----------------------------------------------------------
while(@ARGV<2){
	print "\nUSAGE:\n\tprogram.pl <fastq_file> <number k of k-mer>\n\n";
	die;
}

my $fastq_file=shift @ARGV;
my $lengthkmer=shift @ARGV;

open(FASTQ,$fastq_file) or die "ERROR: not a fastq file!!\n";


my %hash_of_hash_kmer=(); # kmer and its right kmer
my %first_kmer_hash=(); ##the first k-mer of each read
my %right_kmer_hash=(); ##the k-mers in each read except the first one 
my %left_kmer_hash=(); 
my $nextline;
my @linearray;
my $kmer; my $rkmer; 
my $containN=0;
my @kmerarray;
my $currentkmer; my $rightkmer;

##==========the following is to construct the hash of kmer's hashes ===================##
while(<FASTQ>){
	if($_=~/^@.+/){ #@ is the leftmost label of a fastq file
		$nextline=<FASTQ>;
		chomp $nextline;
		if($nextline=~/^[ACGTN]+/){
			$currentkmer=substr($nextline,0,$lengthkmer);
			unless($currentkmer=~/N+/){
				$first_kmer_hash{$currentkmer}=1;
			}
			for(my $i=1;$i<$lengthkmer;$i++){
				$currentkmer=substr($nextline,$i,$lengthkmer);
				$right_kmer_hash{$currentkmer}=1;
			}
			for (my $i=0; $i<=(length($nextline)-2*$lengthkmer); $i++){
				$currentkmer=substr($nextline,$i,$lengthkmer);
				$rightkmer=substr($nextline,$i+$lengthkmer,$lengthkmer);
				$right_kmer_hash{$rightkmer}=1;
				unless(($currentkmer=~/N+/) or ($rightkmer=~/N+/)){
					unless(exists($hash_of_hash_kmer{$currentkmer}{$rightkmer})){
						$hash_of_hash_kmer{$currentkmer}{$rightkmer}=1;
					}else{
						$hash_of_hash_kmer{$currentkmer}{$rightkmer}++;
					}
				}
			}
		}
	}
}


##===========below is label the first_left_kmer of each transcript===================##
my @keys_num_array=();
my %first_left_kmer_hash=();##the first left kmer of all transcripts

foreach $kmer (sort keys (%first_kmer_hash)){
	unless(exists($right_kmer_hash{$kmer})){####
		$first_left_kmer_hash{$kmer}=1;
	}
}
undef %first_kmer_hash;
undef %right_kmer_hash;

print scalar(keys (%first_left_kmer_hash)),"\n";

my %temp_hash=();
my @contigs_array=();

foreach $kmer (sort keys(%first_left_kmer_hash)){##link the basic contigs
	my $contig_seq=$kmer;
	my $next_kmer=$kmer;
	while(1){
		if(exists($hash_of_hash_kmer{$next_kmer})){
			%temp_hash=%{$hash_of_hash_kmer{$next_kmer}};
			my @num_right_kmer=sort { $temp_hash{$b} <=> $temp_hash{$a} } keys(%temp_hash); 
			$next_kmer=shift @num_right_kmer;
			$contig_seq=$contig_seq.$next_kmer;
		}else{
			my $contain_label=0;
			for(my $i=1;$i<$lengthkmer;$i++){
				my $inner_kmer=substr($contig_seq, length($contig_seq)-$lengthkmer-$i, $lengthkmer);
				if(exists($hash_of_hash_kmer{$inner_kmer})){
					$contain_label=1;
					my $new_contig_seq=substr($contig_seq,0,length($contig_seq)-$i);
					%temp_hash=%{$hash_of_hash_kmer{$inner_kmer}};
					my @num_right_kmer=sort { $temp_hash{$b} <=> $temp_hash{$a} } keys(%temp_hash);
					$next_kmer=shift @num_right_kmer;
					$contig_seq=$new_contig_seq.$next_kmer;
					last;
				}
			}
			if($contain_label == 0){
				last;
			}
		}
	}	
	push (@contigs_array,$contig_seq);
}

foreach my $seq (@contigs_array){##delete the used hash_of_hash_kmer's keys and values
	@linearray=split(//,$seq);
	for (my $i=0; $i<=(scalar(@linearray)-2*$lengthkmer); $i++){
		my $currentkmer="";
		my $rightkmer="";
		for(my $j=0;$j<$lengthkmer;$j++){
			$currentkmer=$currentkmer.$linearray[$i+$j];
		}
		for(my $k=$lengthkmer;$k<2*$lengthkmer;$k++){
			$rightkmer=$rightkmer.$linearray[$i+$k];
		}
		delete $hash_of_hash_kmer{$currentkmer}{$rightkmer};
		unless (%{$hash_of_hash_kmer{$currentkmer}}){
			delete $hash_of_hash_kmer{$currentkmer};
		}
	}
} 

###=======================================================================
my @affix_array=();
my @temp_array=@contigs_array;
while(%hash_of_hash_kmer){ ###===find all affix sequences
	my @one_affix_array=();
	foreach my $contig_seq (@temp_array){
		for(my $i=0; $i<=(length($contig_seq)-2*$lengthkmer); $i++){
			$currentkmer=substr($contig_seq,$i,$lengthkmer);
			if(exists($hash_of_hash_kmer{$currentkmer})){
				my $contig_affix=$currentkmer;
				while(1){
					if(exists($hash_of_hash_kmer{$currentkmer})){
						%temp_hash=%{$hash_of_hash_kmer{$currentkmer}};
						my @num_right_kmer=sort { $temp_hash{$b} <=> $temp_hash{$a} } keys(%temp_hash); 
						$currentkmer=shift @num_right_kmer;
						$contig_affix=$contig_affix.$currentkmer;
					}else{
						my $contain_label=0;
						for(my $i=1;$i<$lengthkmer;$i++){
							my $inner_kmer=substr($contig_affix, length($contig_affix)-$lengthkmer-$i, $lengthkmer);
							if(exists($hash_of_hash_kmer{$inner_kmer})){
								$contain_label=1;
								my $new_contig_affix=substr($contig_affix,0,length($contig_affix)-$i);
								%temp_hash=%{$hash_of_hash_kmer{$inner_kmer}};
								my @num_right_kmer=sort { $temp_hash{$b} <=> $temp_hash{$a} } keys(%temp_hash);
								$currentkmer=shift @num_right_kmer;
								$contig_affix=$new_contig_affix.$currentkmer;
								last;
							}
						}
						if($contain_label == 0){
							last;
						}
					}	
				}
				push(@one_affix_array,$contig_affix);
				@linearray=split(//,$contig_affix);
				for (my $i=0; $i<=(scalar(@linearray)-2*$lengthkmer); $i++){##delete used kmers' hash of hashes
					my $currentkmer="";
					my $rightkmer="";
					for(my $j=0;$j<$lengthkmer;$j++){
						$currentkmer=$currentkmer.$linearray[$i+$j];
					}
					for(my $k=$lengthkmer;$k<2*$lengthkmer;$k++){
						$rightkmer=$rightkmer.$linearray[$i+$k];
					}
					delete $hash_of_hash_kmer{$currentkmer}{$rightkmer};
					unless (%{$hash_of_hash_kmer{$currentkmer}}){
						delete $hash_of_hash_kmer{$currentkmer};
					}
				}
				last;
			}					
		}
	}
	push(@temp_array,@one_affix_array);
	push(@affix_array,@one_affix_array);
}
##===========================================================================#

my @new_contigs_array=@contigs_array;
while(@affix_array){
	my @temp_contigs_array=();
	my @temp_affix_array=@affix_array;
	foreach my $mycontig(@new_contigs_array){
		my $temp_contig_seq=$mycontig;
		my $new_contig_seq="";
		foreach my $myaffix(@temp_affix_array){
			my $prefix=substr($myaffix,0, $lengthkmer);
			my $no_suffix=substr($myaffix,0,(length($myaffix)-$lengthkmer));
			my $suffix=substr($myaffix,-$lengthkmer);
			if($temp_contig_seq=~/^(.*)$prefix.*($suffix.*)$/){
				$new_contig_seq=$new_contig_seq.$1.$no_suffix;
				$temp_contig_seq=$2;
				@affix_array=grep{$_ ne $myaffix} @affix_array;
			}elsif($temp_contig_seq=~/^(.*)$prefix/){
				$new_contig_seq=$new_contig_seq.$1.$myaffix;
				$temp_contig_seq="";
				@affix_array=grep{$_ ne $myaffix} @affix_array;
			}
		}
		if(length($new_contig_seq)>0){
			$new_contig_seq=$new_contig_seq.$temp_contig_seq;
			push (@temp_contigs_array,$new_contig_seq);
		}	
	}
	push(@contigs_array,@temp_contigs_array); 
	@new_contigs_array=@temp_contigs_array;
}


##==========================================================================#

open(CONTIG, ">$fastq_file.$lengthkmer.contigs");
foreach my $mycontig (sort{length($b)<=>length($a)}(@contigs_array)){
	print CONTIG $mycontig,"\n";
}

open(AFF, ">$fastq_file.$lengthkmer.affix");
foreach my $myaffix (sort{length($b)<=>length($a)}(@affix_array)){
	print AFF $myaffix,"\n";
}


