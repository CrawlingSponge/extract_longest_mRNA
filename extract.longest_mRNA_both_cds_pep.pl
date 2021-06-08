#!/usr/bin/perl
use strict;
# scripted by YuanJ to extract longest mRNA, and the trans.!

my $dir = "./"; # you need change this path into your species directories
opendir(DIR, "$dir") or die "$!";
my @species = grep(/^\D\D\D$/, readdir(DIR));# desigend the "species abbre" directory with three characters, for example, "mum" represent Mus musculus, common name as mouse.
close DIR;

foreach my $each (@species){
	if(-e "$each.cds.fasta"){next;}
	print "\nRound $each\n";
	my $data_path = "$dir/$each/NCBI";# change this part into your genomic path, for example, the genomic datas of mouse within the directory "./mum/NCBI/*"
	opendir(SDIR, "$data_path") or die "$!";
	my @files = grep(/\S+/, readdir(SDIR));
	close SDIR;

	my $gff; my $pep;  my $sp = $each; my $type; my @fastas = ();
	sort @files;
	foreach my $file (@files){
		if($file =~ m/\.gff$/){$gff = $file;}#the standard file name of gff, tailed with ".gff"
		if($file =~ m/_cds_from_genomic\.fna$/){$pep = $file; push @fastas, $pep;}#push cds, the standard file name of cds, tailed with "_cds_from_genomic.fna"
		if($file =~ m/\.faa$/){$pep = $file; push @fastas, $pep;}#push pep, the standard file name of pep, tailed with ".faa"
	}
	
	foreach my $fa (@fastas){
		print "$gff\n";
		if($fa =~ m/_cds_from_genomic\.fna$/){$pep = $fa; $type = "cds"; print "$pep\n";}#this part just for cds
		if($fa =~ m/\.faa$/){$pep = $fa; $type = "pep"; print "$pep\n";}#this part just for cds

		open(GFF, "$data_path/$gff") and print "step 1: extracting longest mRNA from ori_gff\n" or die "cannot open this file $!";
		#open(MR, ">mrna.all_long.gff") or die "cannot open this file $!";# you can also print the longest mRNA infor
		my $row = 0;
		my (%hash, $max_devalue, $last_gene, @m_line, @m_id, $mrnaid, %hash_m,);
		while(my $line = <GFF>){
			chomp($line);
			my @arr = split(/\s+/, $line); 
			if($arr[2] !~ m/^mRNA$/){next;}
			$row++;
			my @a = split(/;/, $arr[8]);
			
			#my @b = split(/,/, $a[2]); #
			#my $geneid = $b[0]; $geneid =~ s/^Dbxref=//; #
			#my $geneid = $b[0]; $geneid =~ s/^\D+\=//; # this part is geneid, with repreat gene_id, and will filter more genes.
			
			my $geneid = $a[1]; $geneid =~ s/^\D+\=//; # for general gff, actually this geneid means genename in this situation
			
			my $devalue = abs($arr[3] - $arr[4]);
			$hash{$devalue} = $line;
			
			if($row == 1){
				$max_devalue = $devalue;
			}else{
				if($last_gene ne $geneid){
					#print MR "$hash{$max_devalue}\n";
					## this paras mainly for store the mRNA ID for following find the CDS;
					@m_line = split(/\s+/, $hash{$max_devalue});
					@m_id = split(/;/, $m_line[8]);
					$mrnaid = $m_id[0]; $mrnaid =~ s/^ID=rna-//;
					$hash_m{$mrnaid} = 1;
					
					$max_devalue = $devalue;
					%hash = ();## clear all hash value and keys within last genes;
					$hash{$devalue} = $line; ## dont forget this line also need to initiate the hash;
				}else{
					if($devalue >= $max_devalue){
						$max_devalue = $devalue;
					}else{
						$max_devalue = $max_devalue;
					}
				}
			}
			$last_gene = $geneid;
		}
		#print MR "$hash{$max_devalue}\n";
		@m_line = split(/\s+/, $hash{$max_devalue});
		@m_id = split(/;/, $m_line[8]);
		$mrnaid = $m_id[0]; $mrnaid =~ s/^ID=rna-//;
		#print "last_mRNA of GFF: $mrnaid\n";
		$hash_m{$mrnaid} = 1;
		close GFF;
		close MR;


		# step2 filter the cds within longest mRNA and also filter the mRNA without translating to CDS;
		open(GFF, "$data_path/$gff") and print "step 2: filter the trans and mRNA\n" or die "cannot open this file $!";
		my (%hash_cds, %hash_rna);
		while(my $line = <GFF>){
			chomp($line);
			my @arr = split(/\s+/, $line);
			next if $arr[2] ne "CDS";
			my @a = split(/;/, $arr[8]);
			my $rna = $a[1]; $rna =~ s/^Parent=rna-//;
			my $cds = $a[0]; $cds =~ s/^ID=cds-//; 
			if(exists $hash_m{$rna}){
				#print "$rna\t$cds\n";
				$hash_cds{$cds} = 1;
				$hash_rna{$rna} = 1;
			}
		}
		close GFF;


		# step 3 filter cds and rename id, this part did not change the name , but we used this varis to deal cds fasta
		open(CDS, "$data_path/$pep") and print "step 3: filter $type\n" or die "$!";
		my ($pep_id, %hash_pep);
		while(my $line = <CDS>){
			chomp($line);
			if($line =~ m/^>(\S+)\s/){
				$pep_id = $1;
				if($line =~ m/^>\S+_cds_(\S+)_\d+\s+/){
					$pep_id = $1;
				}else{
					next;#non pep like ">lcl|NW_003615950.1_cds_33519" without pep
				}
			}else{
				$hash_pep{$pep_id} .= $line;
			}
		}
		close CDS;

		open(OC, ">$sp.$type.fasta") or die "$!";
		my @ord_pep = sort(keys %hash_pep);
		foreach my $each (@ord_pep){
			if(exists $hash_cds{$each}){
				print OC ">$sp"."_"."$each\n$hash_pep{$each}\n";
			}
		}
		close OC;
	}
}