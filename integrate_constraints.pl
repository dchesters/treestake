

# 
# 
# integrate_constraints.pl, started 2018-11-02
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
##########################################################################################################




my $arg_string  = join ' ', @ARGV;




#####################################
read_command_arguments($arg_string);#
#####################################





###############################################################################################################################################

# from parse_ncbi_tax_db ...

$ignore_subspecies			= 1;	# 0=read taxonomies beyond species level (e.g. include subspecies, varietas, and unnamed bacteria ranks)
						# if 1, these extra ranks are ignored
						# this applies to full species name, assigned to given ncbi tax number
# globals

%ncbi_nodes;	# ncbi taxonomy structure is recorded in this object

%assign_these_new_terminals_to_taxon;


###############
store_nodes();#
###############



###################
parse_namesfile();#
###################




$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;$starting_name_ff =~ s/\s+/_/g;
print "name for the starting node ($starting_node)\nwhich has been specified ($starting_name)\n";
print "traversing NCBI taxonomy tree from this node. recording lineage for each taxon.\n\n";
$starting_name =~ s/^(\w).+$/$1/;

#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################


# hack, some missing species
# in other scripts, the presence of these when they actually are in the ncbi tax DB (e.g. newer db is used), causes error.

#$complete_lineage_for_this_species{"Lemyra_melli"} 		= $complete_lineage_for_this_species{"Lemyra_jankowskii"};
#$complete_lineage_for_this_species{"Promalactis_suzukiella"} 	= $complete_lineage_for_this_species{"Promalactis_jezonica"};
#$complete_lineage_for_this_species{"Hylaeus_dilatatus"} 		= $complete_lineage_for_this_species{"Hylaeus_difficilis"};
#$complete_lineage_for_this_species{"Lethe_dura"} 			= $complete_lineage_for_this_species{"Lethe_diana"};
#$complete_lineage_for_this_species{"Eucryptorrhynchus_brandti"} 	= $complete_lineage_for_this_species{"Eucryptorhynchus_brandti"};




###############################################################################################################################################



for my $i(0 .. $#constraint_files)
	{
	my $current_file = $constraint_files[$i];
	print "\n\nconstraint file $i named $current_file\n";	
	open(IN, $current_file) || die "\ncant open file $current_file\n";
	my $line_count = 0;my $print_limit = 0;
	my $entries_user_taxon = 0; my $entries_not_user_taxon = 0;
	my $lineage_not_retreived =0;
	my $informative_entries = 0;my $non_informative_entries = 0;

	while (my $line = <IN>)
		{
		$line_count++;
		my $skip_line = 0;
		if($line_count == 1)
			{
			print "header:$line";
			if($line =~ /^\s+\d+\s+\d+/)
				{$skip_line = 1; print "phylip header, ignoring line\n"};
			};


		##############################################################################
		unless($skip_line == 1)
			{
			$line =~ s/\n//;$line =~ s/\r//;

			if($line =~ /^([^\s\t]+)[\s\t]+(\S+)$/ )
				{
				my $id = $1; my $constraints = $2;$print_limit++;
				if($print_limit == 6){print ".....\n"}elsif($print_limit <= 5){print "$id\n"};

				my $complete_lineage = $complete_lineage_for_this_species{$id};
				my $lineage_retreived = 0;
				if($complete_lineage =~ /\:\w/)
					{
					$lineage_retreived = 1;
 					}else{
					my $genusname = $id;$genusname =~ s/^([A-Z][a-z]+)_.+/$1/;
					$complete_lineage = $complete_lineage_for_this_species{$genusname};
					if( $complete_lineage =~ /\:\w/ )
						{$lineage_retreived = 1};
					}

				if($limit_taxon_rank =~ /./){}else{$lineage_retreived = 1};
				if($lineage_retreived == 1)
					{
				#	print "complete_lineage:$complete_lineage\n";
					my $test = 0;
					if($complete_lineage =~ / $limit_taxon_rank:$limit_taxon /)
						{$test = 1};
					if($limit_taxon_rank =~ /./){}else{$test = 1};
					
					if($test == 1)	
						{
						$entries_user_taxon++;#print "\tuser taxon\n";
						if($constraints =~ /[01]/)
							{
							$store_all_IDs{$id} = 1;
							$store_constraints{$i}{$id} = $constraints;
							$informative_entries++;
							}else{
							$non_informative_entries++;
							};
						}else{
						$entries_not_user_taxon++;#print "\tsomething else\n";
						};
					
					}else{
					$lineage_not_retreived++
					};

				$missing_chars{$i} = $constraints;$missing_chars{$i} =~ s/[01\?]/-/g;


				}else{
				print "warning, parse error:$line\n";
				if($line =~ /^>[A-Z][a-z]/){$looks_like_fasta++};
				};
			
			};# unless($skip_line == 1)
		##############################################################################


		};
	close IN;

	print "finished reading constraint file, line count:$line_count\n";
	print "lineage_not_retreived:$lineage_not_retreived\n";
	print "entries_user_taxon:$entries_user_taxon entries_not_user_taxon:$entries_not_user_taxon\n";
	print "informative_entries:$informative_entries non_informative_entries:$non_informative_entries\n";

	if($looks_like_fasta >= 50){die "\n\nerror, looks like Fasta input, expecting something else. quitting.\n\n"};

	print "\n";
	};

############################################################################################


# constraints now stored for user taxon.

print "
constraints stored for user taxon. going through all members and concatenating....
";

open(OUT, ">constraints_concat");
open(OUT2, ">constraints_concat.ID_list");

my @all_species = keys %store_all_IDs;@all_species = sort @all_species;
my $final_species_count =0;
my $all_species_count =0;

foreach my $sp(@all_species)
	{
	my $print_string = "$sp     ";my $partitions_present = 0;
	my $seq_only;
	for my $j(0 .. $#constraint_files)
		{
		my $missing = $missing_chars{$j};
		if($store_constraints{$j}{$sp} =~ /[01]/)
			{
			$seq_only .= $store_constraints{$j}{$sp}; $print_string .= $store_constraints{$j}{$sp};$partitions_present++;
			}else{
			$seq_only .= $missing; $print_string .= $missing;
			};

		};
	if($partitions_present >= 1)
		{
		$final_species_count++;

	#	if($print_string =~ /[01]{100}/)
	#		{}else{
			$big_print_string .= "$print_string\n"; # print OUT "$print_string\n"; 
			print OUT2 "$sp\n";
	#		};

		if($final_species_count == 1)
			{
			$ref_length = length($seq_only)
			}else{
			unless($ref_length == length($seq_only)){print "warning, length mismatch\n"};
			};

		if($print_string =~ /^[\w\s]+ \-+$/){die "\nhuh:$print_string\n"};
		unless($print_string =~ /^\w\S+ +[012\-]+$/){die "\nafter concatenating, unexpected structure:$print_string\n"};

		}else{
		print "warning, no partitions for $sp\n"
		};

	};



close OUT2;

$new_constraints_length = 0;
for my $j(0 .. $#constraint_files)
	{
	my $missing = $missing_chars{$j};$new_constraints_length += length($missing);
	};
print "
done, dimensions of new constraint matrix:
 $final_species_count $new_constraints_length
";
print OUT " $final_species_count $new_constraints_length
$big_print_string";

close OUT;


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub read_command_arguments
{

my $arguments = shift;
$arguments =~ s/\n//;$arguments =~ s/\r//;

print "
\n\n
     ***** integrate_constraints.pl *****
		  
\n";

# -constraint_files taxon_constraints_only fasttree_format_constraint -fasta $sp_matrix



if($arguments =~ /-fasta\s+(\S+)/)
	{
	$fasta_file = $1;
	}else{
	print "\nerror reading command arguments  (-fasta)\n";die""
	}

if($arguments =~ /-constraint_files\s+([^\-]+)/)
	{
	$constr_files = $1;
	@constraint_files = split /\s+/ , $constr_files;
	}else{
	print "\nerror 2. reading command arguments \n";die""
	}

if($arguments =~ /-limit_taxon\s+(\w+)\:(\w+)/)
	{
	$limit_taxon_rank = $1;$limit_taxon = $2;
	}else{
	print "\nuser not chosen to limit taxon\n";
	}

if($arguments =~ /-node\s+(\d+)/)
	{
	$starting_node = $1;
	}else{
	print "user did not given NCBI taxonomic number. using default 33208 (Metazoa)
if you have references outside this default, you need to get the appropriate NCBI taxonomy number from http://www.ncbi.nlm.nih.gov/taxonomy/
and input this here with the option -node taxon_number
";
	$starting_node = 33208;
	}




print "
fasta_file:$fasta_file
user given " , scalar @constraint_files , " constraint files: @constraint_files
limit_taxon_rank:$limit_taxon_rank limit_taxon:$limit_taxon

";



}#sub read_command_arguments



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub store_nodes
{

my %rank_hash;

# nodes.dmp contains each taxon, described in a tree-like structure. 
# so each node has a name (eg polyphaga) and the higher group node to which it belongs (beetles). 

open(NODES , "nodes.dmp") || die "cant open nodes.dmp
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";

print "\nreading files from NCBI taxonomy database .... nodes.dmp .... ";


my $line_counter=0;
while (my $line = <NODES>)
	{
	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]+)\t[\|]\t/)
		{
		my $tax_id = $1;my $parent_tax_id = $2;my $rank = $3;
		$rank_hash{$rank}++;
#		print "tax_id:$tax_id parent_tax_id:$parent_tax_id rank:$rank\n";

			$ncbi_nodes{$tax_id}{rank} = $rank;
		#	$ncbi_nodes{$tax_id}{rank_code} = $current_rankcode;
			$ncbi_nodes{$tax_id}{parent} = $parent_tax_id;
			$ncbi_nodes{$parent_tax_id}{child_nodes} .= "\t" . $tax_id;

		}else{
		print "line_counter:$line_counter line:$line";
		die "UNEXPECTED LINE:$line\nquitting\n";
		}
	$line_counter++;
	}

close(NODES);

#my @ranks = keys %rank_hash;@ranks = sort @ranks;
#print "ranks found in nodes.dmp:\n";
#print LOG "ranks found in nodes.dmp:\n";

#foreach my $rank(@ranks){print "$rank\t" , $rank_hash{$rank} , "\n";print LOG "$rank\t" , $rank_hash{$rank} , "\n"};

my @all_nodes = keys %ncbi_nodes;@all_nodes = sort @all_nodes;

print scalar @all_nodes , " nodes have been read.\n";


}




#####################################################################################################
#
#
#
#####################################################################################################


sub parse_namesfile
{

# here just parse the scientific name of each node. ignore synonyms etc

open(NAMES , "names.dmp") || die "cant open names.dmp
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";


print "\nnames.dmp, parsing 'scientific name', ignoring others ... ";

my $names_line_counter=0;
while (my $line = <NAMES>)
	{
# 24	|	Shewanella putrefaciens	|		|	scientific name	|

	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]*)\t[\|]\tscientific name/)
		{
		my $tax_id = $1;my $name = $2;#my $rank = $3;
		# print "tax_id:$tax_id name:$name\n";

		# if you want to remove non-alphanumerical characters from assigned species names:
		$name =~ s/[\(\)\,\[\]\'\#\&\/\:\.\-]/ /g;
		$name =~ s/\s\s+/ /g;$name =~ s/\s+$//;$name =~ s/^\s+//;

		$name =~ s/^Candidatus\s(\w+)$/$1/;

		$ncbi_nodes{$tax_id}{name} = $name;

		$names_line_counter++;#print "$names_line_counter\n";


		}else{
		if($line =~ /^(\d+).+scientific name/){die "UNEXPECTED LINE:\n$line\nquitting\n"}
		}

	}

close(NAMES);

print "$names_line_counter names parsed.\n";




}



#####################################################################################################
#
#
#
#####################################################################################################




sub traverse_nodes
{
my $current_node = $_[0];my $current_node_taxstring = $_[1];# $current_node_taxstring to be deprecated

my $child_nodes = $ncbi_nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;
my @child_nodes_array = split(/\t/, $child_nodes);



if($current_node == $starting_node)
	{
	my $rank = $ncbi_nodes{$starting_node}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $current_node_taxstring = substr $current_node_taxstring, 0,1;
	my $originalname = $ncbi_nodes{$starting_node}{name};$originalname =~ s/[\s\t]/_/g;
	my $child_complete_lineage = $ncbi_nodes{$starting_node}{complete_lineage} . "$rank:$originalname ";
	}



foreach my $child(@child_nodes_array)
	{
	unless($current_node == $child)# one the child nodes of the root node (1), is also 1 
	{
	my $rank = $ncbi_nodes{$child}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	
	
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g; ####################### sep2013 .. fixed mar2015
	
	unless($ncbi_nodes{$current_node}{complete_lineage} =~ /\w/)#22oct2015:otherwise where shared taxon is same as basal node, wont be assigned 
		{$ncbi_nodes{$current_node}{complete_lineage}="root:$starting_name_ff "};

	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name
	
	$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;

	#print "storing complete lineage for taxon:$name_string\n";	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;

	$ncbi_tax_number_for_this_species{$name_string}=$child;

#	print "complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";

	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;

	my $name_assignment_to_taxnumber = "";

	if($ignore_subspecies == 1)
		{
		if ( $rank eq "subspecies" || $rank eq "varietas"   || $nodes{$current_node}{complete_lineage} =~ / species:/)# 
			{
			my $parentoriginalname = $ncbi_nodes{$current_node}{name};
			if($parentoriginalname =~ /^[A-Z][a-z]+\s[a-z]+$/)
				{
				$parentoriginalname =~ s/[\s\t]/_/g;
				$name_assignment_to_taxnumber = $parentoriginalname;
			#	print "node:$child named:$originalname rank:$rank appears to be subspecfic\n";
			#	print "\tassigning parent name instead:$parentoriginalname\n\n";
				}else{$name_assignment_to_taxnumber = $originalname}
			}else{
			$name_assignment_to_taxnumber = $originalname
			}

		}else{
		$name_assignment_to_taxnumber = $originalname
		}

	#print "$name_assignment_to_taxnumber $child $rank $child_complete_lineage\n";


		###########################################
		traverse_nodes($child , $child_taxstring);#
		###########################################
	}}


	
}#sub traverse_nodes





#####################################################################################################
#
#
#
#####################################################################################################


