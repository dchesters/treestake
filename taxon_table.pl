#!/usr/bin/perl



#
#
#
#	taxon_table.pl, parses taxonomic databases in preparation for taxonomically-constrained phylogenetic inference	  
#    	Copyright (C) 2019  Douglas Chesters
#
#	This program is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
#
#	This program is distributed in the hope that it will be useful,
#	but WITHOUT ANY WARRANTY; without even the implied warranty of
#	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#	GNU General Public License for more details.
#
#	You should have received a copy of the GNU General Public License
#	along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#	contact address dc0357548934@live.co.uk
#
#
#
#
#
#
##########################################################################################################################################
#
# 
# 
# 	CHANGE LOG
# 
# 	2019-08-08: first release. 
# 		runs on either NCBI database only, or NCBI and ITIS database. 
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
# 	NOTES
#
#	sections of code taken from relational_constraints.pl, part of a process which will be split into at least 3 steps
#	you need to download the NCBI taxonomy database, go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/
#	download the file taxdump.tar.gz, and unzip in the working directory. 
#	the script will read 2 of the resulting files: names.dmp and nodes.dmp.
# 	it will get from these the higher taxonomies to be used for assignment to species given in your tree.
#	only other input is a fasta file of the queries (-seqfile), the script simply identifies query names from this file 
#
# 	ncbi database has tax_id followed by parent_tax_id
# 	
# 	Insecta = 50557, Endopterygota 	= 33392
# 	
# 	
# 	
# 	
# 	
# 	
# 	
# 	
############################################################################################################################################



my $arg_string  = join ' ', @ARGV; # read options specified by user




# globals:
%check_hash4;
$keyfile;
$reference_file;
%species_names;
$starting_node;
$support_cutoff = "NA";
$outgroup = "NULL";
$outfile_prefix;



#####################################
read_command_arguments($arg_string);#
#####################################



################################################################################################

# from parse_ncbi_tax_db ...

$ignore_subspecies			= 1;	# 0=read taxonomies beyond species level 
						#  (e.g. include subspecies, varietas, and unnamed bacteria ranks)
						# if 1, these extra ranks are ignored
						# this applies to full species name, assigned to given ncbi tax number
# globals

%ncbi_nodes;	# ncbi taxonomy structure is recorded in this object
%assign_these_new_terminals_to_taxon;


if($read_NCBI_database == 1)
  {
  ###############
  store_nodes();#
  ###############

  ###################
  parse_namesfile();#
  ###################
  };


################################################################################################




if($read_NCBI_database == 1)
{

open(OUTPUTFILE, ">$outfile_prefix") || die "\ncant open output\n";

$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;$starting_name_ff =~ s/\s+/_/g;
print "name for the starting node ($starting_node)\nwhich has been specified ($starting_name)\n";
print "traversing NCBI taxonomy tree from this node. recording lineage for each taxon.\n\n";
$starting_name =~ s/^(\w).+$/$1/;
$arbitraryID = 0; # used as node IDs in ITIS, later. 

#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################


print "\n\nproccessed NCBI database.\n";


}; # else continue and use ITIS below




if($read_ITIS_database == 1)
	{
	

	# better to do this, for example in Lepidoptera there were an extra 3034 species in ITIS that wernt in NCBI, 
	#	for genera which were in NCBI (ie not even considering genera not in NCBI)
	# one limitation of ITIS is that it only has rank assigned for the main ranks superfam, fam, genus etc.
	# does not have ranks assigned for tribes, which can limit their application.
 

	##############
	store_ITIS();#
	##############


	@lineage_keys = keys %complete_lineage_for_this_species;
	@lineage_keys = sort @lineage_keys;
#	open(LINEAGES, ">$outfile_prefix") || die "";
	foreach my $linege_key(@lineage_keys)
		{print LINEAGES "$linege_key\t$complete_lineage_for_this_species{$linege_key}\n"};
	close LINEAGES;

	print "\nprinted $#lineage_keys lineages to file $outfile_prefix\n\nFIN.\n";

	}else{
	print "user opted not to parse ITIS database\n";
	};

close OUTPUTFILE;


exit;


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
     ***** taxon_table.pl *****
		  
\n";






#######################################################################

# -ITIS $ITIS_database -NCBI $NCBI_names $NCBI_nodes

if($arguments =~ /-ITIS\s+(\S+)\s+(\w+)/)
	{
	$ITIS = $1; $starting_name_ITIS  = $2; $read_ITIS_database = 1;print "using ITIS database\n";
	};

if($arguments =~ /-NCBI\s+(\S+)\s+(\S+)/)
	{
	$NCBI_names = $1; $NCBI_nodes_file = $2;$read_NCBI_database = 1;print "using NCBI database\n";
	}else{
	print "\nerror reading command arguments \n\n";die""
	};

#######################################################################


if($arguments =~ /-outfile\s+(\S+)/)
	{
	$outfile_prefix = $1;
	}else{
	print "\nerror reading command arguments  (-outfile)\n\n";die""
	}

# if($arguments =~ /-treefile\s+(\S+)/)
#	{
#	$treefile = $1;
#	}else{
#	print "\nerror reading command arguments  (-treefile)\n\n";die""#
#	}


# if($arguments =~ /-references\s+(\S+)/)
#	{
#	$reference_file = $1;
#	}else{
#	}


# currently deactivated, it can stop ncbi lineage being read. noticed in Lep study, where lep not included as a rank to constrain.

if($arguments =~ /-constrain_ranks\s+(\S[^\-]+)/)
	{
	$constrain_ranks1 = $1;
	@constrain_ranks = split /\s+/ , $constrain_ranks1;
#	print "\nuser specified ranks at which to make constraints:\n@constrain_ranks\n";
	}else{
	@constrain_ranks = ("order","suborder","infraorder","superfamily","family","subfamily","genus");
#	print "\nyou did not specify ranks at which to make constraints, applying default:\n@constrain_ranks\n";
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




print "\n
user options have been read.....
";



}#sub read_command_arguments



###############################################################################################################################################





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


#	$NCBI_names = $1; $NCBI_nodes = $2;


open(NODES , $NCBI_nodes_file) || die "cant open nodes file, you specified $NCBI_nodes_file
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";

print "\nreading files from NCBI taxonomy database ....
 first nodes.dmp, you specified $NCBI_nodes_file, 
.... ";


my $line_counter=0;
while (my $line = <NODES>)
	{
	if($line =~ /^(\d+)\t[\|]\t([^\t]+)\t[\|]\t([^\t]+)\t[\|]\t/)
		{
		my $tax_id = $1;my $parent_tax_id = $2;my $rank = $3;
		$rank_hash{$rank}++;
#		print "tax_id:$tax_id parent_tax_id:$parent_tax_id rank:$rank\n";

		my $current_rankcode = 0;
		foreach my $rank_test("order","suborder","infraorder","series","superfamily",
			"family","subfamily","tribe","subtribe","genus",
			"subgenus","species group","species","subspecies")
				{
				$current_rankcode++;
				if($rank eq $rank_test){$ncbi_nodes{$tax_id}{rank_code} = $current_rankcode}
				};


			$ncbi_nodes{$tax_id}{rank} = $rank;
			
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

print scalar @all_nodes , " nodes have been read.

";


}




#####################################################################################################
#
#
#
#####################################################################################################


sub parse_namesfile
{

# here just parse the scientific name of each node. ignore synonyms etc

open(NAMES , $NCBI_names) || die "cant open names.dmp ($NCBI_names)
go to ftp://ftp.ncbi.nih.gov/pub/taxonomy/\nget the file taxdump.tar.gz\nunzip and place in working directory\nquitting.\n";


#	$NCBI_names = $1; $NCBI_nodes = $2;


print "\nnames.dmp ($NCBI_names), 
  parsing 'scientific name', ignoring others ... 
";

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
my $parentname = $ncbi_nodes{$current_node}{name};$parentname =~ s/[\s\t]/_/g;

my $child_nodes = $ncbi_nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;
my @child_nodes_array = split(/\t/, $child_nodes);



if($current_node == $starting_node)
	{
	my $rank = $ncbi_nodes{$starting_node}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $current_node_taxstring = substr $current_node_taxstring, 0,1;
	my $originalname = $ncbi_nodes{$starting_node}{name};$originalname =~ s/[\s\t]/_/g;
	my $child_complete_lineage = $ncbi_nodes{$starting_node}{complete_lineage} . "$rank:$originalname ";
	};



foreach my $child(@child_nodes_array)
	{
	unless($current_node == $child)# one the child nodes of the root node (1), is also 1 
	{
	my $rank = $ncbi_nodes{$child}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g; ####################### sep2013 .. fixed mar2015
	my $rankcode = $ncbi_nodes{$child}{rank_code};$rank_codes{$name_string} = $rankcode;

	
	unless($ncbi_nodes{$current_node}{complete_lineage} =~ /\w/)#22oct2015:otherwise where shared taxon is same as basal node, wont be assigned 
		{$ncbi_nodes{$current_node}{complete_lineage}="root:$starting_name_ff "};

	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name


	###############################################################
	# NOV 2018: user decides which ranks are constrained	
	# JUL 2019: deactivated, it can stop lineage being read correctly
	my $constrain_this_node = 0; # print "\n\n";
	foreach my $user_constrain_rank ( @constrain_ranks )
		{
		if($rank eq $user_constrain_rank)
			{$constrain_this_node = 1};
	#	print "rank:$rank user:$user_constrain_rank constrain:$constrain_this_node\n";
		};
#	if($constrain_this_node == 1){
	$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;
#	};
	###############################################################


# forgot about this, explains some imperfect species level trees
#	unless($rank eq "order") # 31 aug 2016: way to stop seqs only given order ident from being constrained.
#		{
#		$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;
#		};


	#print "storing complete lineage for taxon:$name_string\n";	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;

	$ncbi_tax_number_for_this_species{$name_string}=$child;

#	print "complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";

	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;
	$ncbi_taxnumber_for_taxname{$originalname} = $child;
#	print "ncbi_number:$child tax name:$originalname\n";

	my $name_assignment_to_taxnumber = "";my $print_current_node = 0;

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
			$name_assignment_to_taxnumber = $originalname;$print_current_node=1;
			}

		}else{
		$name_assignment_to_taxnumber = $originalname;$print_current_node=1;
		}

	if($print_current_node == 1)
		{

# 	use same convention as ncbi database which has tax_id followed by parent_tax_id

		# store node IDs, need to retreive genus IDs when writing ITIS lines.
		if($rank eq "genus"){$node_ID_of_genus{$name_string} = $child};

		if($child >= $arbitraryID){$arbitraryID = $child + 1};
		

		my $print_node = 	"$child\t" . 		# child node ID
					"$current_node\t" . 	# parent node id
					"$name_string\t" . 	# child tax
					"$parentname\t" . 	# parent tax
					"$rank\n"; 		# child rank

		print OUTPUTFILE $print_node;
		};

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


sub store_ITIS
{

open(ITIS , $ITIS) || die "cant open ITIS database ($ITIS)
 Go to https://www.itis.gov/.
 Search scientific name (e.g. Arthropoda).
 Click on the correct name then download in DwC-A format. 
 Uncompress the file, place in working directory.

  quitting.\n";


print "\nreading ITIS database  ($ITIS), 
";

my $ITIS_line_counter=0;
while (my $line = <ITIS>)
	{
	if($ITIS_line_counter==0)
		{
		print "ITIS file header:$line";
# taxonID	acceptedNameUsageID	parentNameUsageID	scientificName	scientificNameAuthorship	kingdom	phylum	class	order	superfamily	family	genus	subgenus	specificEpithet	infraspecificEpithet	taxonRank	 taxonomicStatus	modified	namePublishedIn	  scientificNameID	verbatimTaxonRank	taxonRemarks	higherClassification
# 0                  1                      2                       3                 4                             5      6       7       8        9            10     11        12              13                14                   15                 16                   17              18                19                     20                      21            22                                                                                        
		}else{

# 1090146		1090135	Atlantea tulita (Dewitz, 1877)	(Dewitz, 1877)	Animalia	Arthropoda	Insecta	Lepidoptera	PapilionoideaNymphalidae	Atlantea		tulita		species	valid	2019-04-27		1090146			Animalia|Bilateria|Protostomia|Ecdysozoa|Arthropoda|Hexapoda|Insecta|Pterygota|Neoptera|Holometabola|Lepidoptera|Papilionoidea|Nymphalidae|Nymphalinae|Melitaeini|Phyciodina|Atlantea
	
		$line =~ s/\n//;$line =~ s/\r//;
		my @split_line = split /\t/ , $line;
		my $scientificName = $split_line[3];my $lineage = $split_line[22];

		if($line =~ /\t$starting_name_ITIS\t|\|$starting_name_ITIS\|/)
			{
		#	print "\nfound $starting_name_ff in line:$line\n";
			#########################################################################################		

			if($scientificName =~ /^([A-Z][a-z]+)\s([a-z]+)/)
				{
				my $genus = $1; my $species = $2; my $binomial = "$genus" . "_" . "$species";
			#	print "species $binomial belongs to $starting_name_ff\n";

				if( $complete_lineage_for_this_species{$genus} =~ /fam/)
					{
					my $complete_lineage = $complete_lineage_for_this_species{$genus};
					if( $complete_lineage_for_this_species{$binomial} =~ /fam/)
						{
						 # print "\tin NCBI\n";	

						}else{# following, contains genus in NCBI, but not this ITIS species


						$complete_lineage_for_this_species{$binomial} = $complete_lineage; $ITIS_species_added{$binomial}++;# print "\tNot in NCBI\n";

						if($node_ID_of_genus{$genus} =~ /\d/)
							{

							my $print_node = 	"$arbitraryID\t" . 		# child node ID
								"$node_ID_of_genus{$genus}\t" . 	# parent node id
								"$species\t" . 	# child tax
								"$genus\t" . 	# parent tax
								"species\n"; 		# child rank
							print OUTPUTFILE $print_node;

							$arbitraryID = $arbitraryID+1;
							};

						};

					}else{
					$ITIS_genera_not_in_NCBI++; # print "genus not in NCBI\n"
					};


				};# 	if($scientificName =~ /^([A-Z][a-z]+\s[a-z]+)$/)

			#########################################################################################		
			};# 		if($line =~ /$starting_name_ff/)


#		print "scientificName:$scientificName "; # scientificName:(Attems, 1922) scie

		};# 	if($ITIS_line_counter==0){}else{


	$ITIS_line_counter++;
	};

my @ITIS_array = keys %ITIS_species_added;

# ? ->	ITIS species added which were not in NCBI:$#ITIS_array


print "
completed reading ITIS database, 
	line_counter:$ITIS_line_counter
	ITIS species added which were not in NCBI:$#ITIS_array
";



};



#####################################################################################################
#
#
#
#####################################################################################################




