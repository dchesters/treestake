#!/usr/bin/perl

###################################################################################################################################
#
#
#
#
#
# 
#	taxon_constraints.pl
#		  
#    	Copyright (C) 2015-2019  Douglas Chesters
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
#
#
#
# 	Default CITATION:
#		Chesters D. 2017. 
#		Construction of a Species-Level Tree-of-Life for the Insects and Utility in Taxonomic Profiling. 
#		Systematic Biology 66 (3): 426–439. 
#
#
#
#
#
# 	for latest version visit https://sourceforge.net/projects/sophi/
#	uses code from parse_clades.pl / bagpipe_phylo.pl (Papadopoulou et al. 2014 / Chesters et al. 2015)
#
#	INPUTS:
#	fasta file  (-seqfile):
# 		The IDs should each be unique in the file, 
#		there should standard characters only used. 
#		strict binomals are preferred, underscore seperated.  thus: [A-Z][a-z]+_[a-z]+
#	user input backbone tree (-treefile):
# 		should have no branchlengths. script should read tree with no node labels, 
#		or tree with nodes labelled with square brackets (format output by sumtrees)
#	
#	
#	
#	
#	
#	Change Log:
#	21 Jan 2015: version 1.001
#	22 oct 2015: modified some bacterial labels when reading ncbi taxonomy db:/^Candidatus\s(\w+)/$1/
#		these are often not used, at least in hte tree  
#		if cant find taxonomy for binomials, tries genus name
#	01 dec 2015: runs without outgroup specified
#		finds taxonomies for genus name when lineage for species is not found
#	15 jan 2017: minor changes to reading sequence file IDs.
#			checks on user's input tree.
#			a few more instructions added.
#			permits species-level data in which taxonomy cannot be determined,
#			these get missing states for all constraints 
#	27 sep 2017: bit more info printed to taxon constraints log
#	01 oct 2017: needed to account for use of same taxonomic name at different ranks of a single lineage,
#			which made a mess of parsing tax names of backbone tree. noticed with collembola,
#			which wernt being constrained due to this name duplicated in ncbi lineage
#			note sometimes these are legitimate (e.g. subgenus often has name of its genus).
#			might want to consider also using rank in this process?
#	02 nov 2018: in addiiton to usual output, prints just the taxon constraints.
#		will probably develop relational constraints seperatly, in newick context.
#		relation constraints will be applied both to raxml and fasttree.
#		and theres two different ways it could be done 
#		depending on comprehensivness of backbone. taxon constraints just one way.
#		thus this script will now focus just on taxon/monophyly constraints for fasttree.
#		and 2 seperate scripts for relational constraints.
#		a couple of new command linbe options:	
#			-ignore_IDs_with_no_tax_assignment and -print_taxon_constraints_only
#	23 dec 2018: changes made to run better on backbone tree in which there are a 
#		 lot of tips for which taxonomic information not available in NCBI DB	
#	25 dec 2018: improved reading of rank order in fasta IDs
#		renamed read_deep_level_constraints.pl to taxon_constraints.pl
#	19 Jul 2019: code for reading taxonomic database is moved to a seperate script
#		code cleaning and annotation.
#
#
#
#
#
#
#
#
#
####################################################################################################
#
#
#
#	REFERENCES
#
# 	Chesters et al. (2015). 
#		A DNA Barcoding system integrating multigene sequence data. 
#		Methods Ecol Evol, 6: 930–937.
#
#	Chesters (2017). 
#		Construction of a Species-Level Tree-of-Life for the Insects and Utility in Taxonomic Profiling. 
#		Systematic Biology 66 (3): 426–439. 
#
#	Papadopoulou et al. (2015). 
#		Automated DNA-based plant identification for large-scale biodiversity assessment.
#		Mol Ecol Resour. 15(1):136-52.
#
#
#	
#	
#	
##########################################################################################################################################


# back to the drawing board.
# as opposed to selecting a shared taxon, there would often be a few ranks above that that node would represent,
# need to draw these on, then select a criteria for assigning, mybe there would be ambiguity.
# 
# 
# second type of constraint made, monophyly of the assigned name, 
# at first instance root->tip, if contains no other taxa, and >1 of the taxa of rank below.
# 
#	Insecta = 50557, Endopterygota 	= 33392




# 
# 
# 
# SUBS:
# read_command_arguments, store_nodes, parse_namesfile, traverse_nodes, read_fas, record_tree2
# traverse_backbone_tree_and_infer_taxonomic_node_labels {get_terminals_from_this_node2, get_shared_taxa, get_all_taxa}
# traverse_backbone_tree_and_collapse
# 
# 
# 
# 
# 


my $arg_string  = join ' ', @ARGV;



# additional user options:

$verbose 		= 1; # 0 simple screen output, 1=verbose, 2=ulta-verbose

#$tree_ID_format = 2;# i think 1=Order_Genus, 2= Order_Genus_species


$process_backbone_tree_terminal_IDs = 0;
# if == 1, then process ids in backbone tree, e.g. Coleoptera_Apatides_fortis to Apatides_fortis


$fasttree_format_constraints_file 	= "ftcons"; # output.
$fasttree_constraints_tabdelim 		= "relational_constraints_tabdelim";


##########################################################################################################################################


# globals:
%check_hash4;

$treefile;
$keyfile;
$fas_file;
$reference_file;
$output_filename;# 	= "$treefile.diet_clades";
%species_names;
$starting_node;
$support_cutoff = "NA";
$outgroup = "NULL";

#####################################
read_command_arguments($arg_string);#
#####################################

$output_filename = "$fas_file.$treefile";


###############################################################################################################################################



$ignore_subspecies			= 1;	# 0=read taxonomies beyond species level (e.g. include subspecies, varietas, and unnamed bacteria ranks)
						# if 1, these extra ranks are ignored
						# this applies to full species name, assigned to given ncbi tax number
# globals

%ncbi_nodes;	# ncbi taxonomy structure is recorded in this object

%assign_these_new_terminals_to_taxon;
# sub store_tax_heirarchy stores the following items:
# 	ncbi_nodes{}{rank_code}; $ncbi_nodes{}{rank}; $ncbi_nodes{}{parent}; 
#	$ncbi_nodes{}{child_nodes}; $ncbi_nodes{}{name}

#######################
store_tax_heirarchy();#
#######################


#################################################

# traverse nodes starts at first parent node in the input file, which should be the root.
# the input file was made with a root node user specified.


# sub traverse_nodes uses:
#  ncbi_nodes{}{child_nodes}; $ncbi_nodes{}{rank};$ncbi_nodes{}{name};$ncbi_nodes{}{rank_code};
# and writes:
#  complete_lineage_for_this_species ; ncbi_taxnumber_for_taxname	

$starting_name = $ncbi_nodes{$starting_node}{name};$starting_name_ff = $starting_name;$starting_name_ff =~ s/\s+/_/g;
print "name for the starting node ($starting_node)\nwhich has been specified ($starting_name)\n";
$starting_name =~ s/^(\w).+$/$1/;

#################################################
traverse_nodes($starting_node , $starting_name);#
#################################################


my @ranks_array = keys %how_many;
print "heirarchy traversed:
	nodes_traversed:$nodes_traversed
	ranks encountered:$#ranks_array
";




# hack, some missing species
# in other scripts, the presence of these when they actually are in the ncbi tax DB (e.g. newer db is used), causes error.

#$complete_lineage_for_this_species{"Lemyra_melli"} 		= $complete_lineage_for_this_species{"Lemyra_jankowskii"};
#$complete_lineage_for_this_species{"Promalactis_suzukiella"} 	= $complete_lineage_for_this_species{"Promalactis_jezonica"};
#$complete_lineage_for_this_species{"Hylaeus_dilatatus"} 		= $complete_lineage_for_this_species{"Hylaeus_difficilis"};
#$complete_lineage_for_this_species{"Lethe_dura"} 			= $complete_lineage_for_this_species{"Lethe_diana"};
#$complete_lineage_for_this_species{"Eucryptorrhynchus_brandti"} 	= $complete_lineage_for_this_species{"Eucryptorhynchus_brandti"};




###############################################################################################################################################



my %terminals = ();
my $root_identity = "";
my $date = localtime time;

%count_tax_in_new_data;
%query_IDs;	# the fasta file contains names of new entries to be placed onto backbone tree. 
		# each are stored in the hash %query_IDs.


############
read_fas();#
############




#########################################################################################################

	# store backbone phylogeny.

%nodes;		# user phylogeny stored in this hash
$root_node;	# this global will keep the value of the last in the loop, which will be the root
%count_taxa_represented_in_terminals;	# will use this to determine whether taxa are monophyletic.
$print_internal_nodes=0;
# other items stored here: complete_lineage_of_terminal tax_encounted_in_current_lineage 


# read backbone tree into hash. node identities for hash key, then parent/child identities as hash entry
open(BACKBONE_LOG, ">traverse_backbone_LOG") || die "";


#########################
record_tree2($treefile);#
#########################

print "
your backbone tree has been read, TIPS:$count_the_terminal_nodes INTERNAL NODES:$print_internal_nodes.
\n";

#########################################################################################################






#########################################################################################################

	# traverse nodes of the user genome tree, 
	# for use in constraining treesearch in the barcode data


$new_newick_string = "($root_node)";
$new_newick_string2 = "($root_node)";
$newick_print_string = "($root_node)";

%all_taxa_on_backbone_tree;

%monophylys; 	# following sub infers monophylys for each node, 
#		by getting all taxa from node and checking that i) all terminals have it ii) its not found elsewhere
#		first stored as %monophylys then converted to array @mps
#		this information is used in for fasttree constraint table only	


print "\n
calling sub traverse_backbone_tree_and_infer_taxonomic_node_labels, 
which will traverse nodes of the backbone tree, infer taxonomic names for each node
";


# relativly complicated sub.
# does monophyly test for fasttree constraints, and does taxonomic inference for nodes which are assigned to new_newick_string.

#########################################################################
traverse_backbone_tree_and_infer_taxonomic_node_labels($root_node , 0);#
#########################################################################


close BACKBONE_LOG;
print "finished writing backbone analysis log named traverse_backbone_LOG\n";


open(LOG8, ">taxon_constraints_INFO") || die "";
print LOG8 "following is list of all taxa found on backbone tree:\n";
my @alltax9 = keys %all_taxa_on_backbone_tree;@alltax9 = sort @alltax9;

#########################################################################################################



# count_taxa_represented_in_terminals



my @mps = keys %monophylys;@mps = sort @mps;
foreach my $mp(@alltax9)
	{
	my $count = $count_taxa_represented_in_terminals{$mp};
	if($count >= 2)
		{
	
	if(exists($monophylys{$mp}))
		{
		print LOG8 "$mp\t$count\tIS monophyletic\n";
		}else{
		print LOG8 "$mp\t$count\tNOT monophyletic\n";
		}
		};
}
print LOG8 "\n$#mps tax monophylys found on backbone tree, these will be set for next tree search if relevent 
, in addition to the other things\n";

close LOG8;


# taxonomic node labels assigned to backbone tree.
# now traverse tree again and decide at which point to stop (collpase all decendents)
# this must be done where nodes are reached at which not all taxon represented in the fasta file 
# are included as decendents.
# also might need to collapse at where polyphyletic taxa are reached.



print "\ncalling sub traverse_backbone_tree_and_collapse
which will do something.
";

#########################################################################
traverse_backbone_tree_and_collapse($root_node);#
#########################################################################




	# these are RELATIONAL constraints 
my @ft_constraints_again = keys  %ft_constraints;
@ft_constraints_again = sort @ft_constraints_again;
open(CONS_TABLE , ">$fasttree_constraints_tabdelim") || die "\nerror 258\n";
foreach my $constraint(@ft_constraints_again)
	{
	#print "constraint:($constraint)\n";
	$constraint =~ s/\n/\t/g;
	$constraintno++;
	print CONS_TABLE "constraint:$constraintno\t$constraint\n";
	}
close CONS_TABLE;





#print "
#new_newick_string:$new_newick_string\nnew_newick_string2:$new_newick_string2\n
#"; 

#open(NEWICK1 , ">newick_string1") || die "\nerror 193.\n";
#print NEWICK1 "$new_newick_string\n";
#close NEWICK1;
#open(NEWICK2 , ">newick_string2") || die "\nerror 193.\n";#
#print NEWICK2 "$new_newick_string2\n";
#close NEWICK2;








# process newick string by removing the terminals, then ready for adding the list of new taxa
# %query_IDs contains all the fasta IDs (should be standardized to binomials) in the species-dense data

my @new_terminals = keys %query_IDs;@new_terminals = sort @new_terminals;

my @all_taxa = keys %all_taxa_on_backbone_tree;@all_taxa = sort @all_taxa;
print "
all taxa:$#all_taxa
";


open(OUT5 , ">$fasttree_format_constraints_file") || die "\nerrro\n";
print OUT5 " " , scalar @new_terminals , " " , scalar @ft_constraints_again + scalar @mps + 1, "\n";

# fasttree constraint file.
# 5 2
#A        00
#B        0-
#C        10
#D        11
#E        -1

open(OUT6, ">taxon_constraints_only") || die "";
open(OUT7, ">taxon_constraints.fas") || die "";

$char_weight = 10;


my $coutn = 0;my $warning_count=0;


foreach my $new(@new_terminals)
	{

	$coutn++;if($coutn=~ /000$/){print "$coutn of $#new_terminals\n"};
	my $genus = $new;$genus =~ s/^([A-Z][a-z]+)_.+/$1/;my $lineage;

	if(exists($complete_lineage_of_fasta_ID{$new}))
		{
		# complete_lineage_for_this_species, not just species, 
		# but has lineage for every taxon in ncbi, only mod is space->underscore
		$lineage = $complete_lineage_for_this_species{$new};

		}elsif(exists($complete_lineage_for_this_species{$genus}))
			{
			$lineage = $complete_lineage_for_this_species{$genus};
			}else{
			$warning_count++;
			if($warning_count == 50){print "\ntoo many warnings.\n"};
			if($warning_count >= 51){}else{
				print "WARNING ... dont know what this is:$new, cannot constrain it\n";
				}
			};
	#print "complete lineage:$lineage\n";
	#print "\nnew species to set state of constraints tree:($new) genus:($genus) lineage:($lineage)\n";

	print OUT5 "$new        ";print OUT6 "$new        ";print OUT7 ">$new\n";


	# 'relational constraints'

	unless($print_taxon_constraints_only == 1)
	{
	my $constraintno =0;
	foreach my $constraint(@ft_constraints_again)
		{
		$constraint =~ s/[\n\t]$//;# perl behaviour i dont understand. 
				#changing /n /t in a variable earlier, also changed it in the hash it was read from
		$constraintno++;
	#	print "constraint $constraintno\n($constraint)\n";
		my @split_components  = split /[\n\t]/, $constraint;
		unless($#split_components >= 1){die "\nerror. why a constraint with only one taxon\n"};
		my $state_to_assign = "-";
		my $state_A_assigned=0;my $state_B_assigned=0;

		foreach my $i(@split_components)
			{
			# go through each taxa in the current constraint, see if the new sp belongs to any of these

			if($i =~ /(\w+)\:(.)/)# why was it (.+) ???
				{
				my $tax = $1; my $state=$2;
				
				#print "\ttax:$tax state:$state\n";
				if($lineage =~ /$tax\s/)
					{
					$state_to_assign = $state;
					# this wont work:
					#if($state eq "0"){$state_A_assigned=1};if($state eq "1"){$state_B_assigned=1};
					}else{};
				}else{
				die "\nerror 319\n";
				}
			}


		if($unconstrained_IDs{$new} == 1){$state_to_assign = "-"};

		# there will be no MP constraints in which outgroup is not excluded:
		if($new eq $outgroup){$state_to_assign = "0"};


		#print "\t\tstate_to_assign:$state_to_assign\n";
		print OUT5 "$state_to_assign";
		#die;

		# for certain constraints, e.g. order level, all new seqs should have been assigned a state. check.
#		if($constraint =~ /Coleoptera:\d[\n\t]Diptera:\d[\n\t]Hymenoptera:\d[\n\t]Siphonaptera:\d[\n\t]Mecoptera:\d[\n\t]Strepsiptera:\d[\n\t]Amphiesmenoptera:\d/)
#			{#die "\nfound test constrinat... \n";
#			if($state_to_assign eq "-"){die "\nno state assigned.\n"}
#			};

		}# foreach my $constraint(@ft_constraints_again)

	}; # 	unless($print_taxon_constraints_only == 1)




	foreach my $tax(@mps)	# monophyly taxon group constraints.
		{
		#print "tax($tax)\n";

		if($unconstrained_IDs{$new} == 1) # these are members with no lineage info
			{
			print OUT5 "-";print OUT6 "-";print OUT7 "-" x $char_weight;
			}else{
			$count_constraint_each_member{$new}++;
			if($lineage =~ /$tax\s/)
				{
				$record_taxon_constraints_ouput_MPa{$tax}=1;
				print OUT5 "1";print OUT6 "1";print OUT7 "1" x $char_weight;#print "\tcurrent tax belongs to this monophyletic group\n";
				$count_seqs_in_each_monophyly_taxon_constraint{$tax}++;
				}else{
				$record_taxon_constraints_ouput_MPb{$tax}=1;
				print OUT5 "0";print OUT6 "0";print OUT7 "0" x $char_weight;#print "\tcurrent tax doesnt belong to this monophyletic group\n";
				};
			}

		};


	# one extra constraint. although this doesnt seem to work on its own ....

	if($new eq $outgroup)
		{
		print OUT5 "0";print OUT6 "0";print OUT7 "0" x $char_weight;$outgroup_found++;
		}else{
		print OUT5 "1";print OUT6 "1";print OUT7 "1" x $char_weight;
		};

print OUT5 "\n";print OUT6 "\n";print OUT7 "\n";


	};#$new(@new_terminals)


open(OUT57 , ">taxon_constraints_output_LOG") || die "\nerrro\n";
print OUT57 "backbone_taxon	is_taxon_applied_to_new_matrix	count_seqs_constrained\n";
foreach my $tax(@mps)
	{
	my $count = $count_seqs_in_each_monophyly_taxon_constraint{$tax};
	my $out_string = "$tax\t";

	if($record_taxon_constraints_ouput_MPa{$tax}==1 && 
		$record_taxon_constraints_ouput_MPb{$tax}==1)
		{
		print "backbone monophyly $tax\tapplied to new matrix \n";$backbone_monophylys_applied_to_new_matrix++;
		$out_string .= "1\t";
		}else{
		print  "backbone monophyly $tax\tnot applied to new matrix \n";$backbone_monophylys_not_applied_to_new_matrix++;
		$out_string .= "0\t";
		}
	$out_string .= "$count\n";
	print OUT57 $out_string;

	}
close OUT57;


my $count_taxa_found_on_backbone = scalar @alltax9;
my $count_monophyletic_taxa_on_backbone = scalar @mps;
my $count_members_of_new_matrix = scalar @new_terminals;
my $number_of_new_members_with_constraints = scalar keys %count_constraint_each_member;


print "
count_taxa_found_on_backbone:$count_taxa_found_on_backbone
count_monophyletic_taxa_on_backbone:$count_monophyletic_taxa_on_backbone 
backbone_monophylys_applied_to_new_matrix:$backbone_monophylys_applied_to_new_matrix
backbone_monophylys_not_applied_to_new_matrix:$backbone_monophylys_not_applied_to_new_matrix

count_members_of_new_matrix:$count_members_of_new_matrix
fasta_IDs_ignored_cause_no_lineage:$fasta_IDs_ignored_cause_no_lineage
 (you set ignore_IDs_with_no_tax_assignment:$ignore_IDs_with_no_tax_assignment)
number_of_new_members_with_constraints:$number_of_new_members_with_constraints
";

#unless($outgroup_found == 1){die "\nerror assigning outgroup constraint.\n"}


close OUT5;close OUT6;close OUT7;
print "\nscript completed. END.\n";
exit;



print "out of $#new_terminals new seqs, $assigned22 can be unambiguously placed on the backbone tree, 
$not_assigned22 cannot

\n";

if($no_lineage_found >= 1)
	{
	print "no lineage found for $no_lineage_found, therefor also cant assign\n";
	}else{
	print "lineage was found for all ... something\n ";
	};



%store_node_assigned_newIDs;

$newick5 = "($root_node)";


# sub builds up $newick5

#############################################
get_terminals_from_this_node30($root_node);#
#############################################

$newick5 =~ s/^\((.+)\)$/$1/;

$raxml_constraint_tree = "$treefile.rax_constraint";

open(NEWICK5 , ">$raxml_constraint_tree") || die "\nerror 193.\n";
print NEWICK5 "$newick5;\n";
#print "$newick5\n";
close NEWICK5;




# tnt: given up trying to get it running.

open(TNT_FILE_HANDLE , ">$fas_file.$treefile.tnt") || die "\nerror 366\n";

%assigned_already;

# tnt format, 
#	# constraint where 2 taxa 'float'
#	force + [ a b c  (d e) ]
#tnt64 mxram 2000, rseed 1,p default-Second.nex timeout 0:20:00, echo=, taxname=, force + [17 18 19 20], constrain =,
#
#$tnt_command_string

$tnt_command_string = "tnt64 'mxram 2000, rseed 1,p tnt_input.nex timeout 0:20:00, echo=, taxname=, ";


# for fasttree constraint file (uses format similar to phylip)
$constraint_number;
%phylip_constraints;

# will be stored in hash under keys 1..$constraint_number
# ids you can get from my @all_new_IDs = keys %query_IDs;


# sub builds list of commands in $tnt_command_string

#############################################
get_terminals_from_this_node31($root_node);#
#############################################



# fasttree constraint file.
# 5 2
#A        00
#B        0-
#C        10
#D        11
#E        -1



my @all_new_IDs = keys %query_IDs;@all_new_IDs = sort @all_new_IDs;

print "\nmaking constraint file for fasttree.
dimensions: taxa:" , scalar @all_new_IDs,
"
constraints:$constraint_number
\n";

open(OUT_FT , ">$output_filename.ft") || die "\nerror 421\n";
print OUT_FT " " , scalar @all_new_IDs , " $constraint_number\n";

foreach my $id(@all_new_IDs)
	{
	print OUT_FT "$id        ";
		
#	for my $constaint_no(1..$constraint_number)
#		{
		my $state = $phylip_constraints{$id};print OUT_FT "$state";
#		}

	print OUT_FT "\n";
	};




#constrain =, outgroup 0, log tnt_logfile,  hold 2000, xmult=replication 100 drift 10 ratchet 10 hold 10, best, majority *, tchoose / ,ttags= , resample boot replications 500 from 0 , ttags, tsave *tnttrees, save*, tsave / ,export *tnt_output, log / ,quit;

$tnt_command_string .= "constrain =, log tnt_logfile,  hold 2000, xmult=replication 100 drift 10 ratchet 10 hold 10, best, majority *, tchoose / ,ttags= , resample boot replications 500 from 0 , ttags, tsave *tnttrees, save*, tsave / ,export *tnt_output, log / ,quit;'";



#print "\nrunnning tnt\n";
#system($tnt_command_string);

print TNT_FILE_HANDLE "$tnt_command_string";
close TNT_FILE_HANDLE;

print "
FIN.
";
exit;



#$new_newick_string =~ s/\)(\w+)/ $assign_these_new_terminals_to_taxon{$1} . ")" . $1 /ge;

		if($new_newick_string =~ /(........Pierini.......)/)
			{print "dolar 1:$1\n"};


my @taxa_to_which_new_species_are_assigned = keys %assign_these_new_terminals_to_taxon;
@taxa_to_which_new_species_are_assigned = sort @taxa_to_which_new_species_are_assigned;
print "\n" , scalar @taxa_to_which_new_species_are_assigned , 
" taxon names on the backbone tree, to which new terminals will be added\n\n";

print  "new newick string:$new_newick_string;\n";


foreach my $taxon(@taxa_to_which_new_species_are_assigned)
	{
	if($taxon =~ /\w/)
	{
	my $new_terminals_for_this = $assign_these_new_terminals_to_taxon{$taxon};
	$new_terminals_for_this =~ s/\,$//;

	# hack to cut down tree for vewiing:
	$cut_down_for_viewing = 1;
	if($cut_down_for_viewing ==1)
		{
		$new_terminals_for_this =~ s/^([\_A-Za-z]+)\,.+\,([\_A-Za-z]+)$/$1 . "," . $2/e;  #################
		};
	print "\ntax:($taxon)\tnew terminals:($new_terminals_for_this)\n";

	if($new_newick_string =~ /(.......[\(\)\.]$taxon[\(\)\.]......)/)
		{print "looks like:$1\n"};

	my $taxname = $taxon;
	$inculde_tax_branch_labels = 0;
unless($inculde_tax_branch_labels ==1){$taxname = ""};

	if($new_newick_string =~ s/\(\)$taxon([\,\(\)])/ "(" . $new_terminals_for_this . ")" . $taxname . $1/e)
		{
			#if($new_newick_string =~ /(........Pierini.......)/)
			#	{print "dolar 1:$1\n"};
		}elsif($new_newick_string =~ s/(\w)\)$taxon([\,\(\)])/ $1 . "," . $new_terminals_for_this . ")" . $taxname . $2/e)
			{
			}elsif($new_newick_string =~ s/\)\)$taxon([\,\(\)])/ ")," .  $new_terminals_for_this . ")" . $taxname . $1/e){#....tata))Cucujiformia)Polyp....
			}else{
#print  "\nERROR. newick string:$new_newick_string;\n";
		
		die "\nerror 266:cant find $taxon on tree\n"
			
			};

unless($inculde_tax_branch_labels ==1)
	{
#ta)))Nymphalidae,((Ant

while($new_newick_string =~ s/([\,\(\)])$taxon([\,\(\)])/$1$2/g){};

	};



	}
	}; # foreach my $taxon(@taxa_to_which_new_species_are_assigned)






# open(NNS, ">new_newick_string");
print NNS "$new_newick_string;\n";
close NNS;



print "\n\n\nend of script\n";
exit;





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub record_tree2
{


my $tree1= "";
open(IN, $treefile) || die "\n\nerror 1408 cant open file $treefile\n\n";
while (my $line = <IN>)
	{
	if($line =~ /(.+\(.+)/){$tree1 = $1}
	}
close(IN);

print "\nlooking at backbone tree from file named:$treefile\n";

print BACKBONE_LOG "\nparsing tree from file named:$treefile\n";

$tree1 =~ s/ //g;



$tree_parse = 1; 	
if($tree_parse ==1)
	{
	# remove branchlengths, scientific notation, incl negative values for distance trees. example: -8.906e-05
	if($tree1 =~ s/\:\-*\d+\.\d+[eE]\-\d+([\(\)\,])/$1/g){$branchlengths_rm++};
	# remove regular branchlengths: 0.02048
	if($tree1 =~ s/\:\-*\d+\.\d+//g){$branchlengths_rm++};
	# remove 0 length branchlengths
	if($tree1 =~ s/\:\d+//g){$branchlengths_rm++};
	};

if($branchlengths_rm >= 1)
	{
	print "\nit is not possible to constrain branch lengths. your input tree seems to have branch lengths.
	these should be removed prior. this script has tried to do this anyway.\n";
	}else{
	print "\nseems no branchlengths in you backbone tree, good.\n";
	};


if($tree1 =~ /\:/){die "\nstill a colon in your input tree. 
this is usually a branch length, for some reason not sucessfully removed here.
otherwise it might be part of your labels. labels should use only simple characters,
certainly no parentheses nor colons,
\n"};





my $newick_string 	= $tree1;
my $interal_node	= 0;



# new newick parser ....

#while ($newick_string =~ s/\(([^\(\)]+)\)(\d*)/INTERNAL_NODE_$interal_node/) # rooted bipart raxml tree, the 2 last nodes will have no boot values
while ($newick_string =~ s/\(([^\(\)]+)\)([0-9\.]*)/INTERNAL_NODE_$interal_node/) # processed sumtrees output from a load of raxml boots
	{
	my $node = $1;my $boot = $2; my $nodeID = "INTERNAL_NODE_$interal_node"; #print "nodeID:$nodeID node:$node\n";

	# found a pair of parentheses (with no parentheses inside).
	# this string is taken, and replaced with a single node (INTERNAL_NODE_\d)
	# node is numbered (\d in variable $interal_node), 
	# the number being incremented at the end of the loop
	# find what is at the current node (decendents), by splitting at the commas

	if($boot =~ /\d/){$nodes{$nodeID}{support} = $boot}else{$nodes{$nodeID}{support} = 100};#print "\tnodeID:$nodeID boot:$boot\n";

	my @child_nodes = split /\,/ , $node;
	$child_counts{$nodeID} = $#child_nodes;



	for $i(0 .. $#child_nodes)
		{
		#print "$child_nodes[$i]\n";
		my $bl = "NA"; 	
		if($child_nodes[$i] =~ /\:(.+)/)	# branchlength found
			{$bl = $1};
		$child_nodes[$i] =~ s/\:(.+)//;
		#print "node:$child_nodes[$i]\tbl:$bl\n ";

		unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
			{
			if($remove_accessions_from_reference_species == 1)
				{$child_nodes[$i] =~ s/^([A-Z][a-z]+_[a-z]+)_.+$/$1/}
			if($process_backbone_tree_terminal_IDs == 1)
				{$child_nodes[$i] =~ s/.+_([A-Z][a-z]+_[a-z]+)$/$1/}
			};

		# record each decendent of the current node, in a hash nodes{ the current node }{ child number }

		$nodes{$nodeID}{$i} 			= $child_nodes[$i];
		$nodes{$child_nodes[$i]}{parent} 	= $nodeID;
		$nodes{$child_nodes[$i]}{branchlength} 	= $bl;

		# and record whether the current node is a terminal

		unless($child_nodes[$i] =~ /INTERNAL_NODE_/)
			{
			$terminals{$child_nodes[$i]} =1;$count_the_terminal_nodes++;
			my $genus_name = $child_nodes[$i];$genus_name =~ s/^([A-Z][a-z]+)_.+/$1/;
			my $complete_lineage = $complete_lineage_for_this_species{$child_nodes[$i]};
			unless($complete_lineage =~ /\w/)
				{$complete_lineage= $complete_lineage_for_this_species{$genus_name}};
			if( $complete_lineage =~ /\w\s/)
				{
				# here store lineage, instead of cumbersome retrivals each time needed later.
				$complete_lineage_of_terminal{$child_nodes[$i]} = $complete_lineage;
				}else{
				print "\nERROR, cant find a taxonomy for a terminal of your backbone tree; 
				note, in backbone tree only Genus_sp or Genus_species are permitted. no higher taxa.
				details: node $interal_node, child_nodes[i]:($child_nodes[$i]) genus_name:($genus_name) complete_lineage:$complete_lineage\n";
				print "1) check spelling is same as NCBI taxonomy. 2) Try removing the offending terminal from your backbone tree, then re-run. \nerror 756. quitting.\n\n";
			#	die "";
				};



			unless($child_nodes[$i] =~ /^[\w\_]+$/){die "\nerror 870. unexpected node label (re-name of possible?):$child_nodes[$i]\n"};

			print BACKBONE_LOG "child_nodes[i]:$child_nodes[$i] complete_lineage:$complete_lineage\n";

			# need to count how many terminals there are for each taxa,
			# for a given taxa, if all of them are decended from one node
			# it means they are monophlyetic, and will be constrained.

			# argh! ncbi taxonomic heirachy, duplicated taxa:
			# root:Hexapoda class:Collembola no_rank:Collembola suborder:Symphypleona
			# will need to be accounted for ... 
			# some times it is legit, such as genus and subgenus,
			# still, messes up this thing

			my %tax_encounted_in_current_lineage = ();
			while($complete_lineage =~ s/^([^:]+):(\w+)//)
				{
				my $current_rank = $1;my $current_taxname = $2;
				if($tax_encounted_in_current_lineage{$current_taxname} == 1)
					{
					print "tax name $current_taxname is duplicated in ncbi lineage\n";print BACKBONE_LOG "tax name $current_taxname is duplicated in ncbi lineage\n";
					}else{
					$count_taxa_represented_in_terminals{$current_taxname}++; # will use this to determine whether taxa are monophyletic.
					};
				$tax_encounted_in_current_lineage{$current_taxname} = 1;
				};



			}#unless($child_nodes[$i] =~ /INTERNAL_NODE_/)

		}
	#print "node:$interal_node\n\tchild nodes:@child_nodes\n";

	$root_node = $nodeID; # this is a global and will keep the value in the final loop, so the root can be identified later.

	$interal_node++;$print_internal_nodes++;
	}



#print "newick_string:$newick_string\n";die;


unless($interal_node >= 2){die "\nerror reading your phylogeny.\n"}


}#sub record_tree2




#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




sub read_fas
{

open(IN, $fas_file) || die "\nerror ... fasta file you have given ($fas_file) cannot be opened. 
make sure it is named correctly and in the working directory. quitting\n";


print "
opened fasta file $fas_file, reading ...\n";

my $count_diet_IDs = 0;
my $error_quit = 0;

while (my $line = <IN>)
	{
	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^>(.+)/)
		{
		my $id = $1;#print "id:$id\n";
		if($id =~ /^[a-z]+_[a-z]+$/i)
			{
			$count_binomials++;
			}else{ # elsif($id =~ /^[A-Z][a-z]+_[a-z]+_[a-z]+/)
			# 20170115: i dont recall the rationale here. have allowed non binomials for now
				$extended_ids++;
				$warnings_encountered++;

				if($warnings_encountered <= 20)
					{
					print "fasta id:$id. Warning, strict binomial format and one entry per species is preferred." , 
					" effects uncertain. continuing anyway ...\n";
					}else{
					if($warnings_encountered == 21)
						{
						print "too many warnings to print:$warnings_encountered, violation of strict binomial format:fasta id:$id\n";
						};
					};

			};

	#	$index_of_fastaID{$id} = $count_diet_IDs; # used in discontiued subroutine
		$count_diet_IDs++;
		

		# get all tax names of this species, then record +1 for each.
		my $complete_lineage = $complete_lineage_for_this_species{$id};
		unless($complete_lineage =~ /\:\w/)
			{
			my $genusname = $id;$genusname =~ s/^([A-Z][a-z]+)_[a-z]+.*/$1/;
			$complete_lineage = $complete_lineage_for_this_species{$genusname};
			};
		unless($complete_lineage =~ /\:\w/)
			{# doesnt work:
		#	if($id =~ /^(Archaeognatha)|^(Zygentoma)|^(Odonata)|^(Ephemeroptera)|^(Plecoptera)|^(Orthoptera)|^(Phasmatodea)|^(Isoptera)|^(Hemiptera)|^(Hymenoptera)|^(Raphidioptera)|^(Megaloptera)|^(Neuroptera)|^(Coleoptera)|^(Diptera)|^(Lepidoptera)/)
			if($id =~ /^([A-Z][a-z]{2,9}tera)/)
				{
				my $order_name = $1;
				$complete_lineage = $complete_lineage_for_this_species{$order_name};
			#	print "order_name:$order_name complete_lineage:$complete_lineage\n";
				};
			};
		if($complete_lineage =~ /\:\w/)
			{
			$query_IDs{$id}= 1;
			$complete_lineage_of_fasta_ID{$id}= $complete_lineage;
			}else{
			$warnings_encountered_b++;
			if($warnings_encountered_b <= 20)
				{
				print "WARNING 902 ... cant find lineage for fasta id ($id) nor first name ($genusname)," , " which presuming genus\n" , "cannot constrain this\n\n";
				}else{
				if($warnings_encountered_b == 20){print "too many warnings to print:$warnings_encountered_b, cant find lineage for fasta id ($id) nor first name ($genusname)\n"};
				};
			$unconstrained_IDs{$id} = 1;$error_quit++;	
			if($ignore_IDs_with_no_tax_assignment == 1)
				{
				$fasta_IDs_ignored_cause_no_lineage++;
				}else{
				$query_IDs{$id}= 1
				};	
			};		


		while ($complete_lineage =~ s/ \w+\:(\w+) / /)
			{my $tax = $1;$count_tax_in_new_data{$tax}++; #print "tax:$tax\n";
			};		


		}else{
		if($line =~ /^>/){die "\nerror 646.\n"}
		}
	} # while (my $line = <IN>)
close(IN);

if($warnings_encountered_b >= 100)
	{
	print "warning, no taxonomy found for $warnings_encountered_b, cant constrain those.\n";
	};

if($error_quit >= 1)
	{
	print "\n
	ERROR 1001. cannot find taxonomies for $error_quit entries in fasta file (out of total $count_diet_IDs).
";
#	if small number, please remove from your fasta file and re-run.
#	these $error_quit should be indicated above.
#	\n";
#	die "";
	};



print "
count_binomials:$count_binomials 
extended_ids:$extended_ids
";


$autodetect_fasta_IDs = 0;
if($count_diet_IDs == $count_binomials)
	{$autodetect_fasta_IDs = "binomial";print "\nautodetected strict binomial used in your fasta file. good!\n"};
if($count_diet_IDs == $extended_ids)
	{$autodetect_fasta_IDs = "extended";print "\nwarning. you are not using strict binomials in your fasta file, not guarunteed to work\n"};

my @check_species_filtered = keys %query_IDs;

if($count_diet_IDs == scalar @check_species_filtered)
	{
	print "\npass check, all strings in fasta file are unique (species filtered)\n"
	}else{
	print "\nwarning. duplicate IDs found in your fasta file.\n"
	}

print "\n*** $count_diet_IDs fasta IDs *** in your fasta file ($fas_file). ...\n\n";

unless($count_diet_IDs >= 1)
	{die "\nerror ... no fasta ID's found in the file ($fas_file). quitting.\n"}




}#sub read_fas





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub get_terminals_from_this_node
{
my $next = $_[0];my $sum_branchlength = $_[1];


#print "\nNEW NODE ($next)";

my $child_nodes = $child_counts{$next};
my @next1 = ();
my @non_terminal_child_nodes = ();
my @child_nodes_test = ();


my $all_terminals_from_both = "";
my @taxon_for_child_nodes = ();
my $tax_for_all_child_nodes = "";
my $all_child_nodes_are_internal =1;

my $tax_for_parent = "";
					# do this even if child node is terminal
if($next =~ /INTERNAL_NODE_\d/) 	# && $all_child_nodes_are_internal == 1)
	{
	$all_terminals_from_this_node ="";
	get_terminals_from_this_node2($next, 0);
	$tax_for_parent = get_shared_taxonomic_name($all_terminals_from_this_node);
	#print "current node ($next) is internal, shared tax name ($tax_for_parent)\n";
	};


for $i(0 .. $child_nodes)
	{
	$count_terminal_nodes = 0;
	$all_terminals_from_this_node ="";# for current child node only

	if($nodes{$next}{$i} =~ /INTERNAL_NODE_\d/)
		{
		# if current child is internal node, go and get the decending terminals
			# this sub does this for terminal nodes only:
			# $all_terminals_from_this_node .= "$test\t";
			# gets binomial names
		#print "line561. getting terminal binomials for child node $nodes{$next}{$i}\n";
		get_terminals_from_this_node2($nodes{$next}{$i}, 0);
		
		}else{
		# otherwise, just append the node name itself
		my $established_name = $nodes{$next}{$i};
		$all_terminals_from_this_node .= "$established_name\t";
		$all_child_nodes_are_internal =0
		};



		$all_terminals_from_both .= $all_terminals_from_this_node;
		$all_terminals_from_this_node =~ s/\t$//;

				##############################################################
		my $tax = 	get_shared_taxonomic_name($all_terminals_from_this_node);#
				##############################################################

		push @taxon_for_child_nodes, $tax;
		#print "child $i of node $next has ID $nodes{$next}{$i}. $count_terminal_nodestermianls got, shared tax:$tax.\n";

		# here assess if child node is same taxonomy as node.
		# if so, this node is collapsed.
		# this is done by stepping ahead to getting the two grandchild nodes, 
		# then write both these in instead of the child node.


		if($tax eq $tax_for_parent)
			{
			print "this child has same taxon ($tax,$tax_for_parent), collapsing.\n";
			my $grandchild_nodes = $child_counts{$nodes{$next}{$i}};
			for $j(0 .. $grandchild_nodes)
				{
				my $push_node = $nodes{$nodes{$next}{$i}}{$j};
				print "push grandchild to array j:$j, push_node:$push_node\n";
				push @next1 , $push_node;
				if($push_node =~ /INTERNAL_NODE_\d/)	# for makeing the newick string with termainals for new data,
					{
					push @non_terminal_child_nodes, $push_node;
					push @child_nodes_test, $push_node;
					}elsif# first need to remove existing terminals
					($push_node =~ /\w/)
						{
						push @child_nodes_test, $push_node;
						}else{die "\n\nerror 788\n"}	
				};
			}else{
			print "this child has diff taxon ($tax,$tax_for_parent)\n";
			my $push_node = $nodes{$next}{$i};
			push @next1 , $push_node;	
			if($push_node =~ /INTERNAL_NODE_\d/)
				{push @non_terminal_child_nodes, $push_node};
			push @child_nodes_test, $push_node;
			};

		if($i == $child_nodes)
			{
			$all_terminals_from_both=~ s/\t$//;
			print "line793, calling sub get_shared_taxonomic_name\n";
			$tax_for_all_child_nodes = get_shared_taxonomic_name($all_terminals_from_both);
			unless($tax =~ /\w/){die "\nno tax name retreived for this set of terminals\n"}
			};



	}



	# this is the default way to building up the newick string with new nodes.
	# it is overwritten later in situations where the node needs to be collapsed.


# this is collapsed where appropriate, but includes terminal nodes. 
# printed to an additional newick string so user can see better how tree has been collapsed 
my $join_child_nodes_test = join ',', @child_nodes_test;

# removed terminals:
my $join_child_nodes = join ',',@non_terminal_child_nodes;
my $swap_string ="";

#if($#non_terminal_child_nodes <= 0)
#	{
#	print "no decendents. dead end.\n";
##$swap_string = "$tax_for_parent";
#$swap_string = "($join_child_nodes)$tax_for_parent";
#	}else{
$swap_string = "($join_child_nodes)$tax_for_parent";
my $swap_string_test = "($join_child_nodes_test)$tax_for_parent";
#	}



#print "tax_for_parent:$tax_for_parent\n";
#die;
$all_internal_taxa_on_constraint_tree{$tax_for_parent}=1;


if($next =~ /INTERNAL_NODE_\d/ && $all_child_nodes_are_internal == 1)
	{
#	$all_terminals_from_this_node ="";
#	get_terminals_from_this_node2($next, 0);
#	$tax_for_parent = get_shared_taxonomic_name($all_terminals_from_this_node);


	my $all_branches_from_this_node_are_the_same_taxon = 1;
	foreach my $test_taxon(@taxon_for_child_nodes)
		{
		if($tax_for_parent eq $test_taxon)
			{}else{$all_branches_from_this_node_are_the_same_taxon=0}
		};


	# collapse node method 1
	# turns not really what i want to do here.
	# but i keep the code in case it is what i want to do in a future situation
	if($all_branches_from_this_node_are_the_same_taxon ==1)
		{
		print "node $next will be COLLAPSED. tax_for_parent:$tax_for_parent taxon_for_child_nodes:@taxon_for_child_nodes\n";
	#	$swap_string = "$join_child_nodes";
		# collapse node. get identity of parent of current node (which is itself parent)
		# then add the child nodes to this  
		# no. instead you need to get the grandchild node of the matching taxon branch, and pull it up to be a child node.

		#$child_counts{$nodeID} = $#child_nodes;
		#$nodes{$child_nodes[$i]}{parent} 	= $nodeID;

		}


	};

my $join_the_child_nodes;
$collapse_nodes = 0;
if($collapse_nodes == 1)
	{
	# build the newick string, replace current internal node with child nodes
	print "node:$next becomes ($swap_string)\n";
	$new_newick_string =~ s/$next(\W)/$swap_string$1/;
	$newick_print_string=~ s/$next(\W)/$swap_string_test$1/;
	}else{
	@next1 = ();
	for $i2(0 .. $child_nodes)
		{push @next1, $nodes{$next}{$i2}};
	# default:
	$join_the_child_nodes = join ',', @next1;
	$swap_string = "($join_the_child_nodes)$tax_for_parent";
	$new_newick_string =~ s/$next(\W)/$swap_string$1/;
	$newick_print_string=~ s/$next(\W)/$swap_string$1/;

	};

# collapsing is acheived by not looping to child nodes, but skipping strait to grandchild nodes.
# similarly, grandchild nodes are printed to the growing newick string instead of child nodes.


for my $index(0 .. $#next1)
	{
############	my $test = @next1[$index];#print "test:$test\n";####################################### surely this is wrong?
	my $test = $next1[$index];

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node($test , $sum_branchlength+$branchlength_to_child );#	# recurse
		#####################################
		}
	}

return();

}






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################






sub get_shared_taxa
{
my $tax_list = shift;
my @tax_array = split /\t/ , $tax_list;
my $shared_tax_substring = "";
print "$#tax_array taxa input to sub get_shared_taxa\n";

# gets lineage (series of taxa) shared by all members. also option to get only lowest shared taxa

# first count how many have complete lineag available:

my $count_members_with_lineage =0;
my $single_lineage;my $first_lineage;my $found_first_lineage = 0;my $index_of_first_lineage;

foreach my $tax_index(0 .. $#tax_array) 
	{
	my $tax = $tax_array[$tax_index];
	my $test_taxonomy = $complete_lineage_of_terminal{$tax};
	#print "line928. test_species:$test_species test_taxonomy:$test_taxonomy\n";
	if($test_taxonomy =~ /\w+/)
		{
		$single_lineage = $test_taxonomy;$count_members_with_lineage++;
		if($found_first_lineage == 0)
			{$first_lineage = $test_taxonomy; $index_of_first_lineage = $tax_index;$found_first_lineage = 1};
		};
	};

print "count_members_with_lineage:$count_members_with_lineage\n";

my $most_inclusive_name = "NA";
my $most_inclusive_lineage = "";


 if($count_members_with_lineage == 1)
	{
#	die "\nyes 1424\n"
	while($single_lineage =~ s/^([^:]+):(\w+)//)
		{
		my $current_rank = $1;my $current_taxname = $2;
		$most_inclusive_name = $current_taxname;
		$most_inclusive_lineage .= "$current_taxname\t";
		};
	print "for just one, retrun lineage $most_inclusive_lineage\n";
	}elsif($count_members_with_lineage >= 2)
	{
	############################################	
	
	while($first_lineage =~ s/^([^:]+):(\w+)//)
		{
		my $current_rank = $1;my $current_taxname = $2;
		# over counts, should not count taxa at parent node:
		#$number_of_nodes_to_which_this_taxon_has_been_assigned{$current_taxname}++;
		#print "\tcurrent_rank:$current_rank current_taxname:$current_taxname\n";


		my $all_members_have_this_name =1;
		foreach my $k1( ( $index_of_first_lineage+1) .. $#tax_array) 
			{
			my $tax = $tax_array[$k1];
			my $test_taxonomy2 = $complete_lineage_of_terminal{$tax};
			if($test_taxonomy2 =~ /\w/)
				{
				#print "\tspecies:$tax complete taxonomy:$test_taxonomy2\n";
				if($test_taxonomy2 =~ /\:$current_taxname\s/)
					{print ""}else{
					$all_members_have_this_name =0;#print "0";
					}
				};
			}

		#print "\nall_members_have_this_name:$all_members_have_this_name\n";

		if($all_members_have_this_name == 1)
			{
			$most_inclusive_name = $current_taxname;
			$most_inclusive_lineage .= "$current_taxname\t";
			};


		}

	print "for several, retrun lineage $most_inclusive_lineage\n";

	############################################
	};



# print name only:
#return($most_inclusive_name);
# or lineage:
print "\treturning:($most_inclusive_lineage)\n";



return($most_inclusive_lineage, $count_members_with_lineage);

}





#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub get_species_list
{
my $list_of_terminals = shift;

my %hash_of_terminals = ();

if($verbose == 1){print "list_of_terminals:$list_of_terminals\n"};

my @array_of_terminals = split /\t/, $list_of_terminals;

foreach my $termnal(@array_of_terminals)
	{
	$termnal =~ s/^([^_]+)_.+/$1/;
	$hash_of_terminals{$termnal}= 1;
	}

my @array_of_terminals2 = keys %hash_of_terminals;

if($verbose == 1){
#print "list_of_species:@array_of_terminals2\n";
}

my $list_of_species = join ' ', @array_of_terminals2;

return($list_of_species);

}



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
     ***** _constraints.pl *****
		  
\n";


# perl parse_clades.pl -treefile RAxML_bipartitions.diet_group.0.reformatted -seqfile diet_group.0.fas -keyfile key_Mar2013_Magnoliophyta
#$treefile;
#$fas_file;
#$output_filename 	= "$treefile.diet_clades";



if($arguments =~ /-treefile\s+(\S+)/)
	{
	$treefile = $1;
	}else{
	print "\nerror reading command arguments  (-treefile)\n";die""
	}

if($arguments =~ /-support\s+(\S+)/)
	{
	$support_cutoff = $1;
	}else{
	$support_cutoff = 50;
	print "user did not give cutoff for acceptable boot support (not relevent if you are not using bootstrapped trees)\nusing default ($support_cutoff).\n"
	}

if($arguments =~ /-taxon_table\s+(\S+)/)
	{
	$taxon_table = $1;
	}else{
	print "\nerror reading command arguments  (-taxon_table)\n\n";die""
	}




if($arguments =~ /-seqfile\s+(\S+)/)
	{
	$fas_file = $1;
	}else{
	}

if($arguments =~ /-references\s+(\S+)/)
	{
	$reference_file = $1;
	}else{
	}

if($arguments =~ /-outgroup\s+(\S+)/)
	{
	$outgroup = $1;
	unless($outgroup =~ /[A-Z][a-z]+/){die "\noutgroup switch used but cannot parse name \n"}
	}else{
	print "\nyou did not specify the outgroup.\n";
	}

if($arguments =~ /-ignore_IDs_with_no_tax_assignment/)
	{
	$ignore_IDs_with_no_tax_assignment = 1;
	};

if($arguments =~ /-print_taxon_constraints_only/)
	{
	$print_taxon_constraints_only = 1;
	};

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


#$output_filename	= "$treefile.query_clades";


print "\n
user options have been read.....
support_cutoff:			$support_cutoff
taxon node encompassing refs:	$starting_node
treefile:			$treefile
query fasta file:		$fas_file
reference_file:			$reference_file

output will be written to file:	$output_filename


";



}#sub read_command_arguments



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



sub read_keyfile
{
open(IN, $keyfile) || die "\n\nerror 91. cant open $keyfile\n";
my $species_strings_read=0;

while (my $line = <IN>)
	{
	#print "\n$line";

	#MeeRBBegrac Berberis_gracilis 258166 species no_rank:eudicotyledons no_rank:eudicotyledons order:Ranunculales family:Berberidaceae genus:Berberis species:gracilis

	if($line =~ /^(\S+)\s(\S+)\s(\S+)\s(\S+)\s(.+)/)
		{
		my $tobycode = $1;my $species_name = $2;my $ncbi_number = $3; my $rank = $4; my $taxonomic_path = $5;
		#print "\ttobycode:$tobycode species_name:$species_name ncbi_number:$ncbi_number rank:$rank taxonomic_path:$taxonomic_path\n";
		$species_names{$tobycode} = $species_name;$taxonomic_paths{$tobycode} = $taxonomic_path;
		$species_strings_read++;
		$tobycodes{$species_name} = $tobycode;

	#	if($tobycode eq "MLAA4"){print "$tobycode $species_name ... quit\n\n";die ""}

		}
	}

close(IN);

if($species_strings_read == 0){die "\n\nerror 112. seems to be problem reading keyfile, lines in that file dont match regex.\n"}

print "$species_strings_read species strings read from file, plus associated taxonomic information\n";

}




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


# _THIS_SUB_NOT_RUN


sub traverse_nodes_THIS_SUB_NOT_RUN
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
# _THIS_SUB_NOT_RUN
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
# _THIS_SUB_NOT_RUN
	
	$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;

	#print "storing complete lineage for taxon:$name_string\n";	
	$complete_lineage_for_this_species{$name_string}=$child_complete_lineage;

	$ncbi_tax_number_for_this_species{$name_string}=$child;

#	print "complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";

	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;

	my $name_assignment_to_taxnumber = "";
# _THIS_SUB_NOT_RUN

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
# _THIS_SUB_NOT_RUN

	#print "$name_assignment_to_taxnumber $child $rank $child_complete_lineage\n";


		###########################################
		traverse_nodes_THIS_SUB_NOT_RUN($child , $child_taxstring);#
		###########################################
	}}

# _THIS_SUB_NOT_RUN
	
}#sub traverse_nodes





#####################################################################################################
#
#
#
#####################################################################################################




sub get_terminals_from_this_node2
{
my $next = shift;

my $child_nodes = $child_counts{$next};
my @next1 = ();
for $i(0 .. $child_nodes)
	{
	push @next1 , $nodes{$next}{$i};	
	#print "line1423node:$next i:$i child:$nodes{$next}{$i}\n";
	}



for my $index(0 .. $#next1)
	{
	my $test = @next1[$index];#print "test 1430:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{

		#####################################
		get_terminals_from_this_node2($test );#	# recurse
		#####################################

		}else{
		$count_terminal_nodes++;

		if($test eq "")
			{
			print "\nnext:($next)\n";
			die "\nerror 1443. looks like sub get_terminals_from_this_node2 was called with a terminal node.\n"
			}


		### here set the format used in the previous (mtgenome) tree.
		# currently im using Order_Genus in the tree, and parsing the genus name.
#		unless($test =~ s/.+_//)
#			{die "strange format:($test)\n"}

		# 8mar, Coleoptera_Tetraphalerus_bruchi and parsing genus
	#	unless($test =~ s/.+_(.+_.+)/$1/)
	#		{die "error 1229. expecting Order_Genus_species, got strange format:($test)\n"}
		# now process ids in intial reading of the tree, so they just have a regular taxon name

		$all_terminals_from_this_node .= "$test\t";

		#print "line1459. test:$test\t";
		}
	}

return();

}









#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################



# calls get_terminals_from_this_node2 (basic function), get_shared_taxa, get_all_taxa, test_for_monophylies


# 
# sub traverse_backbone_tree_and_infer_taxonomic_node_labels does 
#	taxon monophyly test applied to fasttree only, is relatively simple.
# pseudocode
#
# for each terminal of backbone tree, 
#	retrive all higher taxa and sum each [x]
# 	[count_taxa_represented_in_terminals]
# 
# recursing for each node:
#	get list of terminals derived from the node, retrive taxa shared by all members 
# 	for each higher taxa observed for the node:
# 		test criteria i) they are present for all terminals, ii) their count == x (i.e. not observed elswhere); 
#		for taxa meeting both criteria, stored as putative monophyly.
# 
# foreach barcode member
#	foreach monophyletic taxon
#		if member has no lineage information, code -, else,
# 		if member belongs to taxon, code 1,
# 		else code 0.
#
#
# sub traverse_backbone_tree_and_infer_taxonomic_node_labels also seems to do independent older version of taxonomic inference 
#
#



sub traverse_backbone_tree_and_infer_taxonomic_node_labels
{

my $next = $_[0];my $sum_branchlength = $_[1];


print "\nNEW NODE ($next)\n";

my @non_terminal_child_nodes = ();
my @child_nodes_test = ();
my $all_terminals_from_both = "";
my @taxon_for_child_nodes = ();

my $lowest_taxname_for_node = "";
my $lineage_for_node = "";my $count_of_terminals_with_lineage;
$all_terminals_from_this_node ="";
$taxa_for_terminals_of_node = "";




#################################################################################################################################



	# monophyly test below, applied to fasttree table


# need a list of taxa derived from currnet node. ($all_terminals_from_this_node .= "$test\t";)
# get tax names for these.
# then go through all remaining members of the tree, (%terminals, not incl $all_terminals_from_this_node)
# if none of these have these tax names (get_shared_taxa), they are defined as a constraint.
# then all new members from these, are constrained.			


		# do this even if child node is terminal
if($next =~ /INTERNAL_NODE_\d/) 	# && $all_child_nodes_are_internal == 1)
	{

	  # following is very simple sub, recurses backbone tree from current node and stores all terminals in a list.
	##########################################
	get_terminals_from_this_node2($next, 0);#
	##########################################



	# following sub takes as input the list of terminals derived from node. 
	# and gets lineage (series of taxa) shared by all members of clade. 
	# also option to get only lowest shared taxa

				####################################################
	my @sub_output = 	get_shared_taxa($all_terminals_from_this_node);#
				####################################################
	$lineage_for_node = $sub_output[0];
	$count_of_terminals_with_lineage = $sub_output[1];

	if($lineage_for_node =~ /\w/)
		{
	#print "get_all_taxa 1\n";

	print BACKBONE_LOG "\n$next is internal node, lineage_for_node:$lineage_for_node\n";


	# returns a complete list of all the taxa found in the lineages of the terminals derived from current node.
	#  (might not neccessrily find anything)

					###################################################
	$taxa_for_terminals_of_node = get_all_taxa($all_terminals_from_this_node);#
					###################################################
	unless($taxa_for_terminals_of_node =~ /\w/){die "error 1950 could not parse taxa for terminals of node ($next)\n"};


	if($all_terminals_from_this_node =~ /\w.*\t+.*\w/ && $taxa_for_terminals_of_node =~ /\w/ && 
		$count_of_terminals_with_lineage >= 2)
		{

	#########################################################
	test_for_monophylies($all_terminals_from_this_node);#
	#########################################################
		};

		}else{
		print  "\nERROR 2190. Could not parse lineage for node ($next)\n";
	#	die "";
		};



	}else{#if($next =~ /INTERNAL_NODE_\d/)

	# for tips, you still need to know the lineage. but it wont be a shared lineage.
	# in this case you dont have a list of terminals, theres just one, the curent node

				#########################
	my @sub_output = 	get_shared_taxa($next);#
				#########################

	$lineage_for_node = $sub_output[0];
	$count_of_terminals_with_lineage = $sub_output[1];

	unless($lineage_for_node =~ /\w/)
		{
		print "error 1963 could not parse lineage for node ($next)\n";
	#	die "";
		};

	#print "get_all_taxa 2\n";
					###################################################
	$taxa_for_terminals_of_node = get_all_taxa($next);# bug fix 22 oct 2015
					###################################################
	unless($taxa_for_terminals_of_node =~ /\w/)
		{
		print "error 1969 .could not parse taxa for terminals of node ($next)\n";
		# die "";
		};

	##########################################################################################
	};# if non tip, then if tip

######################################################################################################################################



# parse all names, put into this object
my $lineage_string = $lineage_for_node; # series of taxa which are shared by all members of clade
while( $lineage_string =~ s/([^\t]+)\t$// )
	{
	my $taxon_for_node = $1;$all_taxa_on_backbone_tree{$taxon_for_node}=1; # shared taxa from current node placed into hash.
	#print "taxname1:$taxon_for_node lineage_string:$lineage_string\n";
	};


# for current node, parse lowest (tip level) shared taxa from the lineage:
$lowest_taxname_for_node = $lineage_for_node; # print "lineage_for_node:$lineage_for_node\n"; # lineage_for_node:Hymenoptera	Octostruma

unless($lowest_taxname_for_node =~ s/.+\t([^\t]+)\t$/$1/)
	{
	unless($lowest_taxname_for_node =~ s/^$starting_name_ff	/$starting_name_ff/ )# basal TOL will have shared taxon being cellular orgs, same as root specified by user
		{
		print "next:$next\nlineage_for_node:$lineage_for_node\ntaxa_for_terminals_of_node:$taxa_for_terminals_of_node\nerror 1982\n";
		print "\nERROR 1550, lowest_taxname_for_node:($lowest_taxname_for_node)\n";
		# die "";
		};
	};
# print "lowest_taxname_for_node:$lowest_taxname_for_node\n";# lowest_taxname_for_node:Octostruma

# check for error:
unless($lowest_taxname_for_node =~ /^[A-Za-z\_]+$/)
	{
	print "\nERROR 1638. lowest shared taxa for current node ($next) has unexpected format:($lowest_taxname_for_node).\n";
#	die "";
	}
# and store the LOWEST name in this hash (no intermediates)
$node_label{$next} = $lowest_taxname_for_node;

print "\tlowest_taxname_for_node:$lowest_taxname_for_node\n";

print BACKBONE_LOG "\'lowest_taxname\' for_node:$lowest_taxname_for_node\n";

$all_taxnames_for_terminals_decended_from_this_node{$next} = $taxa_for_terminals_of_node;







#######################################################################################################################################

# following code retreives taxa intermediate of node and parent......


# for each node, 
#	retreive taxa shared by all decendents, including lowest ranking taxon  
# 
# next compare the complete shared lineage of the current node, with that of the parent.
# in the example below, the node in which all decendents are Lucanidae,
# in addition to being the position at which new Lucanids are added
# can also accomodate the above ranks Scarabaeoidea, Polyphaga,
# but not Coleoptera, which is best assigned to the parent node.
# therefor find (if present) all intermediate ranks.
#        
#                    -------
#                    |
#           -----Lucanidae
#           |        |
#           |        ------
# ---Coleoptera      
#           |        -----
#           |        |
#           -----Adephaga
#                    |
#                    -----


my $parent_node = $nodes{$next}{parent};
my $parent_taxon = "";
my $all_potential_taxa ="";


if($node_label{$parent_node} =~ /\w/) # probably stored in previous loop. contains single taxon.
	{

	$parent_taxon = $node_label{$parent_node}; # print "node($next), tax($lowest_taxname_for_node) parent($parent_node) Tax($parent_taxon)\n";


	if($lowest_taxname_for_node eq $parent_taxon)	# compare the lowest tax name of current node
		{						# to lowest tax name of parent node

		# Node and parent both have the same assigned taxon, nothing to do.
		# Although whilst here, need to count the number of nodes each taxon has been assigned to 
		$number_of_nodes_to_which_this_taxon_has_been_assigned{$lowest_taxname_for_node}++;

		}else{
		while($lineage_for_node =~ s/([^\t]+)\t$//)
			{
			# go through whole lineage of current node 
			# (including lowest, which is not trimmed from this particular variable), 
			# and stop when reaching the taxon assigned to parent
			# this direction:Corcyra_cephalonica->Corcyra->Galleriinae->Pyralidae

			my $intermediate_taxon_for_node = $1;
			unless($intermediate_taxon_for_node =~ /^[A-Za-z\_]+$/){die "\nerror 1691. unexpected taxon name:$intermediate_taxon_for_node\n"}
			#print "\ttaxname1:$intermediate_taxon_for_node. comparing to that assigned to parent ($parent_taxon)\n";

			# the first element is the lowest taxon of current node
			# if that is the same as that of the parent node, 
			# there will be nothing put in $all_potential_taxa
			# if there is only one item in $all_potential_taxa, it will be same as lowest taxon
			if($intermediate_taxon_for_node eq $parent_taxon)
				{
				#print "same, break.\n";
				$lineage_for_node = "";# delete lineage string to break loop.
				$all_potential_taxa =~ s/[\s\.]+$//;# if names have been assigned, there will be terminal char
				}else{

				# taxon assigned to parent node has not been reached yet, so append taxon name 
				$all_potential_taxa .= "$intermediate_taxon_for_node.";#print "appending\n";

				# And whilst here, keep counting the number of nodes each taxon has been assigned to 
				$number_of_nodes_to_which_this_taxon_has_been_assigned{$intermediate_taxon_for_node}++;
				};
			};
		}

	}else{#if($node_label{$parent_node} =~ /\w/)
	
	unless($next eq $root_node)# only the root node wont have a parent, so no taxon will be found for that one.
		{
		print "\nerror 1703. no taxon assigned to parent ($parent_node) of node ($next), parent_node:$parent_node\n";
	#	die "";
		};

	};


#######################################################################################################################################



	# one is chosen, and it is pplaced to the currnet node in $new_newick_string.



my $label_string_to_assign = "";
if($all_potential_taxa =~ /\w[\.\s]\w/)# taxa found which are intermediate between those assigned to node and parent
	{
	$label_string_to_assign = $all_potential_taxa	# $all_potential_taxa has fullstop seperated taxon names
	}else{							# so they can be printed by tree viewers.

	$label_string_to_assign = $lowest_taxname_for_node # otherwise, just use lowest shared taxon name
	};
unless($label_string_to_assign =~ /\w\w+/)
	{
	print "\nERROR 2344. why no taxon assigned to current node ($next). \n";
#	die "";
	}


my $child_nodes = $child_counts{$next};
my @next1 = ();

# record the taxonomic name(s) given to this node
$nodes{$next}{node_taxonomic_label}=$label_string_to_assign;
print BACKBONE_LOG "label_string_to_assign:$label_string_to_assign\n\n";


#print "\n\nNEW BACKBONE NODE ($next) assigned name(s):$label_string_to_assign\n";







for $i(0 .. $child_nodes)
	{
	# specify new array of the child nodes as normal for recursing
	my $push_node = $nodes{$next}{$i};
	push @next1 , $push_node;	
	}

#print "node:($next) child nodes:(@next1) label_string_to_assign:($label_string_to_assign) lineage_for_node:($lineage_for_node)\n";

#if($next eq "Apis_florea"){die""};

# get to the tip,
#NEW NODE (Apis_florea)
#node:(Apis_florea) child nodes:() label_string_to_assign:()

# not used any more ?:
# build newick string. this string is an intermediate, 
# it contains empty parentheses ready to be filled with the species-dense data,
#thus not readable by tree viewers

my $swap_string 	= "";
$join_the_child_nodes 	= join ',', @next1;
$swap_string 		= "($join_the_child_nodes)$label_string_to_assign";
$new_newick_string 	=~ s/$next(\W)/$swap_string$1/;

# alternative newick string that is readable at this stage:

if(exists($child_counts{$next}))
	{
	$new_newick_string2 	=~ s/$next(\W)/$swap_string$1/;
	}else{
	# instead of replacing node (which is now a species name) with child nodes (which there arnt any)
	# replace node with lineage
	$new_newick_string2 	=~ s/$next(\W)/$label_string_to_assign$1/;
	}


$all_internal_taxa_on_constraint_tree{$tax_for_parent}	= 1;

# the approach in this sub, is called all the way up to tip nodes, so tax lineages etc can be easier assessed.
# which means it needs to determine the end has been reached below (no child nodes for this tip)
# its a bit different from the usual approach, in which the sub is not called for tip nodes.

if(exists($child_counts{$next}))
	{

	for my $index(0 .. $#next1)
		{
		my $test = $next1[$index];

		#####################################
		traverse_backbone_tree_and_infer_taxonomic_node_labels($test , $sum_branchlength+$branchlength_to_child );#	# recurse
		#####################################

		}

	};

return();

}; # sub traverse_backbone_tree_and_infer_taxonomic_node_labels






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


# sub currently NOT IN USE


sub get_terminals_from_this_node30
{
my $next = $_[0];
my $child_nodes = $child_counts{$next};

#print "\nNEW NODE:($next)\n";

my @next1 = ();


my $add_to_current_node = "";
my $check1=0;
my $backbone_child_nodes_absent_from_new_IDs=0;
my $count_new_IDs_added_to_current_node = 0;

for $i(0 .. $child_nodes)
	{


	# get child nodes as normal
	my $push_node = $nodes{$next}{$i};

	unless(exists($query_IDs{$push_node}))
		{
		unless($push_node =~ /INTERNAL_NODE_/)
			{
			$backbone_child_nodes_absent_from_new_IDs++;
			print "warning 1875. node $push_node is not internal, yet is not found in list of new IDs\n";
			};
		}

	push @next1 , $push_node;	
	#print "child node:$i node ID:$push_node\n";

	my $node_taxonomic_label = $nodes{$push_node}{node_taxonomic_label};

	my @split_node_tax_labels = split /\./ , $node_taxonomic_label;
	foreach my $candidate_tax(@split_node_tax_labels)
		{
		#print "\tchild node is assigned candidate_tax:$candidate_tax\n";# Endopterygota->Acyrthosiphon_pisum->Acyrthosiphon

		unless($candidate_tax =~ /^[A-Za-z\_]+$/){die "\nerror 1864. unexpected taxon name:$candidate_tax\n"}

		if(exists($assign_these_new_terminals_to_taxon{$candidate_tax}))
			{
			my $string_length_new = length($assign_these_new_terminals_to_taxon{$candidate_tax});

			#print "child $i of node $next has tax label:$candidate_tax\n";
			#print "\t\tthis taxon has new IDs to be added, length:$string_length_new\n";

			# these are new IDs, comma seperated. will have comma terminating string.
			my $new_members = $assign_these_new_terminals_to_taxon{$candidate_tax};# . ",";
			$new_members =~ s/\,$//;

			# each new member come across, count in hash.
			# later check duplicates are not met.
			my @split_new_members= split /\,/ , $new_members;
			foreach my $newID(@split_new_members)
				{
				$store_node_assigned_newIDs{$newID}++;$count_new_IDs_added_to_current_node++;
				unless($newID eq $push_node){$add_to_current_node .= $newID . ",";}
				};

			# append, this will be a new child of the current node
			#$add_to_current_node .= $new_members . ",";
			my @check_split= split /\,/ , $add_to_current_node;

			#if($add_to_current_node =~ /pisumAcyr/){die "\n1887.\n"}

			#$check1++;

			}#if(exists($assign_these_new_terminals_to_taxon{$candidate_tax}))

		};#foreach my $candidate_tax(@split_node_tax_labels)

	};

#print "
#count_new_IDs_added_to_current_node:$count_new_IDs_added_to_current_node
#backbone_child_nodes_absent_from_new_IDs:$backbone_child_nodes_absent_from_new_IDs
#";
# check any havnt been assigned elsewhere (shouldnt have)
my @all_new_IDs_assigned_at_this_point = keys %store_node_assigned_newIDs; 
foreach my $newID(@all_new_IDs_assigned_at_this_point)
	{
	if($store_node_assigned_newIDs{$newID}>1)
		{
		print "\nerror 1880. new ID:($newID), to be assigned to node:($next),\n";
		print "which has tax:(...), was assigned previously. quit.\n";die};
	#unless($newID eq ""){};# terminal comma needs keeping
				# so last object will be empty
	};




my $copy = $add_to_current_node;#$copy .= join ',', @next1;
#$copy =~ s/\,$//;
my @check_split= split /\,/ , $copy;
my %check_hash=();
foreach my $check4(@check_split)
	{
#unless($check4 =~ /^[A-Z][a-z]+_[a-z]+$/){die "\nerror 1828:$check4\n"}
	if(exists($query_IDs{$check4}))
		{
	$check_hash{$check4}++;#if($check_hash{$check4}>= 2){die "\nerror 1822, $check4 assigned already\n"}
	#$check_hash4{$check4}++;#if($check_hash4{$check4}>= 2){die "\nerror 1824, $check4 assigned already\n"}
		}
	}

my @check_split2= keys %check_hash , @next1;
my %check_hash2=();
foreach my $check4(@check_split2)
	{

	$check_hash2{$check4}++;#if($check_hash4{$check4}>= 2){die "\nerror 1824, $check4 assigned already\n"}

	}


my @new_nextnode = keys %check_hash2;

#if($check1>= 3){die "\nerror 1816\n"}


#print "node:$next child nodes:@next1 node_taxonomic_label:$node_taxonomic_label\n";
if($add_to_current_node =~ /./)
	{
	#print "new members will be added to this node, ($#check_split)\n"
	};
#if($next1[1] =~ /Trachypachus_holmbergi/i){die}

#	


$howmany+= $backbone_child_nodes_absent_from_new_IDs;

if(
$backbone_child_nodes_absent_from_new_IDs >= 1 && 
$count_new_IDs_added_to_current_node < 2
){
	# this adds to 2, so ok
unless($count_new_IDs_added_to_current_node==1 && $backbone_child_nodes_absent_from_new_IDs == 1)
	{
unless($count_new_IDs_added_to_current_node == 0 && $backbone_child_nodes_absent_from_new_IDs == 1)
{
#$howmany++;
#die "\nyou have a problem. current node ... \n";
}	}
}

# default:
my $join_child_nodes = $add_to_current_node . join ',', @next1;

$swap_string = "($join_child_nodes)";
unless($newick5 =~ s/$next(\W)/$swap_string$1/)
	{print "\nerror 1962, cant build string\n"};

for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];#print "test:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node30($test  );#	# recurse
		#####################################
		}
	}

return();

}








#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################




 #   sub currently NOT IN USE




sub get_terminals_from_this_node31
{

my $next 	= $_[0];
my $child_nodes = $child_counts{$next};


my @next1 = ();
my $add_to_current_node = "";
my $check1=0;
my $backbone_child_nodes_absent_from_new_IDs=0;
my $count_new_IDs_added_to_current_node = 0;

my $node_taxonomic_label = $nodes{$next}{node_taxonomic_label};

# keep this print: 
print "\n*** NEW NODE of backbone *** 
Making TNT / FT constraint no:$constraint_number; node:($next); tax:($node_taxonomic_label)\n";


 # sometimes they have multiple names a la: Diabrotica.Diabroticites.Diabroticina.Luperini.Galerucinae
my @split_node_tax_labels = split /\./ , $node_taxonomic_label;

my $new_IDs_going_to_this_node = "";
my $count_tax_assigned_to_node =0;
my $total_tax_assigned_to_node  =scalar @split_node_tax_labels;
foreach my $candidate_tax(@split_node_tax_labels)
	{
	$count_tax_assigned_to_node++;
	print "\ttax:$count_tax_assigned_to_node of total:$total_tax_assigned_to_node assigned. taxon name:$candidate_tax\n";# Endopterygota->Acyrthosiphon_pisum->Acyrthosiphon

	unless($candidate_tax =~ /^[A-Za-z\_]+$/){die "\nerror 1864. unexpected taxon name:$candidate_tax\n"}

	if(exists($assign_these_new_terminals_to_taxon_TNT{$candidate_tax}))
		{
		print "\tIS new members to assign to this name\n";

		if(exists($assigned_already{$candidate_tax}))
			{
			print "\twarning, taxon($candidate_tax) was come across on previous node. ignoring.\n";
			}else{
			#print "child $i of node $next has tax label:$candidate_tax\n";
			#print "\t\tthis taxon has new IDs to be added\n";

			# these are new IDs, comma seperated. will have comma terminating string.
			my $new_members = $assign_these_new_terminals_to_taxon{$candidate_tax};# . ",";
			$new_IDs_going_to_this_node .= $new_members;
			$new_members =~ s/\,$//;
			my @new_members_to_node = split /\,/ , $new_members;
			print "\tcount new members to this node:" , scalar @new_members_to_node , "\n";
			};

		$assigned_already{$candidate_tax} = 1;
		}else{ 
		# if(exists($assign_these_new_terminals_to_taxon{$candidate_tax}))
		print "\tNO new members to assign to this name\n";
		}

	}; # foreach my $candidate_tax(@split_node_tax_labels)






# constrinat is made of the new taxa, and the list of terminasl in backbone tree, so get the terminals:
$all_terminals_from_this_node = "";
# does this:$all_terminals_from_this_node .= "$test\t";

#######################################
get_terminals_from_this_node2($next);#
#######################################





unless($next eq $root_node)
	{



# 3 components to each constraint: 1)the taxa monophyletic 2) the taxa excluded 3) the taxa that float.
# err on the side of floating!

my %floaters = ();
my %constrained = ();

# all backbone IDs:		%terminals
# backbone IDs to constrain:	$all_terminals_from_this_node


my @all_backbone_IDs = keys %terminals;@all_backbone_IDs = sort @all_backbone_IDs;
foreach my $backbone_ID(@all_backbone_IDs)
	{
	#print "backbone_ID:$backbone_ID\n";	
	if($all_terminals_from_this_node =~ /$backbone_ID/)
		{
		#print "\tconstrain MP\n";
		$constrained{$backbone_ID}=1;# taxa in the backbone tree decended from this node will be constrained MP
		}else{
		$constrained{$backbone_ID}=2;# taxa in backbone tree not decended from this node will be excluded from MP
		#print "\texclude from MP\n";#die;
		};
	};

#$all_terminals_from_this_node

my @all_new_IDs = keys %query_IDs;
foreach my $new_ID(@all_new_IDs)
	{
	my $check_overlapping = $constrained{$new_ID};
	#print "new_ID:$new_ID check_overlapping:$check_overlapping\n";
	unless($check_overlapping =~ /\d/)
	{
	if($new_IDs_going_to_this_node =~ /$new_ID/)
		{
		#print "\tconstrain MP\n";
		$constrained{$new_ID}=1;# new taxa assigned to this node will be constrained MP
		}else{
		#print "\tfloat\n";
		$constrained{$new_ID}=3;# new taxa not assigned to this node will float
		};
	}};



my @everyone = keys %constrained;@everyone = sort @everyone;

my $constrain_string="";
my $float_string="";
my %check_assigned = ();
$constraint_number++;


foreach my $member(@everyone)
	{
	#print "member:$member\n";
	my $contr = $constrained{$member};
	my $fasta_index;
	if($index_of_fastaID{$member} =~ /\d/)
		{
		$fasta_index = $index_of_fastaID{$member};
		if(exists($check_assigned{$fasta_index}))
			{die "\nerror 2218. member:$member fasta_index:$fasta_index\n"};
		$check_assigned{$fasta_index}=1;
		}else{die "\nerror 2190. constr:$contr\tmem:$member\t\n"}
	#print "$contr\t$member\tfile index:$fasta_index\n";

	if($contr == 1){$phylip_constraints{$member} .= "1";$constrain_string .= "$fasta_index "}
	if($contr == 2){$phylip_constraints{$member} .= "0"};
	if($contr == 3){$phylip_constraints{$member} .= "-";$float_string .= "$fasta_index "}
	}

$constrain_string	=~ s/\s$//;
$float_string		=~ s/\s$//;

my $tnt_command;
$tnt_command = "[ $constrain_string ($float_string) ]";

$tnt_command_string .= "force + $tnt_command, ";
#print TNT_FILE_HANDLE "$tnt_command_string";


#		$index_of_fastaID{$id} = $count_diet_IDs;
# tnt format, 
#	# constraint where 2 taxa 'float'
#	force + [ a b c  (d e) ]
#tnt64 mxram 2000, rseed 1,p default-Second.nex timeout 0:20:00, echo=, taxname=, force + [17 18 19 20], constrain =,
#
#$tnt_command_string
	# crashes unless you have single quotes around force command, due to parentheses confusing bash.
#	tnt64 mxram 2000, rseed 1,p tnt_input.nex timeout 0:00:30, echo=, taxname=, 'force + [ 14 46 67 (1 2) ],' constrain =, log tnt_logfile,  hold 2000, xmult=replication 100 drift 10 ratchet 10 hold 10, best, majority *, tchoose / , export *tnt_output, log / ,quit;
#print "backbone constraint (list of all_terminals_from_this_node):$all_terminals_from_this_node\n";
#print "new_IDs_going_to_this_node:$new_IDs_going_to_this_node\n";


# fasttree constraint file.
# 5 2
#A        00
#B        0-
#C        10
#D        11
#E        -1




	};#unless($next eq $root_node)




for $i(0 .. $child_nodes)# get child nodes as normal
	{
	push @next1 , $nodes{$next}{$i};
	#print "child node:$i node ID:$nodes{$next}{$i}\n";
	};

for my $index(0 .. $#next1)
	{
	my $test = $next1[$index];#print "test:$test\n";

	if($test =~ /^INTERNAL_NODE_/)
		{
		#####################################
		get_terminals_from_this_node31($test  );#	# recurse
		#####################################
		}
	}

return();

}


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





sub traverse_backbone_tree_and_collapse
{

my $next = shift;

# this has species names, species groups, genera, families etc. and has the outgroup on first call.
my $all_taxnames_for_terminals = $all_taxnames_for_terminals_decended_from_this_node{$next};
#print "all_taxnames_for_terminals:($all_taxnames_for_terminals)\n";


# taxonomic name assigned to node of backbone tree:
my $label = $node_label{$next};# just one name in this
# or this one:
#my $label = $nodes{$next}{node_taxonomic_label};
#print "\nNEW NODE ($next), label:$label\n";


my @non_terminal_child_nodes = ();my @child_nodes_test = ();
my $all_terminals_from_both = "";my @taxon_for_child_nodes = ();
my $lowest_taxname_for_node = "";my $lineage_for_node = "";
$all_terminals_from_this_node ="";


#my @tax_labels_current_node_of_backbone_tree = split /\./ , $label;



#foreach my $taxlab(@tax_labels_current_node_of_backbone_tree)
#	{
my $taxlab = $label;
# print "label assigned to current node of backbone:$taxlab\n";
my $tax_number = $ncbi_tax_number_for_this_species{$taxlab};

# all taxa at rank below the name:
my $child_nodes = $ncbi_nodes{$tax_number}{child_nodes};$child_nodes =~ s/^\t//;

# print "tax_number:$tax_number has child_nodes:$child_nodes\n";
my @child_nodes_array = split(/\t/, $child_nodes);
my $current_constraint = "";

my $mpA =0;my $mpB =0;

# each child node (taxonomic name) derived from name that was given to current node of tree
foreach my $child(@child_nodes_array)
	{
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g;	
#	print "child taxon ID:$child name:$name_string\n";

	# for each child taxon (1 rank below) derived from the 
	# taxonomic label assigned to the current node of the backbone tree.

	my $found_in_backbone_decendent=0;
	if($all_taxnames_for_terminals =~ /\t$name_string\t/)
		{
		# child taxon is found in the lineage of terminals decended from this backbone node			
		# so a monophyletic constraint may be made			
		$found_in_backbone_decendent =1;#print "\tIS in list of taxnames of decendents\n";

		#print "\t\t\ttax found in a terminal of this node of backbone tree\n"
		}else{
		#print "\tNOT in list of taxnames of decendents\n";
		#print "\t\t\ttax NOT found in a terminal of this node of backbone tree\n";
		};

	# how many times has this taxon been assigned, ANYWHERE on the backbone tree
	my $number_anywhere_in_backbone = $number_of_nodes_to_which_this_taxon_has_been_assigned{$name_string};

	# is this taxon found in the fasta file, if not, ignore.
	my $how_many_in_fasta = $count_tax_in_new_data{$name_string};

	# if child taxon if found in backbone tree and new fasta members, but is not in decendents of current node, 
	# will constrained outside of monophyly
	if($number_anywhere_in_backbone >= 1 && $how_many_in_fasta >= 1 && $found_in_backbone_decendent == 0)
		{$current_constraint .= "$name_string:0\n";$mpA=1;print "mpA:$mpA, constrain $child:$name_string outside MP\n"
		}

	# constrained MP
	if($number_anywhere_in_backbone >= 1 && $how_many_in_fasta >= 1 && $found_in_backbone_decendent == 1)
		{$current_constraint .= "$name_string:1\n";$mpB=1;print "mpB:$mpB, constrain $child:$name_string inside MP\n"
		}

	# no information on where these should be placed (not on backbone). so they float where they want.
	if($number_anywhere_in_backbone <1 && $how_many_in_fasta >= 1)
		{$current_constraint .= "$name_string:-\n"}

#	print "\t$child:$name_string how_many_in_fasta:$how_many_in_fasta\n";
#	print "\tnbr_anywhere_in_backbone:$number_anywhere_in_backbone found_in_backbone_decendent:$found_in_backbone_decendent\n";


	};#foreach my $child(@child_nodes_array)


my $constraint_length = length($current_constraint);

#print "constraint_length:$constraint_length mpA:$mpA mpB:$mpB\n";

# mpA and mpB check that there are taxa defined within and external , i.e. a monophyly to be made.
if($constraint_length < 4000 && $mpA ==1 && $mpB == 1)
	{
	print "Found relational constraint on backbone tree:\n($current_constraint)\n\n";
	$ft_constraints{$current_constraint}=1;

#:Paraneoptera how_many_in_fasta:1 number_in_backbone:1 found_in_backbone_terminal:1
#	33392:Endopterygota how_many_in_fasta:38807 number_in_backbone:7 found_in_backbone_terminal:1
	};





## too many cases of taxon with none in the constraint tree ....
## solutions, if very few, then let them float,
## or, add taxonomy-based constraint to the node, as a new child node, then put the new seqs in there! easy.
## alternative, make a backbone tree with more complete deep level sampling?


my $child_nodes = $child_counts{$next};
my @next1 = ();
for $i(0 .. $child_nodes)
	{
	# specify new array of the child nodes as normal for recursing
	my $push_node = $nodes{$next}{$i};
	push @next1 , $push_node;	
	}

my $swap_string 	= "";
$join_the_child_nodes 	= join ',', @next1;
$swap_string 		= "($join_the_child_nodes)$label_string_to_assign";
#$new_newick_string 	=~ s/$next(\W)/$swap_string$1/;

# alternative newick string that is readable at this stage:

if(exists($child_counts{$next}))
	{
	#$new_newick_string2 	=~ s/$next(\W)/$swap_string$1/;
	}else{
	# instead of replacing node (which is now a species name) with child nodes (which there arnt any)
	# replace node with lineage
	#$new_newick_string2 	=~ s/$next(\W)/$label_string_to_assign$1/;
	}


if(exists($child_counts{$next}))
	{
	for my $index(0 .. $#next1)
		{
		my $test = $next1[$index];

		#####################################
		traverse_backbone_tree_and_collapse($test );#	# recurse
		#####################################
		}
	};

return();

}






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################





sub get_all_taxa
{


my $tax_list = shift;
my @tax_array = split /\t/ , $tax_list;
my $all_tax_names_for_these_terminals = "";

unless($tax_list =~ /./){die "\nerror 2830, why has sub get_all_taxa been called with nothing\n"};

foreach my $tax(@tax_array)
	{
#	my $test_taxonomy = $complete_lineage_for_this_species{$tax};
	my $test_taxonomy = $complete_lineage_of_terminal{$tax};

#	my $genus_name = $tax;$genus_name =~ s/^([A-Z][a-z]+)_.+/$1/;
#	unless($test_taxonomy =~ /\w/){$test_taxonomy = $complete_lineage_for_this_species{$genus_name}};

	#print "\tspecies:$tax complete taxonomy:$test_taxonomy2\n";

	while($test_taxonomy =~ s/^([^:]+):(\w+)//)
		{
		my $current_rank = $1;my $current_taxname = $2;
		unless($all_tax_names_for_these_terminals =~ /\t$current_taxname\t/)
			{$all_tax_names_for_these_terminals .= "	$current_taxname	";
			}
		}
	}
#print "\ncount tax_array:" , scalar @tax_array , ", all_members_have_this_name:$all_members_have_this_name\n";


return($all_tax_names_for_these_terminals);


}






#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################


sub test_for_monophylies
{
my $tax_list = shift; # input to this sub is list of terminals derived from a given node.

unless($tax_list =~ /\w.*\t+.*\w/){die "\nwhy has sub test_for_monophylys been called on single item, not list. quitting ....\n"};
my @tax_array = split /\t/ , $tax_list;
my %all_tax_names_for_these_terminals = ();
my $how_many_terminals_have_lineage;

print BACKBONE_LOG "\tsub test_for_monophylies\n";


#count_taxa_represented_in_terminals


foreach my $tax(@tax_array) # for each terminal derived from certain node of the backbone tree.
	{
#	my $test_taxonomy = $complete_lineage_for_this_species{$tax};
	my $test_taxonomy = $complete_lineage_of_terminal{$tax};

	if($test_taxonomy =~ /\w/)
		{
		$how_many_terminals_have_lineage++;
	#print "\tspecies:$tax\n";# complete taxonomy:$test_taxonomy\n";

		my %tax_encounted_in_current_lineage = ();
		while($test_taxonomy =~ s/^([^:]+):(\w+)//)
			{
			my $current_rank = $1;my $current_taxname = $2;
			if($tax_encounted_in_current_lineage{$current_taxname} == 1)
				{
				print "tax name $current_taxname is duplicated in ncbi lineage\n";
				}else{
				$all_tax_names_for_these_terminals {$current_taxname}++;
			#	print "\t\tcurrent_rank:$current_taxname:" , $all_tax_names_for_these_terminals{$current_taxname}, "\n";
				};
			$tax_encounted_in_current_lineage{$current_taxname} = 1;
			};

		};

	};


#print "\nall_members_have_this_name:$all_members_have_this_name\n";


#count_taxa_represented_in_terminals

my @all_taxnames_from_terminals_of_current_node = keys %all_tax_names_for_these_terminals;
@all_taxnames_from_terminals_of_current_node = sort @all_taxnames_from_terminals_of_current_node;


foreach my $tax(@all_taxnames_from_terminals_of_current_node)
	{
	my $count_total_in_backbone_terminals = $count_taxa_represented_in_terminals{$tax};
	my $count_decended_from_current_node = $all_tax_names_for_these_terminals{$tax};
#	my $how_many_terminals_from_current_node = scalar @tax_array;
	my $recorded_as_monophyletic = 0;
	if($count_total_in_backbone_terminals == $count_decended_from_current_node && 
#		$count_total_in_backbone_terminals == $how_many_terminals_from_current_node
		$count_total_in_backbone_terminals == $how_many_terminals_have_lineage
		)
		{
		$monophylys{$tax}=1; $recorded_as_monophyletic = 1;# print "putative monophyly:$tax\n";
		}

	print BACKBONE_LOG "\t\ttax:$tax recorded_as_monophyletic:$recorded_as_monophyletic " , 
	"count_total_in_backbone_terminals:$count_total_in_backbone_terminals " , 
		"count_decended_from_current_node:$count_decended_from_current_node " , 
		"how_many_terminals_have_lineage:$how_many_terminals_have_lineage\n";
	#	"how_many_terminals_from_current_node:$how_many_terminals_from_current_node\n";

	}




return($all_tax_names_for_these_terminals);








};



#######################################################################################################################################
#
#
#
#
#
#######################################################################################################################################

#####################################################################################################
#
#
#
#####################################################################################################


sub store_tax_heirarchy
{
print "
reading users taxonomic file $taxon_table
";

open(TAXTABLE, $taxon_table) || die "\ncant find taxon table.\n";
my $tax_table_line_count=0;
while (my $line = <TAXTABLE>)
	{

# 	"$child\t" . 		# child node ID
#	"$current_node\t" . 	# parent node id
#	"$name_string\t" . 	# child tax
#	"$parentname\t" . 	# parent tax
#	"$rank\n"; 		# child rank

	$line =~ s/\n//;$line =~ s/\r//;
	if($line =~ /^(.+)\t(.+)\t(.+)\t(.+)\t(.+)/)
		{
		my $child_ID = $1;my $parent_ID = $2;my $child_tax = $3;my $parent_tax = $4;my $child_rank = $5;
		$rank_hash{$child_rank}++;
		if($tax_table_line_count == 0)
			{
			$root_taxon_name = $parent_tax;
			$starting_node = $parent_ID; 			# root node should be first parent node in the file
			$ncbi_nodes{$parent_ID}{name} = $parent_tax; 	# as names are assigned only for child nodes below, 
			}; 						# the root node needs name assigning by parent


#		my $tax_id = $1;my $parent_tax_id = $2;my $rank = $3;
#		print "tax_id:$tax_id parent_tax_id:$parent_tax_id rank:$rank\n";

		my $current_rankcode = 0;
		foreach my $rank_test("order","suborder","infraorder","series","superfamily","family","subfamily","tribe","subtribe","genus",
			"subgenus","species group","species","subspecies")
				{$current_rankcode++;
				if($child_rank eq $rank_test){$ncbi_nodes{$child_ID}{rank_code} = $current_rankcode}
				};

			$ncbi_nodes{$child_ID}{rank} = $child_rank;
			$ncbi_nodes{$child_ID}{parent} = $parent_ID;
			$ncbi_nodes{$parent_ID}{child_nodes} .= "\t" . $child_ID;
	

	
		$child_tax =~ s/[\(\)\,\[\]\'\#\&\/\:\.\-]/ /g;
		$child_tax =~ s/\s\s+/ /g;$child_tax =~ s/\s+$//;$child_tax =~ s/^\s+//;

		$child_tax =~ s/^Candidatus\s(\w+)$/$1/;

		$ncbi_nodes{$child_ID}{name} = $child_tax;


		$tax_table_line_count++;
		}else{
		print "cant parse line:$line\n";
		};

	};

close TAXTABLE;


print "taxon table has been read. 
	root taxon name:$root_taxon_name
 	node stored:$tax_table_line_count

";
}; # sub store_tax_heirarchy


#####################################################################################################
#
#
#
#####################################################################################################






sub traverse_nodes
{
my $current_node = $_[0];my $current_node_taxstring = $_[1];# $current_node_taxstring to be deprecated



 # for the current node, read the list of child nodes
my $child_nodes = $ncbi_nodes{$current_node}{child_nodes};$child_nodes =~ s/^\t//;
my @child_nodes_array = split(/\t/, $child_nodes);
$nodes_traversed++;


 # seems mostly discontinued:
if($current_node == $starting_node)
	{
	my $rank = $ncbi_nodes{$starting_node}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
#	my $current_node_taxstring = substr $current_node_taxstring, 0,1;
#	my $originalname = $ncbi_nodes{$starting_node}{name};$originalname =~ s/[\s\t]/_/g;
#	my $child_complete_lineage = $ncbi_nodes{$starting_node}{complete_lineage} . "$rank:$originalname ";
	};



foreach my $child(@child_nodes_array)
	{
	unless($current_node == $child) # one of the child nodes of the root node (1), is also 1 (NCBI system)
	{

	my $rank = $ncbi_nodes{$child}{rank};$rank =~ s/\s/_/g;$how_many{$rank}++;
	my $name_string = $ncbi_nodes{$child}{name};$name_string =~ s/\s+/_/g; ####################### sep2013 .. fixed mar2015
	my $rankcode = $ncbi_nodes{$child}{rank_code};$rank_codes{$name_string} = $rankcode;


	unless($ncbi_nodes{$current_node}{complete_lineage} =~ /\w/)#22oct2015:otherwise where shared taxon is same as basal node, wont be assigned 
		{$ncbi_nodes{$current_node}{complete_lineage}="root:$starting_name_ff "};

	my $child_complete_lineage = $ncbi_nodes{$current_node}{complete_lineage} . "$rank:$name_string ";#prev assigned $name


	# NOV 2018: user decides which ranks are constrained	
	my $constrain_this_node = 0; # print "\n\n";
	foreach my $user_constrain_rank ( @constrain_ranks )
		{
		if($rank eq $user_constrain_rank)
			{$constrain_this_node = 1};
	#	print "rank:$rank user:$user_constrain_rank constrain:$constrain_this_node\n";
		};
#	if($constrain_this_node == 1){
		$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;
#		};

# forgot about this, explains some imperfect species level trees
#	unless($rank eq "order") # 31 aug 2016: way to stop seqs only given order ident from being constrained.
#		{
#		$ncbi_nodes{$child}{complete_lineage} = $child_complete_lineage;
#		};


	#print "storing complete lineage for taxon:$name_string\n";	
	$complete_lineage_for_this_species{$name_string} = $child_complete_lineage;

	$ncbi_tax_number_for_this_species{$name_string} = $child;

#	print "complete lineage:$ncbi_nodes{$child}{complete_lineage}\n";

	my $originalname = $ncbi_nodes{$child}{name};$originalname =~ s/[\s\t]/_/g;
	$ncbi_taxnumber_for_taxname{$originalname} = $child;
#	print "ncbi_number:$child tax name:$originalname\n";

#	if ( $rank eq "species" ){};

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



